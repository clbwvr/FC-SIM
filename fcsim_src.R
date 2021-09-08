# Required libraries
library(igraph)
library(mgcv)
library(Matrix)
library(matrixcalc)
library(splines)
library(glmnet)
library(fields)
library(gratia)
library(data.table)
library(Rcpp)

# Source cpp functions
sourceCpp("cppsrc.cpp")

# Helper functions
nrm1 = function(x) sum(abs(x))

nrm2 = function(x) sqrt(sum(x^2))

l1prox = function(x,t){
  soft = function(u) sign(u)*max(abs(u)-t,0)
  return(sapply(x,soft))
}

l2prox = function(x,lam){
  nrm = nrm2(x)
  s = max(0, 1 - lam/nrm)
  return(s*x)
}

l12prox = function(x,lam1,lam2){
  l1prox(l2prox(x,lam2),lam1)
}

ilogit = function(x,a=0,b=1) 1/(1+exp(-(a+b*x)))

logseq = function(from,to,length.out) exp(seq(log(from), log(to), length.out = length.out))

inner = function(A,B){
  a=A[lower.tri(A)]
  b=B[lower.tri(B)]
  return(sum(a*b))
}

lower.to.sym = function(x){
  x[upper.tri(x)] = t(x)[upper.tri(x)]
  x
}

matshow = function(m,main=""){
  p = ncol(m)
  m = apply(m, 2, rev)
  col5 <- colorRampPalette(c('#377EB8', 'white', 'red'))
  max_absolute_value=max(abs(c(min(m), max(m)))) 
  color_sequence=seq(-max_absolute_value,max_absolute_value,length.out=101)
  image.plot(m, col=col5(n=100), breaks=color_sequence,main=main,xaxt='n',yaxt='n')
}

# Negative loglikelihood function
l = function(fxb,y,family){
  if(family=="gaussian"){
    return(mean((fxb - y)^2))
  }
  else if(family=="binomial"){
    pi = 1/(exp(-fxb)+1)
    pi[pi==1] = 1-1e-16
    pi[pi==0] = 1e-16
    return(-mean(y * log(pi) + (1-y)*log(1-pi)))
  }
}

# Gradient of negative loglikelihood function
dl = function(fxb,dfxb,tX,y,family){
  if(family=="gaussian"){
    n = length(y)
    s = cppmatvec(tX, ((fxb - y) * dfxb))
    return((2/n) * s)
  }
  else if(family=="binomial"){
    n = length(y)
    p = 1/(exp(-fxb)+1)
    s = cppmatvec(tX, ((p - y) * dfxb))
    return((2/n) * s)
  }
}

# Get groups for node & edge group penalty
# p: # of nodes
get.groups = function(p){
  midx = matrix(NA,p,p)
  midx[lower.tri(midx)] = 1:(p*(p-1)/2)
  group = c()
  betaidx = c()
  for(g in 1:p){
    for(j in 1:p){
      if(j < g){
        group = c(group, g)
        betaidx = c(betaidx, midx[g,j])
      }
      if(j > g){
        group = c(group, g)
        betaidx = c(betaidx, midx[j,g])
      }
    }
  }
  Glist = list()
  for(i in 1:max(group)){
    Glist[[i]] = betaidx[group==i]
  }
  list(Glist=Glist, group=group, betaidx=betaidx)
}

# P-spline fit
pspline = function(y, x, bs="ps", k=7, knots=NULL, family=family) {
  gamfit = gam(y~s(x, k = k, bs=bs),family=family)
  return(gamfit)
}

# Projection onto constraint set of c
proj = function(x, fixed = 1){
  x = x / nrm2(x)
  if(x[fixed] != 0) x = sign(x[fixed]) * x
  return(x)
}

# ADMM Updates
c_update = function(b, w, rho,fixed=1){
  proj(b+w, fixed)
}

z_update=function(y,X,b,z,u,group,Glist,lambda2,rho){
  aj = lapply(unique(group), function(g) (b[Glist[[g]]] + u[group==g]))
  zjk = lapply(aj, l2prox, lam=lambda2/rho)
  zk = unlist(zjk)
  zk[is.na(zk)] = 0
  return(zk)
}

u_update=function(y,X,beta,z,u,group,Glist,rho){
  ujk = lapply(unique(group), function(g) u[group==g] + beta[Glist[[g]]] - z[group==g] )
  uk = unlist(ujk)
  uk[is.na(uk)] = 0
  return(uk)
}

b_update = function(y,X,tX,fhat,fdhat,b,c,z,u,gs,w,lambda1,rho,opt="PG",family="gaussian",maxiter=50,tol=1e-3){
  if(!(opt %in% c("PG","BFGS"))) stop("Optimization method not supported")
  n=nrow(X)
  p=ncol(X)
  datt = data.table(u=u, z=z, g=gs$betaidx)
  agg = datt[, .(u = mean(u), z=mean(z)), by = g]
  ubar = agg$u
  zbar = agg$z
  obj = function(fhat,beta,family=family,lambda1,rho=rho){
    xb = cppmatvec(X, beta)
    fxb = fhat(xb)
    return(l(fxb, y, family=family) + lambda1 * nrm1(beta) + (rho/2)*(nrm2(beta-zbar+ubar)^2) +  (rho/2)*(nrm2(beta-c+w)^2))
  }
  grad = function(fhat,fdhat,beta,lambda1,family=family,rho=rho){
    xb = cppmatvec(X, beta)
    fxb = fhat(xb)
    dfxb = fdhat(xb)
    return(dl(fxb=fxb, dfxb = dfxb, tX=tX, y=y, family=family) + lambda1*sign(beta) + rho*((beta-zbar+ubar) + (beta-c+w)) )
  }
  fobj = function(fhat,beta, family=family){
    xb = cppmatvec(X, beta)
    fxb = fhat(xb)
    return(l(fxb, y, family=family))
  }
  fgrad = function(fhat,fdhat,beta,family=family){
    xb = cppmatvec(X, beta)
    fxb = fhat(xb)
    dfxb = fdhat(xb)
    return(dl(fxb=fxb, dfxb = dfxb, tX=tX, y=y, family=family) )
  }
  
  f1 = function(u) obj(fhat=fhat,beta=u,lambda1=lambda1,family=family,rho=rho)
  f2 = function(u) grad(fhat=fhat,fdhat=fdhat,beta=u,lambda1=lambda1,family=family,rho=rho)
  f3 = function(u) fobj(fhat=fhat,beta=u,family=family)
  f4 = function(u) fgrad(fhat=fhat,fdhat=fdhat,beta=u,family=family)
  
  losses=c()
  # Minimize via proximal gradient descent
  if(opt=="PG"){
    L=1
    for(it in 1:maxiter){
      obj_val = f1(b)
      b_prev = b
      g_prev = f4(b_prev)
      d = 1
      # backtracking line search
      while(d>1e-3){
        gam = 1/L
        step = b_prev - gam * g_prev  
        b = cppl1prox(step, lambda1 * gam)
        d = f1(b) - obj_val - sum(g_prev*(b-b_prev)) - (L/2) * cppnrm2(b-b_prev)^2
        L = L*1.2
      }
      L = L/1.2   
      losses = c(losses, obj_val)
      #stop if change in obj_val is small enough
      if(abs((obj_val - f1(b))/obj_val) < tol) break
    }
  }
  # Minimize via BFGS
  if(opt=="BFGS"){
    obj_val1 = f1(b)
    o1 = optim(b, f1, f2, method="BFGS")
    b = o1$par
    obj_val2 = f1(b)
    obj_vals = c(obj_val1, obj_val2)
  }
  
  return(list("b"=b,"losses"=losses))
}

# Model estimation
fcsim = function(X,y,gs,lambda=NULL,rho=.1,tol=1e-3,beta=NULL,niter=30,bniter=5,delta=0.95,linear=FALSE,family="gaussian",opt="PG"){
  # Set tuning parameter values
  n = length(y)
  if(is.null(lambda)){
    sx = as.matrix(scale(X))
    sy = as.vector(scale(y))
    lambda = .5 * max(abs(colSums(sx*sy)))/n
  }
  pl = length(gs$Glist[[1]])
  lambda1 = delta * lambda
  lambda2 = (1-delta) * sqrt(pl) * lambda
  
  # Initialize parameters
  b = coef(cv.glmnet(X,y,family=family,intercept=FALSE))[-1]
  fixed = which.max(abs(cor(X,y)))
  if(nrm2(b)==0) b[fixed] = 1
  c = w = rep(0, length(b))
  z = u = 0*b[gs$betaidx]
  datt = data.table(u=u, z=z, g=gs$betaidx)
  agg = datt[, .(u = mean(u), z=mean(z)), by = g]
  ubar = agg$u
  zbar = agg$z
  b1n = cppnrm1(b)
  bj2n = sum(sapply(unique(gs$group), function(g) cppnrm2(b[gs$Glist[[g]]]^2)))
  aj = lapply(unique(gs$group), function(g) (b[gs$Glist[[g]]] + u[gs$group==g]))
  zjk = lapply(aj, l2prox, lam=lambda2/rho)
  loss = l(cppmatvec(X, b),y,family) + lambda1 * b1n + lambda2 * bj2n
  losses = c(loss)
  tX = t(X)
  
  # f and b update loop
  for(it in 1:niter){
    
    # update f
    ti = cppmatvec(X, b)
    if(linear){
      fhat = function(u) u 
      fdhat = function(u) 0*u + 1
    } 
    else{
      psp = pspline(y, ti, family=family)
      #fhat = function(u) c(predict(psp,newdata = data.frame(x=u)))
      
      # speed up computation with approxfun
      apsp = approxfun(ti, predict(psp),rule = 2)
      fhat = function(u){apsp(u);}
      fdhat = function(u) approxfun(c(ti), c(fderiv(psp,newdata=data.frame(x=ti))$derivatives$`s(x)`$deriv), rule = 2)(u)
    }
    
    # update b
    if(linear) rho=1e-12
    oldloss = 1e12
    for(i in 1:bniter){
      # ADMM updates
      res = b_update(y,X,tX,fhat,fdhat,b,c,z,u,gs,w,lambda1=lambda1,rho=rho,family=family,maxiter=1,tol=1e-4,opt=opt)
      b = res[['b']]
      c = c_update(b,w,rho,fixed)
      w = w+b-c
      z = z_update(y,X,b,z,u,gs$group,gs$Glist,lambda2,rho)
      u = u_update(u,X,b,z,u,gs$group,gs$Glist,rho)
      b1n = cppnrm1(b)
      bj2n = sum(sapply(unique(gs$group), function(g) cppnrm2(b[gs$Glist[[g]]])))
      loss =  l(fhat(cppmatvec(X, b)),y,family) + lambda1 * b1n + lambda2 * bj2n
      if(abs(loss-oldloss)/oldloss < tol) break
      oldloss = loss
    }
    b1n = cppnrm1(b)
    bj2n = sum(sapply(unique(gs$group), function(g) cppnrm2(b[gs$Glist[[g]]])))
    loss =  l(fhat(cppmatvec(X, b)),y,family) + lambda1 * b1n + lambda2 * bj2n
    if(abs(loss - losses[it])/losses[it] < tol) break
    losses = c(losses, loss)
    if(linear) break
  }
  b=unname(b)
  return(list("fhat"=fhat, "bhat"=b, "t"=cppmatvec(X, b), "pred"=fhat(ti), "losses"=losses))
}
