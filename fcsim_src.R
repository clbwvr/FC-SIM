# Required libraries
library(igraph)
library(mgcv)
library(Matrix)
library(matrixcalc)
library(splines)
library(glmnet)
library(fields)

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
  col5 <- colorRampPalette(c('#377EB8', 'white', 'red'))  #create color ramp starting from blue to red
  color_levels=100 #the number of colors to use
  max_absolute_value=max(abs(c(min(m), max(m)))) #what is the maximum absolute value of raster?
  color_sequence=seq(-max_absolute_value,max_absolute_value,length.out=color_levels+1)
  image.plot(m, col=col5(n=color_levels), breaks=color_sequence,main=main,xaxt='n',yaxt='n')
}

# Negative loglikelihood function
l = function(fxb,y,family){
  if(family=="gaussian"){
    return(mean((y-fxb)^2))
  }
  else if(family=="binomial"){
    pi = 1/(exp(-fxb)+1)
    pi[pi==1] = 1-1e-16
    pi[pi==0] = 1e-16
    return(-mean(y * log(pi) + (1-y)*log(1-pi)))
  }
}

# Gradient of negative loglikelihood function
dl = function(fxb,dfxb,X,y,family){
  if(family=="gaussian"){
    n = length(y)
    s = t(X) %*%  hadamard.prod(fxb - y, dfxb)
    return((2/n) * s)
  }
  else if(family=="binomial"){
    n = length(y)
    exb = exp(fxb)
    s=0
    for(i in 1:n) s = s + (y[i] * dfxb[i] * X[i,])/ ((exb[i]+1)) - (1-y[i]) * (exb[i] * dfxb[i] * X[i,]) / (exb[i] + 1)
    return(-s/n)
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
Proj = function(x, fixed = 1){
  x = x / nrm2(x)
  if(x[fixed] != 0) x = sign(x[fixed]) * x
  return(x)
}

# ADMM Updates
c_update = function(b, w, rho,fixed=1){
  Proj(b+w, fixed)
}

z_update=function(y,X,b,z,u,group,Glist,lambda2,rho){
  aj = lapply(unique(group), function(g) (b[Glist[[g]]] + u[group==g]))
  zjk = lapply(aj, l2prox, lam=lambda2/rho)
  zk = do.call(base::c, zjk)
  zk[is.na(zk)] = 0
  return(zk)
}

u_update=function(y,X,beta,z,u,group,Glist,rho){
  ujk = lapply(unique(group), function(g) u[group==g] + beta[Glist[[g]]] - z[group==g] )
  uk = do.call(base::c, ujk)
  uk[is.na(uk)] = 0
  return(uk)
}

b_update = function(y,X,fhat,fdhat,b,c,z,u,gs,w,lambda1,rho,opt="PG",family="gaussian",maxiter=50,tol=1e-3){
  if(!(opt %in% c("PG","BFGS"))) stop("Optimization method not supported")
  n=nrow(X)
  p=ncol(X)
  ubar = aggregate(u~gs$betaidx,FUN=mean)[,2]
  zbar = aggregate(z~gs$betaidx,FUN=mean)[,2]
  obj = function(fhat,beta,family=family,lambda1,rho=rho){
    xb = X%*%beta
    fxb = fhat(xb)
    return(l(fxb, y, family=family) + lambda1 * nrm1(beta) + (rho/2)*(nrm2(beta-zbar+ubar)^2) +  (rho/2)*(nrm2(beta-c+w)^2))
  }
  grad = function(fhat,fdhat,beta,lambda1,family=family,rho=rho){
    xb = X%*%beta
    fxb = fhat(xb)
    dfxb = fdhat(xb)
    return(dl(fxb=fxb, dfxb = dfxb, X=X, y=y, family=family) + lambda1*sign(beta) + rho*((beta-zbar+ubar) + (beta-c+w)) )
  }
  fobj = function(fhat,beta, family=family){
    xb = X%*%beta
    fxb = fhat(xb)
    return(l(fxb, y, family=family))
  }
  fgrad = function(fhat,fdhat,beta,family=family){
    xb = X%*%beta
    fxb = fhat(xb)
    dfxb = fdhat(xb)
    return(dl(fxb=fxb, dfxb = dfxb, X=X, y=y, family=family) )
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
        b = l1prox(step, lambda1 * gam)
        d = f1(b) - obj_val - sum(g_prev*(b-b_prev)) - (L/2) * nrm2(b-b_prev)^2
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
  ubar = aggregate(u~gs$betaidx,FUN=mean)[,2]
  zbar = aggregate(z~gs$betaidx,FUN=mean)[,2]
  b1n = sum(abs(b))
  bj2n = sum(sapply(unique(gs$group), function(g) (sqrt(sum(b[gs$Glist[[g]]]^2)))))
  aj = lapply(unique(gs$group), function(g) (b[gs$Glist[[g]]] + u[gs$group==g]))
  zjk = lapply(aj, l2prox, lam=lambda2/rho)
  loss =  l(X%*%b,y,family) + lambda1 * b1n + lambda2 * bj2n
  losses =c(loss)
  
  # f and b update loop
  for(it in 1:niter){
    
    # update f
    ti = X%*%b
    if(linear){
      fhat = function(u) u 
      fdhat = function(u) 0*u + 1
    } 
    else{
      psp = pspline(y, ti, family=family)
      fhat = function(u) c(predict(psp,newdata = data.frame(x=u)))
      fdhat = function(u) approxfun(c(ti), c(gratia::fderiv(psp,newdata=data.frame(x=ti))$derivatives$`s(x)`$deriv), rule = 2)(u)
    }
    
    # update b
    if(linear) rho=1e-12
    oldloss = 1e12
    for(i in 1:bniter){
      # ADMM updates
      res = b_update(y,X,fhat,fdhat,b,c,z,u,gs,w,lambda1=lambda1,rho=rho,family=family,maxiter=1,tol=1e-4,opt=opt)
      b = res[['b']]
      c = c_update(b,w,rho,fixed)
      w = w+b-c
      z = z_update(y,X,b,z,u,gs$group,gs$Glist,lambda2,rho)
      u = u_update(u,X,b,z,u,gs$group,gs$Glist,rho)
      b1n = sum(abs(b))
      bj2n = sum(sapply(unique(gs$group), function(g) (sqrt(sum(b[gs$Glist[[g]]]^2)))))
      loss =  l(fhat(X%*%b),y,family) + lambda1 * b1n + lambda2 * bj2n
      if(abs(loss-oldloss)/oldloss < tol) break
      oldloss = loss
    }
    b1n = sum(abs(b))
    bj2n = sum(sapply(unique(gs$group), function(g) (sqrt(sum(b[gs$Glist[[g]]]^2)))))
    loss =  l(fhat(X%*%b),y,family) + lambda1 * b1n + lambda2 * bj2n
    if(abs(loss - losses[it])/losses[it] < tol) break
    losses = c(losses, loss)
    if(linear) break
  }
  b=unname(b)
  return(list("fhat"=fhat, "bhat"=b, "t"=X%*%b, "pred"=fhat(ti), "losses"=losses))
}

