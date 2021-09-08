#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace RcppEigen;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double cppnrm1(NumericVector x){
  double s = sum(abs(x));
  return(s);
}

// [[Rcpp::export]]
double cppnrm2(NumericVector x){
  const Eigen::Map<Eigen::VectorXd> v(as<Eigen::Map<Eigen::VectorXd> >(x));
  double y = v.norm();
  return(y);
}

// [[Rcpp::export]]
NumericVector cppmatvec(NumericMatrix m, NumericVector v){
  const Eigen::Map<Eigen::MatrixXd> tm(as<Eigen::Map<Eigen::MatrixXd> >(m));
  const Eigen::Map<Eigen::VectorXd> tv(as<Eigen::Map<Eigen::VectorXd> >(v));
  Eigen::VectorXd prod = tm*tv;
  SEXP s = Rcpp::wrap(prod);
  NumericVector w(s);
  return(w);
}

// [[Rcpp::export]]
double cppsoft(double x, double lam){
  double bx = abs(x) - lam;
  if(x<0) bx = -bx;
  return(bx);
}

// [[Rcpp::export]]
NumericVector cppl1prox(NumericVector x, double lam){
  NumericVector bx(x.length());
  int n = bx.length();
  for(int i = 0; i < n; i++){
    bx[i] = cppsoft(x[i], lam);
  }
  return(bx);
}

// [[Rcpp::export]]
NumericVector cppl2prox(NumericVector x, double lam){
  if(lam==0){
    return(x);
  }
  else{
    double nrm = cppnrm2(x);
    double t = 1-(lam/nrm);
    double s = 0;
    if(t > 0){
      s = t;
    }
    if(s==0){
      return(0 * x);
    }
    else{
      return(s * x);
    }
  }
}

