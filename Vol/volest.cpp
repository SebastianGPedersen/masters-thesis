#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

/*
* RCPP FOR SINGLE, NONMATRIX, PATH FOR REAL DATA - POSSIBLE ALSO FOR PERFORMANCE BENCHMARKING
*/

// [[Rcpp::export]]
NumericVector vol_est_preA(NumericVector x, double K){
  
  int N = x.length();
  
  NumericVector out(N-K+1);
  
  for(int i = 0; i < out.length(); i++){
    out[i] = (sum(x[Range((i+K/2.0),(i+K-1))]) - sum(x[Range(i,(i+K/2.0-1))]))/K;
  }
  
  return(out);
}

double vol_est_phik(double K){
  return(   (1.0+2.0*pow(K,-2))/12.0   );
}

// [[Rcpp::export]]
double vol_est_RVstar(NumericVector paX, double K, double omega2est, double theta){
  
  int N = paX.length() + K - 1;
  
  double phiK = vol_est_phik(K);
  
  double out = N/(N-K+2.0) * ( 1.0 / (phiK * K)) * sum(paX*paX) - omega2est/(phiK * theta*theta);
  
  return(out);
}

// [[Rcpp::export]]
double vol_est_BVstar(NumericVector paX, double K, double omega2est, double theta){
  
  int N = paX.length() + K - 1;
  
  double phiK = vol_est_phik(K);
  
  double out = N/(N-2.0*K+2.0) * (1.0 / (phiK * K)) * (3.141592653589793238462643383280/2.0) * sum(abs(paX[Range(0,paX.length()-K)]*paX[Range(K,paX.length())]))-omega2est/(theta*theta*phiK);
  
  return(out);
}

// [[Rcpp::export]]
double vol_est_BV(NumericVector paX){
  
  int N = paX.length();
  
  double out = N/(N-1.0) * (3.141592653589793238462643383280/2.0) * sum(   abs(   paX[Range(0, N-1)]*paX[Range(1,N)]   )   );
  
  return(out);
}
