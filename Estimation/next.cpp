#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

/*
 * RCPP FOR SINGLE, NONMATRIX, PATH FOR REAL DATA - POSSIBLE ALSO FOR PERFORMANCE BENCHMARKING
 */

NumericVector kern(NumericVector x){
  // handles only negative
  NumericVector neg = exp(-abs(x[x<=0]));
  return(neg);
}

double parzen(int x, double lag){
  // old school for loop :(
  double z = (x+0.0)/(lag+1.0);
  
  double out;
  if(z < 0.5){
    out = 1.0 - 6.0*pow(z,2.0) + 6.0 * pow(z, 3.0);
  }
  else{
    out = 2.0 * pow((1.0-z),3.0);
  }
  return(out);
}


// [[Rcpp::export]]
NumericVector mu_internal_cpp(NumericVector time, NumericVector dy, NumericVector t,
                              NumericVector rescale, NumericVector start, NumericVector end, double hd , double startmu){
  
  // init res
  int N = t.size();
  NumericVector mu(N);
  
  // Calc'eria (HUSK NUL INDEX!)
  mu[0] = startmu*rescale[0] + 1.0/(hd) * sum(   kern(  (time[Range(start[0]-1, end[0]-1)] - t[0])/(hd)  )*dy[Range(start[0], end[0])] );
  for(int j = 1; j < N; j++){
    mu[j] = mu[j-1]*rescale[j] + 1.0/(hd) * sum(   kern(  (time[Range(start[j]-1, end[j]-1)] - t[j])/(hd)  )*dy[Range(start[j], end[j])] );
  }
  return(mu);
}

// [[Rcpp::export]]
NumericVector sig_internal_cpp(NumericVector time, NumericVector dy, NumericVector t,
                              NumericVector rescale, NumericVector startr, NumericVector endr, double hv, int lag , double startsig){
  
  // init res
  int N = t.size();
  NumericVector sig(N);
  
  // zero indexify
  NumericVector start = startr-1; NumericVector end = endr-1;
  // pre-loop tricks
  int length  = max(end-start);
  NumericVector kdy(length+lag+1);
  
  NumericVector parz(lag+1);
  for(int i = 0; i < lag+1; i++){
    parz[i] = parzen(i, lag);
  }
  
  // ANFANGT!
  // should be updated every single round (i)
  kdy[Range(0, end[0]-start[0]+lag)] = kern((time[Range(start[0]-lag,end[0])] - t[0])/hv)*dy[Range(start[0]+1-lag, end[0]+1)];
  
  sig[0] = startsig*rescale[0] + (1.0/hv) * sum(    pow( kdy[Range(0+lag, end[0]-start[0]+lag)], 2.0)    );
  
  /*
  sig[0] = startsig*rescale[0] + (1.0/hv) * sum(    pow( 
    kern( (time[Range(start[0],end[0])] - t[0])/hv)*dy[Range(start[0]+1,end[0]+1)]
    , 2.0)    );
  */
  
  // lags
  for(int l = 1; l <= lag; l++){
    sig[0] = sig[0] + 2.0*(1.0/hv)*parz(l)*sum(kdy[Range(0+lag, end[0]-start[0]+lag)]*kdy[Range(0+lag-l, end[0]-start[0]+lag-l)]);
  }
  
  // rest
  for(int j = 1; j < N; j++){
    kdy[Range(0, end[j]-start[j]+lag)] = kern((time[Range(start[j]-lag,end[j])] - t[j])/hv)*dy[Range(start[j]+1-lag, end[j]+1)];
    
    sig[j] = sig[j-1]*rescale[j] + (1.0/hv) * sum(    pow( kdy[Range(0+lag, end[j]-start[j]+lag)], 2.0)    );
    // lags
    for(int l = 1; l <= lag; l++){
      sig[j] = sig[j] + 2.0*(1.0/hv)*parz(l)*sum(kdy[Range(0+lag, end[j]-start[j]+lag)]*kdy[Range(0+lag-l, end[j]-start[j]+lag-l)]);
    }
  }
  return(sig);
}


