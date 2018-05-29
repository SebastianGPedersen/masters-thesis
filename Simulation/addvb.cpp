#include <Rcpp.h>
#include <math.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix vol_add_cpp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix vol, Rcpp::NumericVector sigma_t) {
  
  int steps = X.ncol()-1;
  int paths = X.nrow();
  
  Rcpp::NumericMatrix sigma_add(paths,steps+1);
  
  for (int i=0;i < paths;i++) {
    for (int j = 1; j < steps + 1; j++) {
      sigma_add(i,j) = sigma_add(i,j-1) + ((X(i,j)-X(i,j-1)) / pow(vol(i,j-1),0.5))*sigma_t(j-1);
    }
  }

  return sigma_add;
}