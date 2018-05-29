#include <Rcpp.h>
#include <math.h>


Rcpp::NumericVector parzen_func(int lags) {
  
  Rcpp::NumericVector parzen_weights(lags);
  double x;
  
  for (int i=0; i < lags+1; i++) {
    
    x = double(i)/(lags+1.0);
    if (x < 0.5) {
      parzen_weights[i] = 2.0*(1-6*pow(x,2) + 6*pow(x,3));
    }
    else {
      parzen_weights[i] = 2.0*2*pow(1-x,3);
    }
  }
  
  return parzen_weights;
}


// [[Rcpp::export]]
Rcpp::NumericVector sigmas_cpp(Rcpp::NumericVector dY, Rcpp::NumericVector kernels, int lags, double bandwidth) {
  //Assuming parzen kernel and exp kernel and that number of lags is even
  //The runtime of this algorithm is almost independent of the number of lags and calculates sigma at EVERY point
  // Kernels should be same size as dy with value kernels[i] = K((t_{i-1}-t_i)/bandwidth)
  // size of dy should be greater than lags+4
  
  int n = dY.size();

  Rcpp::NumericVector sigmas(n); //The result vector
  Rcpp::NumericVector weights(n); //The weights for the loop
  Rcpp::NumericVector products(n);
  
  //The calculations below are kernels later used in calculations of lag+1, lag+2, lag+3 and lag+4
  Rcpp::NumericVector sum_kernel(lags+4);
  Rcpp::NumericVector extra_scale(4);
  
  extra_scale[0] = kernels[lags];
  for (int i =1; i < 4; i++) {
    extra_scale[i] = extra_scale[i-1]*kernels[lags+i];
  }
  
  sum_kernel[lags+3] = 1/(kernels[lags]*kernels[lags+1]*kernels[lags+2]);
  for (int i = lags+3; i > 0; i--) {
    sum_kernel[i-1] = sum_kernel[i]*kernels[i-1];
  }
  
  
  //Calculate parzen-weights with function above
  Rcpp::NumericVector parzen_weights = parzen_func(lags);
  
  //Correction weights used in final loop below calculated outside for faster correction
  Rcpp::NumericVector first_weights(4);
  Rcpp::NumericVector last_weights(4);
  Rcpp::NumericVector middle_weights(4);
  
  
  //First weights
  first_weights[0] = parzen_weights[lags-1]; // because [i-4] is subtracted
  first_weights[1] = -4*parzen_weights[lags-1] + parzen_weights[lags-2];
  first_weights[2] = 6*parzen_weights[lags-1]-4*parzen_weights[lags-2] + parzen_weights[lags-3];
  first_weights[3] = -4*parzen_weights[lags-1]+6*parzen_weights[lags-2]-4*parzen_weights[lags-3] + parzen_weights[lags-4];
  
  
  //Middle weights - assuming L is even
  int temp1 = 2.0*2*pow(1-(lags/2.0) / (lags+1),3);
  int temp2 = 2.0*2*pow(1-(lags/2.0-1.0) / (lags+1),3);
  int temp3 = 2.0*2*pow(1-(lags/2.0-2.0) / (lags+1),3);
  int temp4 = 2.0*2*pow(1-(lags/2.0-3.0) / (lags+1),3);
  
  middle_weights[0] = parzen_weights[lags/2] - temp1; //Extract the wrong side, add the right
  middle_weights[1] = -4*parzen_weights[lags/2] + 4*temp1-parzen_weights[lags/2-1]+temp2;
  middle_weights[2] = 6*parzen_weights[lags/2]-6*temp1 -4*parzen_weights[lags/2-1]+4*temp2+parzen_weights[lags/2-2]-temp3;
  middle_weights[3] = -4*parzen_weights[lags/2]+4*temp1 + 6*parzen_weights[lags/2-1]-6*temp2 -4*parzen_weights[lags/2-2]+4*temp3+parzen_weights[lags/2-3]-temp4;
  
  
  //Last weights
  last_weights[3] = 2; 
  last_weights[2] = parzen_weights[1] -4*2;
  last_weights[1] = parzen_weights[2] -4*parzen_weights[1]+6*2;
  last_weights[0] = parzen_weights[3]  -4*parzen_weights[2]+6*parzen_weights[1]-4*2;
  
  
  // First calculate weights(L+1,L+2,L+3,L+4), these are needed in the final loop 
  for (int i=0; i < 4; i++) {
    for (int k=0; k<lags+1; k++) {
      weights[lags+i] += (sum_kernel[k+i]*extra_scale[i])*dY[k+i]*parzen_weights[lags-k]; //includes lag = 0 (with twice weight, so will be substracted later)
      //if ((i == 1) and (k == lags)) {Rcpp::Rcout << sum_kernel[k+i]*extra_scale[i] << " || " << dY[k+i]<< " || " << parzen_weights[lags-k] << " || ";}    
    }
  }
  
  //Temporary
  Rcpp::NumericMatrix weights_temp(6,lags+6);
  
  for (int i=0; i < 4; i++) {
    for (int k=0; k<lags+1; k++) {
      weights_temp(i,k+i) = (sum_kernel[k+i]*extra_scale[i])*dY[k+i]*parzen_weights[lags-k]; //includes lag = 0 (with twice weight, so will be substracted later)
    }
  }    
  
  for (int i = 4; i < 6; i++) {
    for (int j = 0; j < lags;j++) {
    weights_temp(i,i+j) = (4*weights_temp(i-4,i+j)-6*weights_temp(i-3,i+j)+4*weights_temp(i-2,i+j)-weights_temp(i-1,i+j));
  }
}
    /*
      first_weights[0]*dY[i-lags-4]+first_weights[1]*dY[i-lags-3] + first_weights[2]*dY[i-lags-2] + first_weights[3]*dY[i-lags-1] +
      middle_weights[0]*dY[i-lags/2-4]+middle_weights[1]*dY[i-lags/2-3] + middle_weights[2]*dY[i-lags/2-2] + middle_weights[3]*dY[i-lags/2-1] +
      last_weights[0]*dY[i-lags/2-4]+last_weights[1]*dY[i-lags/2-3] + last_weights[2]*dY[i-lags/2-2] + last_weights[3]*dY[i-lags/2-1];
    
    sigmas[i] = (sigmas[i-1] + 1/bandwidth * (weights[i]-1)*dY[i])*kernels[i]*kernels[i]; 
  }
  */
  /*
  //For third degree polynomial: f(x) = 4*f(x-k)-6*f(x-2k)+4*f(x-3k)-f(x-4k)
  //Using the equation above and then corrects the endpoints and middle (where the polynomium in parzen kernel changes)
  for (int i = lags+4; i < n; i++) {
    weights[i] = (4*weights[i-1]-6*weights[i-2]+4*weights[i-3]-weights[i-4]) + 
      first_weights[0]*dY[i-lags-4]+first_weights[1]*dY[i-lags-3] + first_weights[2]*dY[i-lags-2] + first_weights[3]*dY[i-lags-1] +
      middle_weights[0]*dY[i-lags/2-4]+middle_weights[1]*dY[i-lags/2-3] + middle_weights[2]*dY[i-lags/2-2] + middle_weights[3]*dY[i-lags/2-1] +
      last_weights[0]*dY[i-lags/2-4]+last_weights[1]*dY[i-lags/2-3] + last_weights[2]*dY[i-lags/2-2] + last_weights[3]*dY[i-lags/2-1];
      
    sigmas[i] = (sigmas[i-1] + 1/bandwidth * (weights[i]-1)*dY[i])*kernels[i]*kernels[i]; 
  }
  */
  //return sigmas;
  //return weights;
  return weights_temp;
}
