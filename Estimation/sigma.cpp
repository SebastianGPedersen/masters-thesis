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
Rcpp::NumericVector sigmas_cpp(Rcpp::NumericVector KdY, int lags) {
  //Assuming parzen kernel and exp kernel and that number of lags is even
  //The runtime of this algorithm is almost independent of the number of lags and calculates sigma at EVERY point
  // Kernels should be same size as dy with value kernels[i] = K((t_{i-1}-t_i)/bandwidth)
  // size of dy should be greater than lags+4
  
  int n = KdY.size();

  Rcpp::NumericVector sigmas(n); //The result vector
  Rcpp::NumericVector weights(n); //The weights for the loop
  Rcpp::NumericVector products(n);
  
  
  //Calculate parzen-weights with function above
  Rcpp::NumericVector parzen_weights = parzen_func(lags);
  
  //Correction weights used in final loop below calculated outside for faster correction
  Rcpp::NumericVector first_weights(4);
  Rcpp::NumericVector last_weights(4);
  Rcpp::NumericVector middle_weights(4);
  
  
  //First weights
  first_weights[0] = parzen_weights[lags];
  first_weights[1] = -4*parzen_weights[lags] + parzen_weights[lags-1];
  first_weights[2] = 6*parzen_weights[lags]-4*parzen_weights[lags-1] + parzen_weights[lags-2];
  first_weights[3] = -4*parzen_weights[lags]+6*parzen_weights[lags-1]-4*parzen_weights[lags-2] + parzen_weights[lags-3];
  
  
  //Middle weights - assuming L is even
  //The correct weights
  double lags_db = double(lags);
  double temp1 = 2.0*2*pow(1-(lags_db/2.0) / (lags_db+1),3);
  double temp2 = 2.0*2*pow(1-(lags_db/2.0-1.0) / (lags_db+1),3);
  double temp3 = 2.0*2*pow(1-(lags_db/2.0-2.0) / (lags_db+1),3);
  double temp4 = 2.0*2*pow(1-(lags_db/2.0-3.0) / (lags_db+1),3);
  
  //Extract the wrong side, add the right
  middle_weights[0] = parzen_weights[lags/2] - temp1; 
  middle_weights[1] = -4*parzen_weights[lags/2] + 4*temp1+parzen_weights[lags/2-1]-temp2;
  middle_weights[2] = 6*parzen_weights[lags/2]-6*temp1 -4*parzen_weights[lags/2-1]+4*temp2+parzen_weights[lags/2-2]-temp3;
  middle_weights[3] = -4*parzen_weights[lags/2]+4*temp1 + 6*parzen_weights[lags/2-1]-6*temp2 -4*parzen_weights[lags/2-2]+4*temp3+parzen_weights[lags/2-3]-temp4;
  
  
  //Last weights
  last_weights[3] = 2; 
  last_weights[2] = parzen_weights[1] -4*2;
  last_weights[1] = parzen_weights[2] -4*parzen_weights[1]+6*2;
  last_weights[0] = parzen_weights[3]  -4*parzen_weights[2]+6*parzen_weights[1]-4*2;
  
  
  // First calculate weights(L+1,L+2,L+3,L+4), and sigmas(L+1,L+2,L+3) - these are needed in the final loop
  //First sigma (divided to avoid if statement)
  for (int k=0; k<lags+1; k++) {
    weights[lags] += KdY[k]*parzen_weights[lags-k]; //includes lag = 0 (with twice weight, so will be substracted later)
    sigmas[lags] += KdY[k]*KdY[k];
  }
  
  long double temp;
  for (int lag =1; lag < lags+1;lag++) {
    temp = 0;
    for (int k = 0; k < lags+1-lag; k++){
      temp += KdY[k]*KdY[k+lag];
    }
    sigmas[lags] += temp*parzen_weights[lag];
  }

  
  //Next 3 weights and sigmas
  for (int i=1; i < 4; i++) {
    for (int k=0; k<lags+1; k++) {
      weights[lags+i] += KdY[k+i]*parzen_weights[lags-k]; //includes lag = 0 (with twice weight, so will be substracted later)
    }
    sigmas[lags+i] = sigmas[lags+i-1] + (weights[lags+i]-KdY[lags+i])*KdY[lags+i];
  }
  
  /*
  //-------------------Temporary-----------------
  Rcpp::NumericMatrix weights_temp(6,lags+6);
  
  //The first 4 weights and sigmas
  for (int i=0; i < 4; i++) {
    for (int k=0; k<lags+1; k++) {
      weights_temp(i,k+i) = KdY[k+i]*parzen_weights[lags-k]; //includes lag = 0 (with twice weight, so will be substracted later)
    }
  }    

  //The next weights
  for (int i = 4; i < 6; i++) {
    for (int j = 0; j < lags+4;j++) {
    weights_temp(i,i-4+j) = (4*weights_temp(i-1,i-4+j)-6*weights_temp(i-2,i-4+j)+4*weights_temp(i-3,i-4+j)-weights_temp(i-4,i-4+j));
    }
  }
  weights_temp(4,0) += first_weights[0]*KdY[0];
  weights_temp(4,1) += first_weights[1]*KdY[1];
  weights_temp(4,2) += first_weights[2]*KdY[2];
  weights_temp(4,3) += first_weights[3]*KdY[3];
  
  weights_temp(4,lags/2) += middle_weights[0]*KdY[lags/2];
  weights_temp(4,lags/2+1) += middle_weights[1]*KdY[lags/2+1];
  weights_temp(4,lags/2+2) += middle_weights[2]*KdY[lags/2+2];
  weights_temp(4,lags/2+3) += middle_weights[3]*KdY[lags/2+3];

  weights_temp(4,lags+1) += last_weights[0]*KdY[lags+1];
  weights_temp(4,lags+2) += last_weights[1]*KdY[lags+2];
  weights_temp(4,lags+3) += last_weights[2]*KdY[lags+3];
  weights_temp(4,lags+4) += last_weights[3]*KdY[lags+4];
  
  //Rcpp::Rcout << first_weights[0] << " || " << KdY[0] << " || " << parzen_weights[lags] << " || ";
  */
  
  //For third degree polynomial: f(x) = 4*f(x-k)-6*f(x-2k)+4*f(x-3k)-f(x-4k)
  //Using the equation above and then corrects the endpoints and middle (where the polynomium in parzen kernel changes)
  
  for (int i = lags+4; i < n; i++) {
    weights[i] = (4*weights[i-1]-6*weights[i-2]+4*weights[i-3]-weights[i-4]) + 
      first_weights[0]*KdY[i-lags-4]+first_weights[1]*KdY[i-lags-3] + first_weights[2]*KdY[i-lags-2] + first_weights[3]*KdY[i-lags-1] +
      middle_weights[0]*KdY[i-lags/2-4]+middle_weights[1]*KdY[i-lags/2-3] + middle_weights[2]*KdY[i-lags/2-2] + middle_weights[3]*KdY[i-lags/2-1] +
      last_weights[0]*KdY[i-3]+last_weights[1]*KdY[i-2] + last_weights[2]*KdY[i-1] + last_weights[3]*KdY[i];
      
    sigmas[i] = sigmas[i-1] + (weights[i]-KdY[i])*KdY[i];
  }
  
  
  return sigmas;
  //return weights;
  //return weights_temp;
}
