#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; using namespace arma;

//' Functions for using conclique based GOF test.
//' 
//' @param data A list containing two elements, 
//'        sums and nums which contain the sum of the data in each neighborhood as well as the number of locations 
//'        in the neighborhood for each point in the conclique. 
//' @param params A list of parameter values, rho, kappa, and eta, that parameterize the Gaussian MRF.
//' @name cdf

//' @rdname cdf
//' @export
// [[Rcpp::export]]
arma::vec gaussian_single_param_cdf(List data, List params) {
  
  NumericVector datum = data["data"];
  vec res(datum.length());
  
  RNGScope scope;
  
  double rho = params["rho"];
  double kappa = params["kappa"];
  double eta = params["eta"];
  
  List sums = data["sums"];
  List nums = data["nums"];
  
  vec sums_vec = sums[0];
  vec nums_vec = nums[0];
  
  vec mean_structure = kappa + eta * (sums_vec - nums_vec * kappa);
  
  for(int i = 0; i < mean_structure.n_elem; ++i) {
    NumericVector prob = pnorm(datum, mean_structure(i), rho);
    res(i) = prob(i);
  }
  return(res);
}

//' Functions for using conclique based GOF test.
//' 
//' @param data A list containing two elements, 
//'        sums and nums which contain the sum of the data in each neighborhood as well as the number of locations 
//'        in the neighborhood for each point in the conclique. 
//' @param params A list of parameter values, rho, kappa, and eta, that parameterize the Gaussian MRF.
//' @name cdf

//' @rdname cdf
//' @export
// [[Rcpp::export]]
arma::vec binary_single_param_cdf(List data, List params) {
  
  NumericVector datum = data["data"];
  vec res(datum.length());
  
  RNGScope scope;
  
  double kappa = params["kappa"];
  double eta = params["eta"];
  
  List sums = data["sums"];
  List nums = data["nums"];
  
  vec sums_vec = sums[0];
  vec nums_vec = nums[0];
  
  vec mean_structure = log(kappa) - log(1 - kappa) + eta * (sums_vec - nums_vec * kappa);
  
  for(int i = 0; i < mean_structure.n_elem; ++i) {
    NumericVector prob = pbinom(datum, 1, exp(mean_structure(i))/(1 + exp(mean_structure(i))));
    res(i) = prob(i);
  }
  return(res);
}



