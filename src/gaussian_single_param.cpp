#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; using namespace arma;

//' Functions for using conclique based GOF test for Gaussian MRF with a single dependence parameter.
//' 
//' @param data A list containing two elements, 
//'        sums and nums which contain the sum of the data in each neighborhood as well as the number of locations 
//'        in the neighborhood for each point in the conclique. 
//' @param params A list of parameter values, rho, kappa, and eta, that parameterize the Gaussian MRF.
//' @name gaussian_single_param


//' @rdname gaussian_single_param
// [[Rcpp::export]]
arma::vec gaussian_single_param_sampler(List data, List params) {
  RNGScope scope;
  
  double rho = params["rho"];
  double kappa = params["kappa"];
  double eta = params["eta"];
  
  vec sums = data["sums"];
  vec nums = data["nums"];
  
  sums = sums[1];
  nums = nums[1];
  
  vec mean_structure = kappa + eta * (sums - nums * kappa);
  vec res(rnorm(mean_structure.n_elem));
  res = res * rho + mean_structure;
  
  return(res);
}

//' @rdname gaussian_single_param
// [[Rcpp::export]]
arma::mat gaussian_single_param_cdf(List data, List params) {
  RNGScope scope;
  
  double rho = params["rho"];
  double kappa = params["kappa"];
  double eta = params["eta"];
  
  vec sums = data["sums"];
  vec nums = data["nums"];
  mat datum = data["data"];
  
  vec mean_structure = kappa + eta * (sums - nums * kappa);
  NumericMatrix res(datum.n_rows, datum.n_cols);
  for(int i = 0; i < mean_structure.n_elem; ++i) {
    colvec column = datum.col(i);
    NumericVector prob = pnorm(as<NumericVector>(wrap(column)), mean_structure(i), rho);
    res(_, i) = prob;
  }
  mat res_mat(res.begin(), res.nrow(), res.ncol());
  return(res_mat);
}



