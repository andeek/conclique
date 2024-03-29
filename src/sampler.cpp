#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; using namespace arma;

//' Functions for sampling from conditional distributions.
//' 
//' @param data A list containing two elements, 
//'        sums and nums which contain the sum of the data in each neighborhood as well as the number of locations 
//'        in the neighborhood for each point in the conclique. 
//' @param params A list of parameter values, rho, kappa, and eta, that parameterize the Gaussian MRF.
//' @name sampler


//' @rdname sampler
//' @export
// [[Rcpp::export]]
arma::vec gaussian_single_param_sampler(List data, List params) {
  
  List sums = data["sums"];
  List nums = data["nums"];
  
  vec sums_vec = sums[0];
  vec nums_vec = nums[0];
  
  vec res(sums_vec.n_elem);
  
  RNGScope scope;
  
  double rho = params["rho"];
  double kappa = params["kappa"];
  double eta = params["eta"];
  
  vec mean_structure = kappa + eta * (sums_vec - nums_vec * kappa);
  res = rnorm(mean_structure.n_elem);
  res = res * rho + mean_structure;
  
  return(res);
}

//' @rdname sampler
//' @export
// [[Rcpp::export]]
arma::vec gaussian_two_param_sampler(List data, List params) {
  
  List sums = data["sums"];
  List nums = data["nums"];
  
  vec sums_vec_1 = sums[0];
  vec nums_vec_1 = nums[0];
  vec sums_vec_2 = sums[1];
  vec nums_vec_2 = nums[1];
  
  vec res(sums_vec_1.n_elem);
  
  RNGScope scope;
  
  double rho = params["rho"];
  double kappa = params["kappa"];
  double eta_1 = params["eta_1"];
  double eta_2 = params["eta_2"];
  
  vec mean_structure = kappa + eta_1 * (sums_vec_1 - nums_vec_1 * kappa) + eta_2 * (sums_vec_2 - nums_vec_2 * kappa);
  res = rnorm(mean_structure.n_elem);
  res = res * rho + mean_structure;
  
  return(res);
}


//' @rdname sampler
//' @export
// [[Rcpp::export]]
arma::vec binary_single_param_sampler(List data, List params) {
  List sums = data["sums"];
  List nums = data["nums"];
  
  vec sums_vec = sums[0];
  vec nums_vec = nums[0];
  
  vec res(sums_vec.n_elem);
  
  RNGScope scope;
  
  double kappa = params["kappa"];
  double eta = params["eta"];
  
  vec mean_structure = log(kappa) - log(1 - kappa) + eta * (sums_vec - nums_vec * kappa);
  
  int i;
  
  for(i = 0; i < res.n_elem; ++i) {
    res(i) = as<double>(rbinom(1, 1, exp(mean_structure(i))/(1 + exp(mean_structure(i)))));
  }
  return(res);
}

//' @rdname sampler
//' @export
// [[Rcpp::export]]
arma::vec binary_two_param_sampler(List data, List params) {
  
  List sums = data["sums"];
  List nums = data["nums"];
  
  vec sums_vec_1 = sums[0];
  vec nums_vec_1 = nums[0];
  vec sums_vec_2 = sums[1];
  vec nums_vec_2 = nums[1];
  
  vec res(sums_vec_1.n_elem);
  
  RNGScope scope;
  
  double kappa = params["kappa"];
  double eta_1 = params["eta_1"];
  double eta_2 = params["eta_2"];
  
  vec mean_structure = log(kappa) - log(1 - kappa) + eta_1 * (sums_vec_1 - nums_vec_1 * kappa) + eta_2 * (sums_vec_2 - nums_vec_2 * kappa);
  
  int i;
  
  for(i = 0; i < res.n_elem; ++i) {
    res(i) = as<double>(rbinom(1, 1, exp(mean_structure(i))/(1 + exp(mean_structure(i)))));
  }
  return(res);
}

//' @rdname sampler
//' @export
// [[Rcpp::export]]
arma::vec binary_two_param_reg_sampler(List data, List params) {
  
  List sums = data["sums"];
  List nums = data["nums"];
  vec u = data["u"];
  
  vec sums_vec_1 = sums[0];
  vec nums_vec_1 = nums[0];
  vec sums_vec_2 = sums[1];
  vec nums_vec_2 = nums[1];
  
  vec res(sums_vec_1.n_elem);
  
  RNGScope scope;
  
  double beta_0 = params["beta_0"];
  double beta_1 = params["beta_1"];
  double eta_1 = params["eta_1"];
  double eta_2 = params["eta_2"];
  
  vec reg = beta_0 + beta_1 * u;
  vec kappa = exp(reg) / (1 + exp(reg));
  
  vec mean_structure = reg + eta_1 * (sums_vec_1 - nums_vec_1 % kappa) + eta_2 * (sums_vec_2 - nums_vec_2 % kappa);
  
  int i;
  
  for(i = 0; i < res.n_elem; ++i) {
    res(i) = as<double>(rbinom(1, 1, exp(mean_structure(i))/(1 + exp(mean_structure(i)))));
  }
  return(res);
}

