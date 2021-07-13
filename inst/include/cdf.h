#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; using namespace arma;


#ifndef CDF_H
#define CDF_H
arma::vec gaussian_single_param_cdf(List data, List params); 
arma::vec binary_single_param_cdf(List data, List params); 
#endif