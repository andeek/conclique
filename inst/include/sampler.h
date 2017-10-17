#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; using namespace arma;


#ifndef SAMPLER_H
#define SAMPLER_H
arma::vec gaussian_single_param_sampler(List data, List params); 
arma::vec binary_single_param_sampler(List data, List params); 
arma::vec binary_two_param_sampler(List data, List params); 
arma::vec binary_two_param_reg_sampler(List data, List params); 
#endif