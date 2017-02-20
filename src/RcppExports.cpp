// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// gaussian_single_param_cdf
arma::vec gaussian_single_param_cdf(List data, List params);
RcppExport SEXP conclique_gaussian_single_param_cdf(SEXP dataSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussian_single_param_cdf(data, params));
    return rcpp_result_gen;
END_RCPP
}
// run_conclique_gibbs
arma::mat run_conclique_gibbs(List conclique_cover, List neighbors, arma::mat inits, std::string conditional_sampler, List params, int n_iter);
RcppExport SEXP conclique_run_conclique_gibbs(SEXP conclique_coverSEXP, SEXP neighborsSEXP, SEXP initsSEXP, SEXP conditional_samplerSEXP, SEXP paramsSEXP, SEXP n_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type conclique_cover(conclique_coverSEXP);
    Rcpp::traits::input_parameter< List >::type neighbors(neighborsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inits(initsSEXP);
    Rcpp::traits::input_parameter< std::string >::type conditional_sampler(conditional_samplerSEXP);
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(run_conclique_gibbs(conclique_cover, neighbors, inits, conditional_sampler, params, n_iter));
    return rcpp_result_gen;
END_RCPP
}
// run_sequential_gibbs
arma::mat run_sequential_gibbs(List neighbors, arma::mat inits, std::string conditional_sampler, List params, int n_iter);
RcppExport SEXP conclique_run_sequential_gibbs(SEXP neighborsSEXP, SEXP initsSEXP, SEXP conditional_samplerSEXP, SEXP paramsSEXP, SEXP n_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighbors(neighborsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inits(initsSEXP);
    Rcpp::traits::input_parameter< std::string >::type conditional_sampler(conditional_samplerSEXP);
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(run_sequential_gibbs(neighbors, inits, conditional_sampler, params, n_iter));
    return rcpp_result_gen;
END_RCPP
}
// gaussian_single_param_sampler
arma::vec gaussian_single_param_sampler(List data, List params);
RcppExport SEXP conclique_gaussian_single_param_sampler(SEXP dataSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussian_single_param_sampler(data, params));
    return rcpp_result_gen;
END_RCPP
}
// binary_single_param_sampler
arma::vec binary_single_param_sampler(List data, List params);
RcppExport SEXP conclique_binary_single_param_sampler(SEXP dataSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(binary_single_param_sampler(data, params));
    return rcpp_result_gen;
END_RCPP
}
// binary_two_param_sampler
arma::vec binary_two_param_sampler(List data, List params);
RcppExport SEXP conclique_binary_two_param_sampler(SEXP dataSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(binary_two_param_sampler(data, params));
    return rcpp_result_gen;
END_RCPP
}
// spatial_residuals
arma::vec spatial_residuals(arma::vec data, List neighbors, std::string conditional_cdf, List params);
RcppExport SEXP conclique_spatial_residuals(SEXP dataSEXP, SEXP neighborsSEXP, SEXP conditional_cdfSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type neighbors(neighborsSEXP);
    Rcpp::traits::input_parameter< std::string >::type conditional_cdf(conditional_cdfSEXP);
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(spatial_residuals(data, neighbors, conditional_cdf, params));
    return rcpp_result_gen;
END_RCPP
}
