#include <RcppArmadillo.h>
#include "cdf.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; using namespace arma;

typedef arma::vec (*cdfPtr)(List, List);

//' Get spatial residuals from data given a neighborhood structure and MRF model.
//' 
//' @param data A vector containing data values for each location in the lattice.
//' @param neighbors A matrix N*N by (max # neighbors) + 1, where the first column is the location id of each location in the lattice. This could be the result from get_neighbors().
//' @param conditional_cdf A function that has two inputs: 
//'        \itemize{
//'          \item{data, and}
//'          \item{params.}
//'        }
//'        The data is a list containing data, sums, and nums, which contains the sampled values for each location, the neighboring locations of each location,
//'        and the number of neighbors each location has, respectively. params is a list of parameter values. This function returns the inverse cdf 
//'        at a value between 0 and 1 from the conditional distribution
//' @param params A list of parameters to be passed to the conditional_density function 
//' @export
// [[Rcpp::export]]
arma::vec spatial_residuals(arma::vec data, List neighbors, std::string conditional_cdf, List params) {
  bool r_func = true;
  Environment global = Environment::global_env();
  
  std::map<std::string, cdfPtr> cdf_map;
  cdf_map["gaussian_single_param"] = gaussian_single_param_cdf;
  
  // defined samplers in the package
  std::map<std::string, cdfPtr>::iterator it;
  it = cdf_map.find(conditional_cdf);
  if (it != cdf_map.end())
    r_func = false;
    
  int N = neighbors.length();
  List sums(N);
  List nums(N);
  int n, m;
  
  for(n = 0; n < N; ++n) { // neighborhood structures
    mat neigh = neighbors[n];
    int q = neigh.n_cols - 1;
    uvec cols = regspace<uvec>(1,  1,  q);
    uvec idx = conv_to<uvec>::from(vectorise(neigh));
    
    mat dat = data.elem(idx - 1);
    dat.reshape(data.n_rows, q);
    
    vec sums_inner(dat.n_rows);
    vec nums_inner(dat.n_rows);
    for(m = 0; m < dat.n_rows; ++m) { // rowSums
      sums_inner(m) = sum(dat.row(m)); // need NA handling for non regular neighborhoods
      nums_inner(m) = dat.row(m).n_elem;
    } // end m
    sums[n] = sums_inner;
    nums[n] = nums_inner;
  } // end n
  
  List sums_nums;
  sums_nums["data"] = data;
  sums_nums["sums"] = sums;
  sums_nums["nums"] = nums;

  if(r_func) {
    Function cdf = global[conditional_cdf];
    NumericVector resid = cdf(sums_nums, params);
    vec res(resid.begin(), resid.length());
    return(res);
  } else {
    vec res = cdf_map[conditional_cdf](sums_nums, params);
    return(res);
  }

}
