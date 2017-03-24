//include <RcppArmadillo.h>
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
//'        There are three built in cdfs:
//'         \itemize{
//'           \item{"gaussian_single_param" - a Gaussian cdf with a single dependence parameter,}
//'           \item{"binary_single_param" - a binary cdf with a single dependence parameter, and}
//'           \item{"binary_two_param" - a binary cdf with two dependence parameters.}
//'         }
//'        If the user chooses to write their own cdf in R, they must pass the name of the cdf function that is available in the global environment as this parameter.
//'        The input "data" is a list containing at least two elements, 
//'        sums and nums which contain the sum of the data in each neighborhood as well as the number of locations 
//'        in the neighborhood for each point in the conclique. In addition, the data can contain two additional elements, u and v, 
//'        which are vectors that contain the horizontal and vertical location of each point in space.
//'        The input "params" is a list of parameter values. This function returns the inverse cdf at a value between 0 and 1 from the conditional distribution
//' @param params A list of parameters to be passed to the conditional_density function 
//' @param ncols A integer of the number of columns in the original grid. Only necessary if dealing with "u" and "v" in cdf function.
//' @export
// [[Rcpp::export]]
arma::vec spatial_residuals(arma::vec data, List neighbors, std::string conditional_cdf, List params, int ncols = 0) {
  int N = neighbors.length();
  List sums(N);
  List nums(N);
  int n, m;
  vec loc_u(N), loc_v(N);
  
  if(ncols > 0) {
    // location information
    mat neigh = neighbors[0]; //this should always work because there is always a first element of this list
    vec s = neigh.col(0) - 1;
    vec row = conv_to<vec>::from(floor(s / ncols));
    loc_v = row + 1;
    loc_u = s - row*data.n_cols + 1;
  }

  
  bool r_func = true;
  Environment global = Environment::global_env();
  
  std::map<std::string, cdfPtr> cdf_map;
  cdf_map["gaussian_single_param"] = gaussian_single_param_cdf;
  
  // defined samplers in the package
  std::map<std::string, cdfPtr>::iterator it;
  it = cdf_map.find(conditional_cdf);
  if (it != cdf_map.end())
    r_func = false;
  
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
  
  // store in a list
  List sums_nums_loc = List::create(Named("data") = data,
                                    Named("sums") = sums,
                                    Named("nums") = nums);
  if(ncols > 0) {
    sums_nums_loc.push_back(loc_u, "u");
    sums_nums_loc.push_back(loc_v, "v");
  }

  if(r_func) {
    Function cdf = global[conditional_cdf];
    NumericVector resid = cdf(sums_nums_loc, params);
    vec res(resid.begin(), resid.length());
    return(res);
  } else {
    vec res = cdf_map[conditional_cdf](sums_nums_loc, params);
    return(res);
  }

}
