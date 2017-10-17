#include <RcppArmadillo.h>
#include "sampler.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; using namespace arma;

typedef arma::vec (*samplerPtr)(List, List);

//' Run a conclique-based Gibbs sampler to sample spatial data given a lattice and neighborhood structure.
//' 
//' @param conclique_cover A list of class "conclique_cover" encoding the locations in each conclique for 
//'        the conclique cover
//' @param neighbors A list of matrices encoding the neighbors for each point, where the first column of each matrix is the location id of each location in the lattice. This could be the result from get_neighbors().
//'                  Include multiple matrices within the list for more complicated neighborhood structures in your corresponding sampler, for example a (two) dependency parameter for N/S, E/W.
//' @param inits Initial values for the lattice, formatted as a grid.
//' @param conditional_sampler The string name of a function that has two inputs: 
//'        \itemize{
//'          \item{data, and}
//'          \item{params.}
//'        }
//'        There are three built in samplers:
//'         \itemize{
//'           \item{"gaussian_single_param" - a Gaussian sampler with a single dependence parameter,}
//'           \item{"binary_single_param" - a binary sampler with a single dependence parameter, and}
//'           \item{"binary_two_param" - a binary sampler with two dependence parameters.}
//'           \item{"binary_two_param_reg" - a binary sampler with two dependence parameters and regression on the horizontal coordinate.}
//'         }
//'        If the user chooses to write their own sampler in R, they must pass the name of the sampler that is available in the global environment as this parameter.
//'        The input "data" is a list containing at least two elements, 
//'        sums and nums which contain the sum of the data in each neighborhood as well as the number of locations 
//'        in the neighborhood for each point in the conclique. In addition, the data can contain two additional elements, u and v, 
//'        which are vectors that contain the horizontal and vertical location of each point in space.
//'        The input "params" is a list of parameter values. This function returns 
//'        a value sampled from the specified conditional distribution given the data and parameters passed.
//' @param params A list of parameters to be passed to the conditional_sampler function     
//' @param n_iter Number of times to run the Gibbs sampler
//' @export
// [[Rcpp::export]]
arma::mat run_conclique_gibbs(List conclique_cover, List neighbors, arma::mat inits, std::string conditional_sampler, List params, int n_iter = 100) {
  mat result(n_iter, inits.n_elem); 
  mat data = inits;

  int Q = conclique_cover.length();
  int N = neighbors.length();
  List sums(N);
  List nums(N);
  int i, j, n, m;
  
  // location information
  mat neigh = neighbors[0]; //this should always work because there is always a first element of this list
  vec s = neigh.col(0) - 1;
  vec row = conv_to<vec>::from(floor(s / data.n_cols));
  vec loc_u(s.n_elem);
  vec loc_v(s.n_elem);
  loc_v = row + 1;
  loc_u = s - row*data.n_cols + 1;
  
  bool r_func = true;
  Environment global = Environment::global_env();
  Function sampler("identity"); //initialize sampler

  std::map<std::string, samplerPtr> sampler_map;
  sampler_map["gaussian_single_param"] = gaussian_single_param_sampler;
  sampler_map["binary_single_param"] = binary_single_param_sampler;
  sampler_map["binary_two_param"] = binary_two_param_sampler;
  sampler_map["binary_two_param_reg"] = binary_two_param_reg_sampler;
  
  // defined samplers in the package
  std::map<std::string, samplerPtr>::iterator it;
  it = sampler_map.find(conditional_sampler);
  if (it != sampler_map.end())
    r_func = false;

  if(r_func)
    sampler = global[conditional_sampler];

  for(i = 0; i < n_iter; ++i) { // iterations
    for(j = 0; j < Q; ++j) { // concliques
      uvec conc = conclique_cover[j];

      for(n = 0; n < N; ++n) { // neighborhood structures
        mat neigh = neighbors[n];
        int q = neigh.n_cols - 1;
        uvec cols = regspace<uvec>(1,  1,  q);
        vec pre_idx = vectorise(neigh.submat(conc - 1, cols) - 1);
        uvec idx = conv_to<uvec>::from(pre_idx);
        uvec where_nan = find_nonfinite(pre_idx);
        
        mat dat = data.elem(idx);
        dat.elem(where_nan).fill(datum::nan); // keep track of NAs for non regular neighborhoods
        dat.reshape(idx.n_elem/q, q);
        
        vec sums_inner(dat.n_rows);
        vec nums_inner(dat.n_rows);
        for(m = 0; m < dat.n_rows; ++m) { // rowSums
          uvec where = find_finite(dat.row(m));
          rowvec vals = dat.row(m);
          
          sums_inner(m) = sum(vals.elem(where)); // NA handling for non regular neighborhoods
          nums_inner(m) = where.n_elem;
        } // end m
        sums[n] = sums_inner;
        nums[n] = nums_inner;
      } // end n
      
      vec u = loc_u.elem(conc - 1);
      vec v = loc_v.elem(conc - 1);
      
      // store in a list
      List sums_nums_loc = List::create(Named("sums") = sums,
                                        Named("nums") = nums,
                                        Named("u") = u,
                                        Named("v") = v);

      if(r_func) {
        NumericVector new_data = sampler(sums_nums_loc, params);
        vec new_data_vec(new_data.begin(), new_data.length());
        data.elem(conc - 1) = new_data_vec;
      } else {
        vec new_data_vec = sampler_map[conditional_sampler](sums_nums_loc, params);
        data.elem(conc - 1) = new_data_vec;
      }
    } // end j
    data.reshape(1, data.n_elem);
    result.row(i) = data;
  } // end i
  return(result);
  
}

//' Run a sequential Gibbs sampler to sample spatial data given a lattice and neighborhood structure.
//' 
//' @param inits Initial values for the lattice, formatted as a grid.
//' @param neighbors A matrix N*N by (max // neighbors) + 1, where the first column is the location id of each location in the lattice. This could be the result from get_neighbors().
//'                  If NULL, will be calculated within the function.
//' @param conditional_sampler The string name of a function that has two inputs: 
//'        \itemize{
//'          \item{data, and}
//'          \item{params.}
//'        }
//'        There are three built in samplers:
//'         \itemize{
//'           \item{"gaussian_single_param" - a Gaussian sampler with a single dependence parameter,}
//'           \item{"binary_single_param" - a binary sampler with a single dependence parameter, and}
//'           \item{"binary_two_param" - a binary sampler with two dependence parameters.}
//'           \item{"binary_two_param_reg" - a binary sampler with two dependence parameters and regression on the horizontal coordinate.}
//'         }
//'        If the user chooses to write their own sampler in R, they must pass the name of the sampler that is available in the gloabl environment as this parameter.
//'        The input "data" is a list containing two elements, 
//'        sums and nums which contain the sum of the data in each neighborhood as well as the number of locations 
//'        in the neighborhood for each point in the conclique. The input "params" is a list of parameter values. This function returns 
//'        a value sampled from the specified conditional distribution given the data and parameters passed.
//' @param params A list of parameters to be passed to the conditional_sampler function 
//' @param n_iter Number of times to run the Gibbs sampler
//' @export
// [[Rcpp::export]]
arma::mat run_sequential_gibbs(List neighbors, arma::mat inits, std::string conditional_sampler, List params, int n_iter = 100) {
  mat result(n_iter, inits.n_elem);
  mat data = inits;
  
  int N = neighbors.length();
  int M = data.n_elem;
  List sums(N);
  List nums(N);
  uvec rows = uvec(1);
  int i, j, n, m, q;
  
  // location information
  mat neigh = neighbors[0]; //this should always work because there is always a first element of this list
  vec s = neigh.col(0) - 1;
  vec row = conv_to<vec>::from(floor(s / data.n_cols));
  vec loc_u(s.n_elem);
  vec loc_v(s.n_elem);
  loc_v = row + 1;
  loc_u = s - row*data.n_cols + 1;
  
  bool r_func = true;
  Environment global = Environment::global_env();
  Function sampler("identity"); //initialize sampler
  
  std::map<std::string, samplerPtr> sampler_map;
  sampler_map["gaussian_single_param"] = gaussian_single_param_sampler;
  sampler_map["binary_single_param"] = binary_single_param_sampler;
  sampler_map["binary_two_param"] = binary_two_param_sampler;
  sampler_map["binary_two_param_reg"] = binary_two_param_reg_sampler;
  
  // defined samplers in the package
  std::map<std::string, samplerPtr>::iterator it;
  it = sampler_map.find(conditional_sampler);
  if (it != sampler_map.end())
    r_func = false;
  
  if(r_func)
    sampler = global[conditional_sampler];
    
  for(i = 0; i < n_iter; ++i) { // iterations
      for(j = 0; j < M; ++j) { // each data point in sequence
        rows(0) = j;
        
        for(n = 0; n < N; ++n) { // neighborhood structures
          mat neigh = neighbors[n];
          q = neigh.n_cols - 1;
          uvec cols = regspace<uvec>(1,  1,  q);
          
          vec pre_idx = vectorise(neigh.submat(rows, cols) - 1);
          uvec idx = conv_to<uvec>::from(pre_idx);
          uvec where_nan = find_nonfinite(pre_idx);
          
          mat dat = data.elem(idx);
          dat.elem(where_nan).fill(0); // NA for irregular lattices
          
          sums[n] = accu(dat);
          nums[n] = dat.n_elem - where_nan.n_elem;
        } // end n
        
        int u = loc_u(j);
        int v = loc_v(j);
        
        // store in a list
        List sums_nums_loc = List::create(Named("sums") = sums,
                                          Named("nums") = nums,
                                          Named("u") = u,
                                          Named("v") = v);
      
      if(r_func) {
        double new_data = as<double>(sampler(sums_nums_loc, params));
        data(j) = new_data;
      } else {
        double new_data = sampler_map[conditional_sampler](sums_nums_loc, params)(0);
        data(j) = new_data;
      }
      
      
    } // end j
    data.reshape(1, data.n_elem);
    result.row(i) = data;
  } // end i
  
  return(result);
  
}

