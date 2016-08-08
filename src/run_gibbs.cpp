#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; using namespace arma;

//' Run a conclique-based Gibbs sampler with a single dependency parameter to sample spatial data given a lattice and neighborhood structure.
//' 
//' @param lattice The simplified igraph object storing the lattice and dependency structure, ordered by 
//'        location. 
//' @param conclique_cover A list of class "conclique_cover" encoding the locations in each conclique for 
//'        the conclique cover
//' @param neighbors A list of matrices encoding the neighbors for each point, where the first column of each matrix is the location id of each location in the lattice. This could be the result from get_neighbors().
//'                  Include multiple matrices within the list for more complicated neighborhood structures in your corresponding sampler, for example a (two) dependency parameter for N/S, E/W.
//' @param inits Initial values for the lattice, formatted as a grid.
//' @param conditional_sampler A function that has three inputs: 
//'        \itemize{
//'          \item{value,}
//'          \item{data, and}
//'          \item{params.}
//'        }
//'        Value is a vector of values between 0 and 1. The data is a list containing two element, 
//'        sums and nums which contain the sum of the data in each neighborhood as well as the number of locations 
//'        in the neighborhood for each point in the conclique. params is a list of parameter values. This function returns 
//'        a value sampled from the specified conditional distribution given the data and parameters passed.
//' @param params A list of parameters to be passed to the conditional_density function 
//' @param directional Indication of if neighborhood needs to be split into North/South, East/West directions. Defaults to FALSE.
//' @param grid A grid storing the locations of each point in the lattice. Only necessary is direction = TRUE.       
//' @param n.iter Number of times to run the Gibbs sampler
//' 
// [[Rcpp::export]]
arma::mat run_conclique_gibbs(List conclique_cover, List neighbors, arma::mat inits, Function conditional_sampler, List params, int n_iter = 100) {
  mat result(n_iter, inits.n_elem);
  mat data = inits;
  int Q = conclique_cover.length();
  
  for(int i = 0; i < n_iter; ++i) { // iterations
    for(int j = 0; j < Q; ++j) { // concliques
      uvec conc = conclique_cover[j];
      int N = neighbors.length();
      List sums(N);
      List nums(N);
      for(int n = 0; n < N; ++n) { // neighborhood structures
        mat neigh = neighbors[n];
        int q = neigh.n_cols - 1;
        uvec cols = regspace<uvec>(1,  1,  q);
        uvec idx = conv_to<uvec>::from(vectorise(neigh.submat(conc - 1, cols)));
        
        mat dat = data.elem(idx - 1);
        dat.reshape(idx.n_elem/q, q);
        
        vec sums_inner(dat.n_rows);
        vec nums_inner(dat.n_rows);
        for(int m = 0; m < dat.n_rows; ++m) { // rowSums
          sums_inner(m) = sum(dat.row(m)); // need NA handling for non regular neighborhoods
          nums_inner(m) = dat.row(m).n_elem;
        } // end m
        sums[n] = sums_inner;
        nums[n] = nums_inner;
      } // end n
      
      List sums_nums;
      sums_nums["sums"] = sums;
      sums_nums["nums"] = nums;
  
      NumericVector new_data = conditional_sampler(Rcpp::Named("data", sums_nums), Rcpp::Named("params", params));
      vec new_data_vec(new_data.begin(), new_data.length());
      data.elem(conc) = new_data_vec;
      
    } // end j
    result.row(i) = data;
  } // end i
  
  return(result);
  
}
  

