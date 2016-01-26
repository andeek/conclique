#' Run a conclique-based Gibbs sampler to sample spatial data given a lattice and neighborhood structure.
#' 
#' @param data A vector containing data values for each location in the lattice.
#' @param neighbors A matrix N*N by (max # neighbors) + 1, where the first column is the location id of each location in the lattice. This could be the result from get_neighbors().
#' @param conditional_dsn A function that has two inputs: 
#'        \itemize{
#'          \item{data, and}
#'          \item{params.}
#'        }
#'        The data is a list containing data, sums, and nums, which contains the sampled values for each location, the neighboring locations of each location,
#'        and the number of neighbors each location has, respectively. params is a list of parameter values. This function returns the inverse cdf 
#'        at a value between 0 and 1 from the conditional distribution
#' @param params A list of parameters to be passed to the conditional_density function 
#' @export
#' @importFrom igraph V
#' @importFrom igraph as_ids
#' @importFrom igraph adjacent_vertices
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by_
#' @importFrom dplyr do_
#' @importFrom dplyr mutate_
#' @importFrom dplyr arrange_
#' @importFrom tidyr spread_
spatial_residuals <- function(data, neighbors, conditional_dsn, params) {
  cdf_func <- match.fun(conditional_dsn)
  
  data_sums <- rowSums(matrix(data[neighbors[,-1]], ncol = ncol(neighbors) - 1), na.rm = TRUE)
  num_neighbors <- rowSums(!is.na(matrix(data[neighbors[,-1]], ncol = ncol(neighbors) - 1)))
  
  
  cdf_func(data = list(data = data, sums = data_sums, nums = num_neighbors), params = params)
}
