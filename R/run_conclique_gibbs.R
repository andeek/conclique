#' Run a conclique-based Gibbs sampler to sample spatial data given a lattice and neighborhood structure.
#' 
#' @param lattice The simplified igraph object storing the lattice and dependency structure, ordered by 
#'        location. 
#' @param conclique_cover A list of class "conclique_cover" encoding the locations in each conclique for 
#'        the conclique cover
#' @param neighbors A matrix N*N by (max # neighbors) + 1, where the first column is the location id of each location in the lattice. This could be the result from get_neighbors().
#'                  If NULL, will be calculated within the function.
#' @param inits Initial values for the lattice, formatted as a grid.
#' @param conditional_sampler A function that has three inputs: 
#'        \itemize{
#'          \item{value,}
#'          \item{data, and}
#'          \item{params.}
#'        }
#'        Value is a vector of values between 0 and 1. The data is a list containing two element, 
#'        sums and nums which contain the sum of the data in each neighborhood as well as the number of locations 
#'        in the neighborhood for each point in the conclique. params is a list of parameter values. This function returns 
#'        a value sampled from the specified conditional distribution given the data and parameters passed.
#' @param params A list of parameters to be passed to the conditional_density function 
#' @param n.iter Number of times to run the Gibbs sampler
#' @export
#' @importFrom igraph get.graph.attribute
run_conclique_gibbs <- function(lattice, conclique_cover, neighbors = NULL, inits, conditional_sampler, params, n.iter = 100) {
  
  stopifnot("igraph" %in% class(lattice) & "conclique_cover" %in% class(conclique_cover))
  dimvector <- get.graph.attribute(lattice, "dimvector")
  sampler_func <- match.fun(conditional_sampler)
  data <- array(dim = c(n.iter + 1, prod(dimvector)))
  data[1, ] <- inits
  
  if(is.null(neighbors)) neighbors <- get_neighbors(lattice)
  
  Q <- length(conclique_cover)
  for(i in 1:n.iter) {
    #initialize neighboring data
    data_sums <- rowSums(matrix(data[i, neighbors[,-1]], ncol = ncol(neighbors) - 1), na.rm = TRUE)
    num_neighbors <- rowSums(!is.na(matrix(data[i, neighbors[,-1]], ncol = ncol(neighbors) - 1)))
    
    for(j in 1:Q) {
      data_sums_conc <- data_sums[conclique_cover[[j]]]
      num_neighbors_conc <- num_neighbors[conclique_cover[[j]]]
      data[i + 1, conclique_cover[[j]]] <- sampler_func(data = list(sums = data_sums_conc, nums = num_neighbors_conc), params = params)
    }
  }
  
  return(data[-1, ])
}