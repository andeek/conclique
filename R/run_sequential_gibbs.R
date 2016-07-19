#' Run a sequential Gibbs sampler to sample spatial data given a lattice and neighborhood structure.
#' 
#' @param lattice The simplified igraph object storing the lattice and dependency structure, ordered by 
#'        location. 
#' @param inits Initial values for the lattice, formatted as a grid.
#' @param neighbors A matrix N*N by (max # neighbors) + 1, where the first column is the location id of each location in the lattice. This could be the result from get_neighbors().
#'                  If NULL, will be calculated within the function.
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
#' @param directional Indication of if neighborhood needs to be split into North/South, East/West directions. Defaults to FALSE.
#' @param grid A grid storing the locations of each point in the lattice. Only necessary is direction = TRUE. 
#' @param n.iter Number of times to run the Gibbs sampler
#' @export
#' @importFrom igraph get.graph.attribute
#' @importFrom igraph V
#' @importFrom igraph as_ids
#' @importFrom igraph adjacent_vertices
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by_
#' @importFrom dplyr do_
#' @importFrom dplyr mutate_
#' @importFrom dplyr arrange_
#' @importFrom tidyr spread_
run_sequential_gibbs <- function(lattice, neighbors = NULL, inits, conditional_sampler, params, directional = FALSE, grid = NULL, n.iter = 100) {
  
  stopifnot("igraph" %in% class(lattice))
  dimvector <- get.graph.attribute(lattice, "dimvector")
  sampler_func <- match.fun(conditional_sampler)
  
  if(is.null(neighbors)) neighbors <- get_neighbors(lattice, directional, grid)
  result <- array(dim = c(n.iter, prod(dimvector)))
  data <- inits
  
  for(i in 1:n.iter) {
    for(j in 1:length(data)) {
      if(directional) {
        res <- lapply(neighbors, function(x) {
          q <- ncol(x) - 1
          idx <- x[j, -1]
          dat <- data[idx]
          
          list(sum(dat, na.rm = TRUE), sum(!is.na(dat)))
        })
        res <- do.call(rbind, res)
        
        data_sums <- res[, 1]
        num_neighbors <- res[, 2]

      } else {
        idx <- neighbors[j, -1]
        dat <- data[idx]
        
        data_sums <- sum(dat, na.rm = TRUE)
        num_neighbors <- sum(!is.na(dat))
      }
      
      data[j] <- sampler_func(data = list(sums = data_sums, nums = num_neighbors), params = params)
    }
    result[i, ] <- data
  }
  
  return(result)
}