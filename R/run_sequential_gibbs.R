#' Run a sequential Gibbs sampler to sample spatial data given a lattice and neighborhood structure.
#' 
#' @param lattice The simplified igraph object storing the lattice and dependency structure, ordered by 
#'        location. 
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
#' @importFrom igraph V
#' @importFrom igraph as_ids
#' @importFrom igraph adjacent_vertices
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by_
#' @importFrom dplyr do_
#' @importFrom dplyr mutate_
#' @importFrom dplyr arrange_
#' @importFrom tidyr spread_
run_sequential_gibbs <- function(lattice, inits, conditional_sampler, params, n.iter = 100) {
  
  stopifnot("igraph" %in% class(lattice))
  dimvector <- get.graph.attribute(lattice, "dimvector")
  sampler_func <- match.fun(conditional_sampler)
  data <- array(dim = c(n.iter + 1, prod(dimvector)))
  data[1, ] <- inits
  
  data.frame(vertex = as_ids(V(lattice))) %>%
    group_by_("vertex") %>%
    do_(~data.frame(neighbors = as_ids(adjacent_vertices(lattice, .$vertex)[[1]]))) %>%
    group_by_(~vertex) %>%
    mutate_(key = ~paste0("neighbor_", 1:n())) %>%
    spread_("key", "neighbors") %>%
    arrange_(~vertex) %>%
    data.matrix() -> neighbors
  
  for(i in 1:n.iter) {
    #initialize neighboring data
    data_sums <- rowSums(matrix(data[i, neighbors[,-1]], ncol = ncol(neighbors) - 1), na.rm = TRUE)
    num_neighbors <- rowSums(!is.na(matrix(data[i, neighbors[,-1]], ncol = ncol(neighbors) - 1)))
    
    for(j in 1:ncol(data)) {
      data[i + 1, j] <- sampler_func(data = list(sums = data_sums[j], nums = num_neighbors[j]), params = params)
    }
  }
  
  return(data[-1, ])
}