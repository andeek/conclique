#' Run a conclique-based Gibbs sampler to sample spatial data given a lattice and neighborhood structure.
#' 
#' @param lattice The simplified igraph object storing the lattice and dependency structure, ordered by 
#'        location. 
#' @param conclique_cover A list of class "conclique_cover" encoding the locations in each conclique for 
#'        the conclique cover
#' @param inits Initial values for the lattice, formatted as a grid.
#' @param inv_conditional_dsn A function that has three inputs: 
#'        \itemize{
#'          \item{value,}
#'          \item{data, and}
#'          \item{params.}
#'        }
#'        Value is a vector of values between 0 and 1. The data is a list containing two element, 
#'        sums and nums which contain the sum of the data in each neighborhood as well as the number of locations 
#'        in the neighborhood for each point in the conclique. params is a list of parameter values. This function returns the inverse cdf 
#'        at a value between 0 and 1 from the conditional distribution.
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
run_conclique_gibbs <- function(lattice, conclique_cover, inits, inv_conditional_dsn, params, n.iter = 100) {
  
  stopifnot("igraph" %in% class(lattice) & "conclique_cover" %in% class(conclique_cover))
  dimvector <- get.graph.attribute(lattice, "dimvector")
  inv_func <- match.fun(inv_conditional_dsn)
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
  
  Q <- length(conclique_cover)
  for(i in 1:n.iter) {
    #initialize neighboring data
    data_sums <- rowSums(matrix(data[i, neighbors[,-1]], ncol = ncol(neighbors) - 1), na.rm = TRUE)
    num_neighbors <- rowSums(!is.na(matrix(data[i, neighbors[,-1]], ncol = ncol(neighbors) - 1)))
    
    for(j in 1:Q) {
      data_sums_conc <- data_sums[conclique_cover[[j]]]
      num_neighbors_conc <- num_neighbors[conclique_cover[[j]]]
      U <- runif(length(data_sums_conc))
      data[i + 1, conclique_cover[[j]]] <- inv_func(value = U, data = list(sums = data_sums_conc, nums = num_neighbors_conc), params = params)
    }
  }
  
  return(data[-1, ])
}