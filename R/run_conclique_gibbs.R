#' Run a conclique-based Gibbs sampler to sample spatial data given a lattice and neighborhood structure.
#' 
#' @param lattice The simplified igraph object storing the lattice and dependency structure, ordered by 
#'        location
#' @param conclique_cover A list of class "conclique_cover" encoding the locations in each conclique for 
#'        the conclique cover
#' @param inits Initial values for the lattice, formatted as a grid.
#' @param inv_conditional_dsn A function that has three inputs: 
#'        \itemize{
#'          \item{value,}
#'          \item{data, and}
#'          \item{params.}
#'        }
#'        Value is a vector of values between 0 and 1. The data is a data frame containing two columns, 
#'        vertex and sum_neighbor, which encode the values for the neighboring locations of each location 
#'        in the conclique and params is a list of parameter values. This function returns the inverse cdf 
#'        at a value between 0 and 1 from the conditional distribution
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
#' @importFrom dplyr filter_
run_conclique_gibbs <- function(lattice, conclique_cover, inits, inv_conditional_dsn, params, n.iter = 100) {
  stopifnot("igraph" %in% class(lattice) & "conclique_cover" %in% class(conclique_cover))
  dimvector <- get.graph.attribute(lattice, "dimvector")
  inv_func <- match.fun(inv_conditional_dsn)
  data <- array(dim = c(n.iter + 1, dimvector))
  data[1, , ] <- inits

  Q <- length(conclique_cover)
  for(i in 1:n.iter) {
    data.frame(vertex = as_ids(V(lattice))) %>%
      group_by_("vertex") %>%
      do_(~data.frame(sum_neighbor = sum(data[i, , ][as_ids(adjacent_vertices(lattice, .$vertex)[[1]])]),
                      num_neighbor = length(as_ids(adjacent_vertices(lattice, .$vertex)[[1]])))) -> neighbor_data
    
    for(j in 1:Q) {
      neighbor_data %>% filter_(.dots = paste0("vertex %in% c(", paste(conclique_cover[[j]], collapse = ", "), ")")) -> neighbor_data_conc
      U <- runif(length(conclique_cover[[j]]))
      data[i + 1, , ][conclique_cover[[j]]] <- inv_func(value = U, data = neighbor_data_conc, params = params)
    }
  }
  
  return(data[-1, , ])
}