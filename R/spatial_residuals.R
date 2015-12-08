#' Run a conclique-based Gibbs sampler to sample spatial data given a lattice and neighborhood structure.
#' 
#' @param data_grid The grid of data values for each location in the lattice
#' @param lattice The simplified igraph object storing the lattice and dependency structure, ordered by 
#'        location
#' @param conditional_dsn A function that has two inputs: 
#'        \itemize{
#'          \item{data, and}
#'          \item{params.}
#'        }
#'        The data is a data frame containing four columns, vertex, data_value, sum_neighbor, and num_neighbor which encode the values for each location and the neighboring locations of each location 
#'        and the number of neighbors each location has. params is a list of parameter values. This function returns the inverse cdf 
#'        at a value between 0 and 1 from the conditional distribution
#' @param params A list of parameters to be passed to the conditional_density function 
#' @export
#' @importFrom igraph get.graph.attribute
#' @importFrom igraph V
#' @importFrom igraph as_ids
#' @importFrom igraph adjacent_vertices
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by_
#' @importFrom dplyr do_
spatial_residuals <- function(data_grid, lattice, conditional_dsn, params) {
  stopifnot("igraph" %in% class(lattice))
  dimvector <- get.graph.attribute(lattice, "dimvector")
  cdf_func <- match.fun(conditional_dsn)
  vertices <- as_ids(V(lattice))
  
  data.frame(vertex = vertices) %>%
    group_by_("vertex") %>%
    do_(~data.frame(data_value = data_grid[.$vertex],
                    sum_neighbor = sum(data_grid[as_ids(adjacent_vertices(lattice, .$vertex)[[1]])]),
                    num_neighbor = length(as_ids(adjacent_vertices(lattice, .$vertex)[[1]])))) -> neighbor_data
  
  cdf_func(data = neighbor_data, params = params)
}
