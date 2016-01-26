#' Get the neighbors of each location in the lattice.
#' 
#' @param lattice The simplified igraph object storing the lattice and dependency structure, ordered by 
#'        location. 
#' @return A data frame N*N by (max # neighbors) + 1, where the first column is the location id of each location in the lattice.
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
get_neighbors <- function(lattice) {
  
  stopifnot("igraph" %in% class(lattice))
            
  data.frame(vertex = as_ids(V(lattice))) %>%
    group_by_("vertex") %>%
    do_(~data.frame(neighbors = as_ids(adjacent_vertices(lattice, .$vertex)[[1]]))) %>%
    group_by_(~vertex) %>%
    mutate_(key = ~paste0("neighbor_", 1:n())) %>%
    spread_("key", "neighbors") %>%
    arrange_(~vertex) %>%
    data.matrix()
}
