#' Get the neighbors of each location in the lattice.
#' 
#' @param lattice The simplified igraph object storing the lattice and dependency structure, ordered by 
#'        location. 
#' @param directional Indication of if neighborhood needs to be split into North/South, East/West directions. Defaults to FALSE.
#' @param grid A grid storing the locations of each point in the lattice. Only necessary is direction = TRUE. 
#' @return A list containting data frames N*N by (max # neighbors) + 1, where the first column is the location id of each location in the lattice.
#' @export
#' @importFrom igraph V
#' @importFrom igraph as_ids
#' @importFrom igraph adjacent_vertices
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by_
#' @importFrom dplyr do_
#' @importFrom dplyr mutate_
#' @importFrom dplyr arrange_
#' @importFrom dplyr n
#' @importFrom tidyr spread_
get_neighbors <- function(lattice, directional = FALSE, grid = NULL) {
  
  stopifnot("igraph" %in% class(lattice))
            
  data.frame(vertex = as_ids(V(lattice))) %>%
    group_by_("vertex") -> vertices
  
  if(directional) {
    vertices %>%
      do_(~data.frame(neighbors = adjacent_ew(lattice, .$vertex, grid))) %>%
      group_by_(~vertex) %>%
      mutate_(key = ~paste0("neighbor_", 1:n())) %>%
      spread_("key", "neighbors") %>%
      arrange_(~vertex) %>%
      data.matrix() -> res_ew
    
    vertices %>%
      do_(~data.frame(neighbors = adjacent_ns(lattice, .$vertex, grid))) %>%
      group_by_(~vertex) %>%
      mutate_(key = ~paste0("neighbor_", 1:n())) %>%
      spread_("key", "neighbors") %>%
      arrange_(~vertex) %>%
      data.matrix() -> res_ns
    
    res <- list(ew = res_ew, ns = res_ns)
  } else {
    vertices %>%
      do_(~data.frame(neighbors = as_ids(adjacent_vertices(lattice, .$vertex)[[1]]))) %>%
      group_by_(~vertex) %>%
      mutate_(key = ~paste0("neighbor_", 1:n())) %>%
      spread_("key", "neighbors") %>%
      arrange_(~vertex) %>%
      data.matrix() -> res
    res <- list(res)
  }
  
  return(res)
    
}


#' @importFrom igraph V
#' @importFrom igraph as_ids
#' @importFrom igraph adjacent_vertices
#' @noRd
adjacent_ns <- function(lattice, vertex, grid) {
  
  stopifnot("igraph" %in% class(lattice))
  
  adj_ids <- as_ids(adjacent_vertices(lattice, vertex)[[1]])
  adj_ids[col(grid)[adj_ids] == col(grid)[vertex]]
  
}

#' @importFrom igraph V
#' @importFrom igraph as_ids
#' @importFrom igraph adjacent_vertices
#' @noRd
adjacent_ew <- function(lattice, vertex, grid) {
  
  stopifnot("igraph" %in% class(lattice))
  
  adj_ids <- as_ids(adjacent_vertices(lattice, vertex)[[1]])
  adj_ids[row(grid)[adj_ids] == row(grid)[vertex]]
}
