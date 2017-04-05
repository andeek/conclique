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
    adj <- adjacent_vertices(lattice, as_ids(V(lattice)))
    max_length <- do.call(max, lapply(adj, FUN = length))
    adj_aug <- lapply(adj, FUN = function(row){ r <- as_ids(row); if(max_length > length(r)) r <- c(r, rep(NA, max_length - length(r))); as.numeric(r) })
    neighs <- do.call(rbind, adj_aug)
    res <- data.frame(vertex = as.numeric(names(adj_aug)))
    res <- cbind(res, neighs)
    names(res)[-1] <- paste0("neighbor_", 1:max_length)
    
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
