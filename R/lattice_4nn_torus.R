#neigbor structure functions -----------------------
rowdiff <- function(y, z, mat) abs(row(mat)[y] - row(mat)[z])
coldiff <- function(y, z, mat) abs(col(mat)[y] - col(mat)[z])

neighbors_ew <- function(y, z, mat){ 
  if(col(mat)[y] == 1) {
    #left boundary
    (rowdiff(y, z, mat) == 0 & col(mat)[z] == ncol(mat)) | (rowdiff(y, z, mat) == 0 & coldiff(y, z, mat) == 1)
  } else if(col(mat)[y] == ncol(mat)) {
    #right boundary
    (rowdiff(y, z, mat) == 0 & col(mat)[z] == 1) | (rowdiff(y, z, mat) == 0 & coldiff(y, z, mat) == 1)
  } else {
    #interior
    rowdiff(y, z, mat) == 0 & coldiff(y, z, mat) == 1 
  }
  
}
neighbors_ns <- function(y, z, mat){ 
  if(row(mat)[y] == 1) {
    #top boundary
    (coldiff(y, z, mat) == 0 & row(mat)[z] == nrow(mat)) | (rowdiff(y, z, mat) == 1 & coldiff(y, z, mat) == 0 )
  } else if(row(mat)[y] == nrow(mat)) {
    #bottom boundary
    (coldiff(y, z, mat) == 0 & row(mat)[z] == 1) | (rowdiff(y, z, mat) == 1 & coldiff(y, z, mat) == 0 )
  } else {
    #interior
    rowdiff(y, z, mat) == 1 & coldiff(y, z, mat) == 0 
  }  
}
is_neighbor_4nn_torus <- function(y, z, mat) {
  neighbors_ns(y, z, mat) | neighbors_ew(y, z, mat)
}


#' Create a regular lattice of specified size with four-nearest neighborhood dependence structure on a torus
#' 
#' @param dimvec A 2-dimensional vector specifying the size of the lattice
#' @export
#' @importFrom dplyr rowwise
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate_
#' @importFrom dplyr filter_
#' @importFrom igraph V
#' @importFrom igraph as_ids
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph simplify
#' @importFrom igraph permute.vertices
#' @importFrom igraph set.graph.attribute
#' @importFrom igraph remove.vertex.attribute
#' @importFrom lazyeval interp
#' @examples
#' \dontrun{
#'     lattice <- lattice_4nn_torus(c(6,6))
#' }
lattice_4nn_torus <- function(dimvec) {
  stopifnot(length(dimvec) == 2)
  
  grid <- matrix(1:prod(dimvec), nrow = dimvec[1])
  
  expand.grid(index1 = 1:(prod(dimvec)), index2 = 1:(prod(dimvec))) %>%
    rowwise() %>%
    mutate_(edge = interp(~is_neighbor_4nn_torus(index1, index2, grid), .values = list(index1 = as.name("index1"),
                                                                                                 index2 = as.name("index2"),
                                                                                                 grid = as.name("grid")))) %>%
    filter_("edge") -> lattice_as_graph
  
  lattice <- simplify(graph_from_data_frame(lattice_as_graph, directed = FALSE, vertices = NULL))
  lattice <- permute.vertices(lattice, as.numeric(as_ids(V(lattice)))) #reorder for printing
  lattice <- set.graph.attribute(lattice, "dimvector", dimvec) #important step for assigning the concliques to a grid
  lattice <- remove.vertex.attribute(lattice, "name")
  
  return(lattice)
}



