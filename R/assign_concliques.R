#' Assign conclique cover to regular lattice positions
#' 
#' @param lattice The simplified igraph object storing the lattice and dependency structure, ordered by location
#' @param conclique_cover A list of class "conclique_cover" encoding the locations in each conclique for the conclique cover
#' @export
#' @import igraph
#' @examples
#' \dontrun{
#'     lattice <- igraph::make_lattice(c(6,6))
#'     concliques <- min_conclique_cover(lattice)
#'     assign_concliques(lattice, concliques)
#' }


assign_concliques <- function(lattice, conclique_cover) {
  stopifnot("igraph" %in% class(lattice) & "conclique_cover" %in% class(conclique_cover))
  dimvector <- get.graph.attribute(lattice, "dimvector")
  
  grid <- matrix(1:(prod(dimvector)), nrow = dimvector[1])
  for(i in 1:length(conclique_cover)) {
    grid[conclique_cover[[i]]] <- i
  }
  
  class(grid) <- "conclique_grid"
  return(grid)
}