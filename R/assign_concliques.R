#' Assign conclique cover to regular lattice positions
#' 
#' @param lattice The simplified igraph object storing the lattice and dependency structure, ordered by location
#' @param grid A grid storing the locations of each point in the lattice. Either lattice or grid must be specified.
#' @param conclique_cover A list of class "conclique_cover" encoding the locations in each conclique for the conclique cover
#' @export
#' @importFrom igraph get.graph.attribute
#' @examples
#' \dontrun{
#'     lattice <- igraph::make_lattice(c(6,6))
#'     concliques <- min_conclique_cover(lattice)
#'     assign_concliques(lattice, concliques)
#' }
assign_concliques <- function(lattice = NULL, grid = NULL, conclique_cover) {
  #check for either lattice or grid
  stopifnot("conclique_cover" %in% class(conclique_cover) & (!is.null(lattice) | !is.null(grid)))
  if(is.null(grid)) {
    stopifnot("igraph" %in% class(lattice))
    dimvector <- get.graph.attribute(lattice, "dimvector")
    grid <- matrix(1:(prod(dimvector)), nrow = dimvector[1])
  } 
  
  for(i in 1:length(conclique_cover)) {
    grid[conclique_cover[[i]]] <- i
  }
  
  class(grid) <- "conclique_grid"
  return(grid)
}