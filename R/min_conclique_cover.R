#' Create a minimal conclique cover for a regular lattice
#' 
#' @param lattice The simplified igraph object storing the lattice and dependency structure, ordered by location
#' @param max_concliques The current list of the largest concliques, defaults to empty initial value
#' @export
#' @import igraph
#' @examples
#' \dontrun{
#'     lattice <- igraph::make_lattice(c(6,6))
#'     concliques <- min_conclique_cover(lattice)
#' }
min_conclique_cover <- function(lattice, max_concliques = list()) {
  stopifnot("igraph" %in% class(lattice))
  
  n <- length(max_concliques)
  included_nodes <- do.call(c, max_concliques)
  missing <- difference(V(lattice), included_nodes)
  
  if(length(missing) != 0) {
    max_concliques[[n + 1]] <- missing[largest_ivs(induced_subgraph(lattice, missing))[[1]]]
    min_conclique_cover(lattice, max_concliques)
  } else {
    class(max_concliques) <- "conclique_cover"
    return(max_concliques)
  }
  
}