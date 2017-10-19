#' Get goodness of fit statistics for the spatial residuals of a MRF model.
#' 
#' @param residuals A vector containing residuals for each location in the lattice, 
#'        output from spatial_residuals.
#' @param conclique_cover A list of class "conclique_cover" encoding the locations in each conclique for 
#'        the conclique cover
#' @param statistic Which goodness of fit statistic to use, Kolmogorov-Smirnov ("ks") or 
#'        Cramer von Mises ("cvm"). Kolmogorov-Smirnov is the default choice. 
#' @param aggregate How to aggregate GOF statistics across concliques, mean ("mean") or 
#'        max ("max"). Mean is the default. Can also be a name of a user defined aggregation function.
#' @export
gof_statistics <- function(residuals, conclique_cover, statistic = c("ks", "cvm"), aggregate = c("mean", "max")) {
  stopifnot("conclique_cover" %in% class(conclique_cover))
  stat <- statistic[1]
  aggregate <- aggregate[1]
  resids <- lapply(conclique_cover, FUN = function(conc) { sort(residuals[conc]) })
  N <- length(residuals)
  
  if(stat == "ks") {
    s <- lapply(resids, FUN = function(conc) { 
      n <- length(conc)
      tmp <- seq_len(n)/n
      sqrt(N)*max(max(abs(conc - tmp)), max(abs(conc - tmp + 1/n))) 
      })
  } else if (stat == "cvm") {
    s <- lapply(resids, FUN = function(conc) {
      n <- length(conc)
      conc <- c(0, conc, 1)
      res <- array(NA, n + 1)
      for(i in 2:(n + 2)) {
        res[i-1] <- ((i - 2)/n)^2*(conc[i] - conc[i - 1]) - (i - 2)/n*(conc[i]^2 - conc[i - 1]^2) + 1/3
      }
      sqrt(N)*sum(res)
    })
  } else {
    #attempt to evaluate user-input statistic
  }
  do.call(aggregate, s)
}