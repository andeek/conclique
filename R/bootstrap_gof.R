#' Get reference distribution for goodness of fit statistics for the spatial residuals of a MRF model.
#' 
#' @param data A vector containing data values for each location in the lattice.
#' @param conclique_cover A list of class "conclique_cover" encoding the locations in each conclique for 
#'        the conclique cover
#' @param neighbors A matrix N*N by (max # neighbors) + 1, where the first column is the location id of each location in the lattice. This could be the result from get_neighbors().
#' @param inits Initial values for the lattice, formatted as a grid.
#' @param conditional_sampler The string name of a function that has two inputs: 
#'        \itemize{
#'          \item{data, and}
#'          \item{params.}
#'        }
#'        There are three built in samplers:
#'         \itemize{
#'           \item{"gaussian_single_param" - a Gaussian sampler with a single dependence parameter,}
#'           \item{"binary_single_param" - a binary sampler with a single dependence parameter, and}
#'           \item{"binary_two_param" - a binary sampler with two dependence parameters.}
#'         }
#'        If the user chooses to write their own sampler in R, they must pass the name of the sampler that is available in the gloabl environment as this parameter.
#'        The input "data" is a list containing two elements, 
#'        sums and nums which contain the sum of the data in each neighborhood as well as the number of locations 
#'        in the neighborhood for each point in the conclique. The input "params" is a list of parameter values. This function returns 
#'        a value sampled from the specified conditional distribution given the data and parameters passed.
#' @param conditional_cdf A function that has two inputs: 
#'        \itemize{
#'          \item{data, and}
#'          \item{params.}
#'        }
#'        There are three built in cdfs:
#'         \itemize{
#'           \item{"gaussian_single_param" - a Gaussian cdf with a single dependence parameter,}
#'           \item{"binary_single_param" - a binary cdf with a single dependence parameter, and}
#'           \item{"binary_two_param" - a binary cdf with two dependence parameters.}
#'         }
#'        If the user chooses to write their own cdf in R, they must pass the name of the cdf function that is available in the global environment as this parameter.
#'        The input "data" is a list containing at least two elements, 
#'        sums and nums which contain the sum of the data in each neighborhood as well as the number of locations 
#'        in the neighborhood for each point in the conclique. In addition, the data can contain two additional elements, u and v, 
#'        which are vectors that contain the horizontal and vertical location of each point in space.
#'        The input "params" is a list of parameter values. This function returns the inverse cdf at a value between 0 and 1 from the conditional distribution
#' @param params A list of fitted parameters to be passed to the conditional_sampler and the conditional_cdf function     
#' @param B The number of bootstrap samples to take, defaults to 1000
#' @param statistic Which goodness of fit statistic to use, Kolmogorov-Smirnov ("ks") or 
#'        Cramer von Mises ("cvm"). Kolmogorov-Smirnov is the default choice. 
#' @param aggregate How to aggregate GOF statistics across concliques, mean ("mean") or 
#'        max ("max"). Mean is the default. Can also be a name of a user defined aggregation function
#' @param quantiles (Optional) A vector of quantiles to return from the reference distribution. If NULL (default), no quantiles are returned
#' @param plot.include TRUE/FALSE parameterizing if a plot of the distribution is returned.
#' 
#' @import ggplot2
#' @importFrom stats quantile
#' @importFrom stats density
#' @export
bootstrap_gof <- function(data, conclique_cover, neighbors, inits, conditional_sampler, conditional_cdf, params, B = 1000, statistic = c("ks", "cvm"), aggregate = c("mean", "max"), quantiles = NULL, plot.include = FALSE) {
  burnin <- 500
  thin <- 5
  iter <- B*thin + burnin
  y_star <- run_conclique_gibbs(conclique_cover, neighbors, inits, conditional_sampler, params, iter)
  y_star <- y_star[(burnin:iter)[burnin:iter %% thin == 1], ]
  resids_star <- apply(y_star, 1, spatial_residuals, neighbors, conditional_cdf, params) %>% t()
  gof_stat_star <- apply(resids_star, 1, gof_statistics, conclique_cover, statistic, aggregate)
  
  resids <- spatial_residuals(data, neighbors, conditional_cdf, params)
  gof_stat <- gof_statistics(resids, conclique_cover, statistic, aggregate)
  
  res <- list()
  
  res$t <- gof_stat
  res$p.value <- (sum(gof_stat_star >= gof_stat) + 1)/(B + 1)
  if(!is.null(quantiles)) {
    res$quantiles <- quantile(gof_stat_star, quantiles)
  }
  
  if(plot.include) {
    res$plot <- ggplot() +
      geom_density(aes(gof_stat_star), fill = "grey20", alpha = .5) +
      geom_vline(aes(xintercept = gof_stat), colour = "red") +
      geom_text(aes(x = gof_stat + diff(range(gof_stat_star))/100, y = max(density(gof_stat_star)$y)*.75, label = paste0("T = ", round(gof_stat, 4), "\np = ", round(res$p.value, 4))), family = "serif", hjust = 0) +
      xlab("GOF Statistic")
  }

  return(res)
}



