#' Functions for using conclique based GOF test for Gaussian MRF with a single dependence parameter.
#' 
#' @param data A list containing two elements, 
#'        sums and nums which contain the sum of the data in each neighborhood as well as the number of locations 
#'        in the neighborhood for each point in the conclique. 
#' @param params A list of parameter values, rho, kappa, and eta, that parameterize the Gaussian MRF.
#' @name gaussian_single_param
NULL


#' @rdname gaussian_single_param
#' @export
gaussian_single_param_sampler <- function(data, params) {
  rho <- params$rho
  kappa <- params$kappa
  eta <- params$eta
  
  mean_structure <- kappa + eta*(data$sums - data$nums*kappa)
  rnorm(length(mean_structure))*rho + mean_structure
}


#' @rdname gaussian_single_param
#' @export
gaussian_single_param_cdf <- function(data, params) {
  rho <- params$rho
  kappa <- params$kappa
  eta <- params$eta
  
  mean_structure <- kappa + eta*(data$sums - data$nums*kappa)
  pnorm(data$data, mean = mean_structure, sd = rho)
}