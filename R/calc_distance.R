#' A routine to extract outputs from a stochastic SEIR model
#'
#' This function implements a stochastic SEIR model that includes three
#' parameters: transmission rate per unit time ($\beta$), average incubation
#' period (1/$\sigma$), and recovery rate ($\gamma$).
#'
#' @param data Function such as a stochstic SEIR model
#' @param sim Parameter values for the stochastic model
#' @return Sum of squred distance between the model and the data
#' @export
calc_distance <- function(data, sim) {
  distance <- sqrt(sum((data - sim)^2))
  return (distance)
}
