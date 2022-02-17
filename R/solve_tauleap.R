#' A routine to extract outputs from a stochastic SEIR model
#'
#' This function implements a stochastic SEIR model that includes three
#' parameters: transmission rate per unit time ($\beta$), average incubation
#' period (1/$\sigma$), and recovery rate ($\gamma$).
#'
#' @param func Function such as a stochastic model implementing tauleap
#' @param params rate of transitioning from exposed to infectious state
#' @param times A sequence of times for which you want to see the output
#' @param y A vector of state variables Parameter values for the stochastic model
#' @return A list of changes in S, E, I, and R
#' @export
solve_tauleap <- function(func, y, times, params, delta) {
  # times indicate the times for which we want to see outputs
  out <- data.frame(matrix(NA, nrow = length(times), ncol = (length(y)+1)))
  out[1, ] <- c(times[1], y)
  row <- 2

  substeps <- round((times[2]-times[1])/delta)
  for (t in 1:(length(times)-1)) {
    for (t2 in 1:substeps) {
      y <- y + unlist(func(y, params, delta))
    }
    out[row, ] <- c(t, y)
    row <- row + 1
  }
  names(out) <- c("time", names(y))
  return (out)
}
