#' An implementation of a stochastic SEIR model
#'
#' This function implements a stochastic SEIR model that includes three
#' parameters: transmission rate per unit time \code{$\beta$} and
#' recovery rate \code{$\gamma$}.
#'
#' @param params Parameters \code{$\beta$} - transmission rate per unit of time,
#' \code{$\gamma$} rate of transitioning from infectious to recovered state
#' @return A list of changes in S, E, I, and R
#' @export
ode_sir <- function(t, y, params) {
  S <- y[1]; I <- y[2]; R <- y[3]; CI <- y[4];
  beta <- params["beta"]
  gamma <- params["gamma"]

  N <- S + I + R
  muSI <- beta * I / N
  muIR <- gamma

  dS <- - muSI*S
  dI <- muSI*S - muIR*I
  dR <- muIR*I
  dCI <- muSI*S

  return(list(c(dS, dI, dR, dCI)))
}
