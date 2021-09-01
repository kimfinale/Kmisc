#' An implementation of a stochastic SEIR model
#'
#' This function implements a stochastic SEIR model that includes three
#' parameters: transmission rate per unit time \code{$\beta$}, average incubation
#' period \code{1/$\sigma$}, and recovery rate \code{$\gamma$}.
#'
#' @param params Parameters \code{$\beta$} - transmission rate per unit of time,
#' \code{$\sigma$} rate of transitioning from exposed to infectious state,
#' \code{$\gamma$} rate of transitioning from infectious to recovered state
#' @return A list of changes in S, E, I, and R
#' @export
ode_seir <- function(t, y, params) {
  S <- y[1]; E <- y[2]; I <- y[3]; R <- y[4]; CE <- y[5];  CI <- y[6]  ## now extract the parameters
  beta <- params["beta"]
  sigma <- params["sigma"]
  gamma <- params["gamma"]

  N <- S + E + I + R
  muSE <- beta * I / N
  muEI <- sigma
  muIR <- gamma

  dS <- - muSE*S
  dE <-  muSE*S - muEI*E
  dI <-  muEI*E - muIR*I
  dR <-  muIR*I
  dCE <-  muSE*S
  dCI <-  muEI*E

  return(list(c(dS, dE, dI, dR, dCE, dCI)))
}
