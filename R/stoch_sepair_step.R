#' An implementation of a stochastic SEIR model
#'
#' The function \code{stoch_seir} implements a stochastic SEIR model that includes three
#' parameters: transmission rate per unit time ($\beta$), average incubation
#' period (1/$\sigma$), and recovery rate ($\gamma$).
#'
#' @param params Parameters \code{$\beta$} - transmission rate per unit of time,
#' \code{$\sigma$} rate of transitioning from exposed to infectious state,
#' \code{$\gamma$} rate of transitioning from infectious to recovered state
#' @return A vector of changes in S, E, I, and R
#' @export
stoch_seir <- function (y, params, delta) {
  beta <- params["beta"]
  sigma <- params["sigma"]
  gamma <- params["gamma"]

  S <- y[1]; E <- y[2]; I <- y[3]; R <- y[4]; CE <- y[5]; CI <- y[6]

  N <- S + E + I + R
  muSE <- beta * I / N
  muEI <- sigma
  muIR <- gamma

  N_SE <- rbinom(1, S, 1 - exp(- muSE * delta))
  N_EI <- rbinom(1, E, 1 - exp(- muEI * delta))
  N_IR <- rbinom(1, I, 1 - exp(- muIR * delta))

  dSdt <- - N_SE
  dEdt <- N_SE - N_EI
  dIdt <- N_EI - N_IR
  dRdt <- N_IR
  dCEdt <- N_SE
  dCIdt <- N_EI

 return (list(c(dSdt, dEdt, dIdt, dRdt, dCEdt, dCIdt)))
}
