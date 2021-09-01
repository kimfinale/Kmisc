#' An implementation of a stochastic SEIR model
#'
#' This function implements a stochastic SEIR model that includes three
#' parameters: transmission rate per unit time \code{$\beta$}, average incubation
#' period \code{1/$\sigma$}, and recovery rate \code{$\gamma$}.
#'
#' @param params Parameters \code{$\beta$} - transmission rate per unit of time,
#' \code{$\delta$} rate of transitioning from exposed to infectious state,
#' \code{$\gamma$} rate of transitioning from infectious to recovered state
#' @return A list of changes in S, E, I, and R
#' @export
ode_sepair <- function(t, y, params) {

  S <- y["S"]
  E <- y["E"]
  P <- y["P"]
  A <- y["A"]
  I <- y["I"]
  R <- y["R"]
  CE <- y["CE"] # cumulative infected used to extract over a period (e.g., daily)
  CI <- y["CI"] # cumulative symptom onset to extract number over a period (e.g., daily)


  ## set beta based on the predefined Rt
  ## first set the duration of infectiousness correct
  ## account for the relative infectiousness of P and A state
  durP <- (1 / params[["delta"]] - 1 / params[["epsilon"]])
  durI <- (1 / params[["gamma"]])
  fa <- params[["frac_a"]]
  rhoa <- params[["rho_a"]]
  rhop <- params[["rho_p"]]
  R0_dur <- ((1 - fa) + fa * rhoa) * durI + rhop * durP


  ## now extract the parameters
  epsilon <- params["epsilon"]
  delta <- params["delta"]
  gamma <- params["gamma"]

  rho_p <- params["rho_p"] # relative infectiousness of pre-symptomatic state
  rho_a <- params["rho_a"] # relative infectiousness of asymptomatic state
  frac_a <- params["frac_a"] # fraction of asymptomatic state

  N <- S + E + P + A + I + R
  # beta <- params["beta"]
  beta_t <- get_Rt(t = t) / R0_dur #transmission rate at time t, reflecting already S/N

  dS <- 0
  dE <- 0
  dP <- 0
  dA <- 0
  dI <- 0
  dR <- 0
  dCE <- 0
  dCI <- 0

  rate_from_p <- 1 / (1/delta - 1/epsilon) # rate of transitioning from state P
  rate_StoE <- beta_t * (rho_p * P + rho_a * A + I)
  rate_EtoP <- epsilon * E
  rate_PtoA <- rate_from_p * frac_a * P
  rate_PtoI <- rate_from_p * (1 - frac_a) * P
  rate_AtoR <- gamma * A
  rate_ItoR <- gamma * I

  dS <- - rate_StoE
  dE <-  rate_StoE - rate_EtoP
  dP <-  rate_EtoP - rate_PtoA - rate_PtoI
  dA <-  rate_PtoA - rate_AtoR
  dI <-  rate_PtoI - rate_ItoR
  dR <-  rate_AtoR + rate_ItoR
  dCE <-  rate_StoE
  dCI <-  rate_PtoI

  return (list(c(dS, dE, dP, dA, dI, dR, dCE, dCI)))
}
