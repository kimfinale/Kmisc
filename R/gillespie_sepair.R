#' An Gillespie algorithm (direct method) implementation of a stochastic
#' SEPAIR model
#' The function \code{gillespie_sepair} implements a stochastic SEPAIR model
#' based on #' a Gillespie direct method
#' parameters: transmission rate per unit time ($beta$), mean incubation
#' period (1/$delta$), mean latent period (1/$epsilon$)
#' and recovery rate ($gamma$).
#'
#' @param params Parameters \code{$beta$} - transmission rate per unit of time,
#' \code{delta} rate of transitioning from exposed to onset of symptoms (E->P->I)
#' \code{gamma} rate of transitioning from infectious to recovered state
#' @return A vector of changes in S, E, I, and R
#' @export
gillespie_sepair <- function(y, params) {

  S <- y["S"]
  E <- y["E"]
  P <- y["P"]
  A <- y["A"]
  I <- y["I"]
  R <- y["R"]
  CE <- y["CE"]
  CI <- y["CI"]

  beta <- params["beta"]
  epsilon <- params["epsilon"]
  delta <- params["delta"]
  gamma <- params["gamma"]

  rho_p <- params["rho_p"] # relative infectiousness of pre-symptomatic state
  rho_a <- params["rho_a"] # relative infectiousness of asymptomatic state
  frac_a <- params["frac_a"] # fraction of asymptomatic state
  frac_detect <- params["frac_detect"] # fraction of A and I state that are detected

  N <- S + E + P + A + I + R
  dS <- 0
  dE <- 0
  dP <- 0
  dA <- 0
  dI <- 0
  dR <- 0
  dCE <- 0
  dCI <- 0
  # mean residence time in P state is the mean incubation period (1/delta)
  #  - mean latent period(1/epsilon)
  rate_from_p <- 1 / (1 / delta - 1 / epsilon)
  event_occurred <- FALSE
  tau <- 0
  if ((E + P + A + I) > 0 & S > 0) { ## no need to proceed if no one is infectious
    rate_StoE <- beta * S * (rho_p * P + rho_a * A + I) / N
    rate_EtoP <- epsilon * E
    rate_PtoA <- rate_from_p * frac_a * P
    rate_PtoI <- rate_from_p * (1 - frac_a) * P
    rate_AtoR <- gamma * A
    rate_ItoR <- gamma * I

    rate_all <- c(rate_StoE, rate_EtoP, rate_PtoA, rate_PtoI, rate_AtoR, rate_ItoR) # event rates
    tau <- rexp(1, rate = sum(rate_all)) # time to the next event
    event <- sample(length(rate_all), 1, prob = rate_all) # next event
    if (event == 1) {
      dS <- - 1
      dE <- 1
      dCE <- 1
    }
    else if (event == 2) {
      dE <- - 1
      dP <- 1
    }
    else if (event == 3) {
      dP <- - 1
      dA <- 1
    }
    else if (event == 4) {
      dP <- - 1
      dI <- 1
      dCI <- 1
    }
    else if (event == 5) {
      dA <- - 1
      dR <- 1
    }
    else if (event == 6) {
      dI <- - 1
      dR <- 1
    }
    event_occurred <- TRUE
  }
  return (list(y = c(S + dS, E + dE, P + dP, A + dA, I + dI,
                     R + dR, CE + dCE, CI + dCI),
               tau = tau,
               event_occurred = event_occurred))
}
