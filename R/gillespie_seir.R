#' An Gillespie algorithm (direct method) implementation of a stochastic SEIR model
#'
#' The function \code{gillespie_seir} implements a stochastic SEIR model based on #' a Gillespie algorithm
#' parameters: transmission rate per unit time ($beta$), average incubation
#' period (1/$sigma$), and recovery rate ($gamma$).
#'
#' @param params Parameters \code{$beta$} - transmission rate per unit of time,
#' \code{sigma} rate of transitioning from exposed to infectious state,
#' \code{gamma} rate of transitioning from infectious to recovered state
#' @return A vector of changes in S, E, I, and R
#' @export
gillespie_seir <- function(y, params) {

  S <- y["S"]
  E <- y["E"]
  I <- y["I"]
  R <- y["R"]

  beta <- params["beta"]
  delta <- params["delta"]
  gamma <- params["gamma"]

  N <- S + E + I + R
  dS <- 0
  dE <- 0
  dI <- 0
  dR <- 0
  event_occurred <- FALSE
  if ((I + E) > 0 & S > 0) { ## no need to proceed if no one is infectious or no one
    rate_StoE <- beta * S * I / N
    rate_EtoI <- delta * E
    rate_ItoR <- gamma * I

    rate_all <- c(rate_StoE, rate_EtoI, rate_ItoR) # event rates
    tau <- rexp(1, rate = sum(rate_all)) # time to the next event
    event <- sample(length(rate_all), 1, prob = rate_all) # next event
    if (event == 1) {
      dS <- - 1
      dE <- 1
    }
    else if (event == 2) {
      dE <- - 1
      dI <- 1
    }
    else if (event == 3) {
      dI <- - 1
      dR <- + 1
    }
    event_occurred <- TRUE;
  }
  return (list(y = c(S + dS, E + dE, I + dI, R + dR),
               tau = tau,
               event_occurred = event_occurred))
}
