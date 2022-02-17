#' An implementation of a stochastic SEIR model
#'
#' The function \code{stoch_seir} implements a stochastic SEIR model that includes three
#' parameters: transmission rate per unit time ($\beta$), average incubation
#' period (1/$\sigma$), and recovery rate ($\gamma$).
#'
#' @param params Parameters \code{$\beta$} - transmission rate per unit of time,
#' \code{$\sigma$} rate of transitioning from exposed to infectious state,
#' \code{$\gamma$} rate of transitioning from infectious to recovered state
#' \code{$\delta$} \code{1/\delta} = average incuation period
#' \code{$\ba$} relative infectiousness of asymptomatic state compared to I
#' \code{$\bp$} relative infectiousness of pre-symptomatic state compared to I
#' \code{$\fa$} fraction of asymptomatic states out of A and I states
#' @return A vector of changes in S, E, I, and R
#' @export
tauleap_sepair <- function (y, params, delta) {
#
#   S <- y["S"]
#   E <- y["E"]
#   P <- y["P"]
#   A <- y["A"]
#   I <- y["I"]
#   R <- y["R"]
#   CE <- y["CE"] # cumulative infected used to extract over a period (e.g., daily)
#   CI <- y["CI"] # cumulative symptom onset to extract number over a period (e.g., daily)
  S <- y[1]; E <- y[2]; P <- y[3]; A <- y[4]; I <- y[5]; R <- y[6];
  CE <- y[7]; CI <- y[8]
  ## now extract the parameters
  epsilon <- params["epsilon"]
  delta <- params["delta"]
  gamma <- params["gamma"]
  ## set beta based on the predefined Rt
  ## first set the duration of infectiousness correct
  ## account for the relative infectiousness of P and A states
  durP <- (1 / params[["delta"]] - 1 / params[["epsilon"]])
  durI <- (1 / params[["gamma"]])
  fa <- params[["fa"]]# fraction of asymptomatic state
  bp <- params["bp"] # relative infectiousness of pre-symptomatic state
  ba <- params["ba"] # relative infectiousness of asymptomatic state
  R0_dur <- ((1 - fa) + fa * ba) * durI + bp * durP

  N <- S + E + P + A + I + R
  beta <- params["beta"]
  if (Rt_varying) {
    beta <- get_Rt(t = t) / R0_dur #transmission rate at time t, reflecting already S/N
  }

  dS <- 0
  dE <- 0
  dP <- 0
  dA <- 0
  dI <- 0
  dR <- 0
  dCE <- 0
  dCI <- 0

  muSE <- beta * I / N
  muEI <- sigma
  muIR <- gamma

  muP <- 1 / (1/delta - 1/epsilon) # rate of transitioning from state P
  muSE <- beta * (bp * P + ba * A + I)
  muEP <- epsilon * E
  muPA <- muP * fa
  muPI <- muP * (1 - fa)
  mAR <- gamma
  muIR <- gamma

  N_SE <- rbinom(1, S, 1 - exp(- muSE * delta))
  N_EP <- rbinom(1, E, 1 - exp(- muEP * delta))
  N_P <- rmultinorm(1, P,
          c(exp(-muP*delta), 1-exp(-muP*delta)*f, 1-exp(-muP*delta)*(1-f)))
  N_PA <- NP[2]
  N_PI <- NP[3]
  N_AR <- rbinom(1, A, 1 - exp(- muIR * delta))
  N_IR <- rbinom(1, I, 1 - exp(- muIR * delta))

  dSdt <- - N_SE
  dEdt <- N_SE - N_EP
  dPdt <- N_EP - N_PA - N_PI
  dAdt <- N_PA - N_AR
  dIdt <- N_PI - N_IR
  dRdt <- N_AR + N_IR
  dCEdt <- N_SE
  dCIdt <- N_PI

 return (list(c(dSdt, dEdt, dPdt, dAdt, dIdt, dRdt, dCEdt, dCIdt)))
}
