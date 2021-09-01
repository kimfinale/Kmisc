#' Simulate the model implemented in the Gillespie's direct method
#'
#' @param fun Model
#' @param tend Simulation time
#' @param nrun The number of simulation runs
#' @param y0 Initial states
#' @param params Parameters for the model
#' @return Simulation results
#' @export
gillespie_run <- function(func, tend, nrun, y0, params, report_dt = 1) {
  sim_res <- list()
  ## set beta based on the predefined Rt
  ## first set the duration of infectiousness correct
  ## account for the relative infectiousness of P and A state
  durP <- (1 / params[["delta"]] - 1 / params[["epsilon"]])
  durI <- (1 / params[["gamma"]])
  fa <- params[["frac_a"]]
  rhoa <- params[["rho_a"]]
  rhop <- params[["rho_p"]]
  R0_dur <- ((1-fa) + fa*rhoa) * durI + rhop * durP

  for (i in seq_len(nrun)) {
    # cat("i =", i, "\n")
    res <- data.frame(time = 0, t(y0))
    t <- 0
    yt <- y0
    dt <- report_dt # output is recorded at t close to dt,
    while (t < tend) {
      params["beta"] <- get_Rt(t = t) / R0_dur
      yt_1 <- c(t, t(yt)) # holds values before report_dt arrives
      sim <- func(y = yt, params = params)
      yt <- sim$y
      t <- t + sim$tau
      if (t > dt){
        res <- rbind(res, yt_1)
        dt <- dt + report_dt # dt
      }
      if (!sim$event_occurred) break
    }
    sim_res[[i]] <- res
  }
  return (sim_res)
}
