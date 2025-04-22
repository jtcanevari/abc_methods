library(dplyr)

sir_stoch <- function(params, init, tstart, tsteps){
  
  beta <- params$beta; gamma <- params$gamma
  
  # Initialize storage for results
  results <- list()
  
  S <- numeric(tsteps)
  I <- numeric(tsteps)
  R <- numeric(tsteps)
  
  S[1] <- init$S
  I[1] <- init$I
  R[1] <- init$R
  
  N <- S[1] + I[1] + R[1]
    
  for (t in 2:tsteps) {
    new_infections <- pmin(rpois(1, lambda = beta * S[t-1] * I[t-1] / N), S[t-1])
    new_recoveries <- pmin(rpois(1, lambda = gamma * I[t-1]), I[t-1])
    
    S[t] <- S[t-1] - new_infections
    I[t] <- I[t-1] + new_infections - new_recoveries
    R[t] <- R[t-1] + new_recoveries
  
  }
    
  results <- data.frame(
    time = tstart:(tstart+tsteps-1),
    S = S,
    I = I,
    R = R
  )

  return(results)
}

# # Parameters
# parameters <- list(beta = 0.3, gamma = 0.1)
# initial_state <- list(S = 99000, I = 1000, R = 0)
# tstart <- 0
# tsteps <- 100
# n_sims <- 100  # Number of simulations
# 
# # Run simulations
# simulations <- lapply(1:n_sims, function(x) {
#   sir_stoch(params = parameters, init = initial_state, tstart = tstart, tsteps = tsteps)
# })
# 
# # Combine results to compute mean and quantiles
# sim_results <- do.call(rbind, lapply(seq_along(simulations), function(i) {
#   sim <- simulations[[i]]
#   sim$sim <- i
#   return(sim)
# }))
# 
# 
# # Compute mean and quantiles for prediction intervals
# 
# # Calculate quantiles
# summary_stats <- sim_results %>%
#   group_by(time) %>%
#   summarise(
#     mean_I = mean(I),
#     median_I = median(I),
#     lower_I = quantile(I, 0.025),
#     upper_I = quantile(I, 0.975)
#   )
# 
# # Plotting
# plot(summary_stats$time, summary_stats$mean_I, type = "l", col = "blue", 
#      ylim = range(c(summary_stats$lower_I, summary_stats$upper_I)), 
#      xlab = "Time", ylab = "Number of Infected Individuals", 
#      main = "Stochastic SIR Model Simulations with 95% Prediction Intervals")
# 
# # Add shaded area for 95% prediction interval
# polygon(c(summary_stats$time, rev(summary_stats$time)), 
#         c(summary_stats$upper_I, rev(summary_stats$lower_I)), 
#         col = rgb(0.2, 0.5, 0.5, 0.5), border = NA)

