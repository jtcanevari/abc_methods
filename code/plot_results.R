library(tidyverse)

# Helper to convert particle array into tidy dataframe
get_particle_df <- function(particles, weights) {
  n_steps <- dim(particles)[1]
  n_particles <- dim(particles)[2]
  
  purrr::map_dfr(1:n_steps, function(s) {
    tibble(
      step = s,
      weight = weights[s, ],
      beta = particles[s, , 1],
      gamma = particles[s, , 2]
    )
  })
}

# Posterior density plots
plot_posteriors <- function(result) {
  df_particles <- get_particle_df(result$particles, result$weights)
  
  particle_df_long <- df_particles %>%
    pivot_longer(cols = c(beta, gamma), names_to = "parameter", values_to = "value")
  
  ggplot(particle_df_long, aes(x = value, weight = weight)) +
    geom_density(fill = "steelblue", alpha = 0.5) +
    facet_grid(parameter ~ step, scales = "free_x") +
    labs(title = "Posterior Distributions Across ABC-SMC Steps",
         x = "Parameter value", y = "Density") +
    theme_minimal()
}

# Epidemic fit: observed vs simulated
plot_epidemic_fit <- function(obs_data, result, step = NULL, n_curves = 10) {
  if (is.null(step)) step <- length(result$x_particles)
  
  sims <- result$x_particles[[step]]
  sampled <- sample(sims, n_curves, replace = FALSE)
  
  plot(obs_data$time, obs_data$I, type = 'l', col = 'black', lwd = 2,
       xlab = "Time", ylab = "Infected", main = paste("Observed vs Simulated Epidemics (Step", step, ")"))
  
  for (sim in sampled) {
    lines(sim$time, sim$I, col = rgb(1, 0, 0, 0.3))
  }
  
  legend("topright", legend = c("Observed", "Simulated"),
         col = c("black", rgb(1, 0, 0, 0.3)), lwd = c(2, 1))
}

# Epsilon evolution plot
plot_epsilon <- function(result) {
  plot(result$epsilons, type = 'b', pch = 19,
       xlab = "ABC-SMC Step", ylab = "Tolerance (epsilon)",
       main = "Epsilon Over SMC Steps")
}

# Acceptance rate (if j.it.matrix is available)
plot_acceptance <- function(result) {
  if (!"j.it.s" %in% names(result)) {
    warning("result$j.it.s not found.")
    return()
  }
  
  acceptance <- 1 / rowMeans(result$j.it.s)
  plot(acceptance, type = 'b', pch = 19,
       xlab = "Step", ylab = "Acceptance rate",
       main = "Acceptance Rate Per Step")
}
