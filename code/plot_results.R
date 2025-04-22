library(tidyverse)

#' Plot posterior distributions over ABC-SMC steps
#'
#' Generates density plots for each parameter across ABC-SMC generations.
#' Uses weighted particles from the ABC result object.
#'
#' @param result A list returned by \code{abc_smc()}, containing \code{particles} and \code{weights}.
#'
#' @return A ggplot object showing parameter distributions over time.
#' @export
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

#' Plot observed vs. simulated epidemic curves
#'
#' Overlays epidemic curves from accepted ABC particles on top of the observed data.
#' Useful for visually assessing how well simulated epidemics match the observed trajectory.
#'
#' @param obs_data A dataframe of observed epidemic data with columns \code{time} and \code{I}.
#' @param result A list returned by \code{abc_smc()}, containing the list \code{x_particles}
#'   (simulated epidemic outputs per accepted particle).
#' @param step Integer indicating which ABC-SMC step to plot. Defaults to the last step.
#' @param n_curves Number of simulated epidemic curves to overlay. Default is 10.
#'
#' @return A base R plot (not ggplot) of observed and simulated epidemic curves.
#' @export
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

#' Plot relative summary errors across particles
#'
#' Computes and visualizes the relative error in summary statistics
#' (e.g. peak size, time to peak, final size) for all particles in each step.
#'
#' @param obs_stats Numeric vector of summary statistics from observed data.
#' @param result The result list returned by \code{abc_smc()}.
#' @param summary_names Optional names for the summary statistics.
#'
#' @return A ggplot showing relative errors by summary statistic and ABC step.
#' @export
plot_summary_errors <- function(obs_stats, result, summary_names = NULL) {
  if (is.null(summary_names)) {
    summary_names <- paste0("S", seq_along(obs_stats))
  }
  
  error_df <- purrr::map_dfr(seq_along(result$x_particles), function(s) {
    sims <- result$x_particles[[s]]
    purrr::map_dfr(sims, function(sim) {
      sim_stats <- calculate_summary_stats(sim)
      tibble(step = s, rel_error = (sim_stats - obs_stats) / obs_stats)
    })
  })
  
  error_df <- error_df %>%
    pivot_longer(cols = starts_with("rel_error"), names_to = "summary", values_to = "value") %>%
    mutate(summary = rep(summary_names, each = nrow(error_df) / length(summary_names)))
  
  ggplot(error_df, aes(x = factor(step), y = value, fill = summary)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~summary, scales = "free_y") +
    labs(title = "Relative Error in Summary Statistics Across Steps",
         x = "Step", y = "Relative error") +
    theme_minimal()
}
