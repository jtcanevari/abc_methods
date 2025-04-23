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
#' Computes and visualises the relative error in summary statistics
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
  
  n_steps <- length(result$x_particles)
  n_summaries <- length(obs_stats)
  
  # Build tidy long-format error dataframe
  error_df <- purrr::map_dfr(1:n_steps, function(s) {
    sims <- result$x_particles[[s]]
    purrr::map_dfr(seq_along(sims), function(j) {
      sim_stats <- calculate_summary_stats(sims[[j]])
      tibble(
        step = s,
        summary = summary_names,
        rel_error = (sim_stats - obs_stats) / obs_stats
      )
    })
  })
  
  ggplot(error_df, aes(x = factor(step), y = rel_error)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~summary, scales = "free_y") +
    labs(title = "Relative Error in Summary Statistics Across Steps",
         x = "Step", y = "Relative error") +
    theme_minimal()
}


#' Plot epidemic fit with posterior predictive envelope
#'
#' Simulates epidemics from particles in the final ABC-SMC step and summarizes
#' infection trajectories using median and quantile envelopes.
#'
#' @param obs_data A dataframe of observed epidemic data with columns \code{time} and \code{I}.
#' @param result The result object returned by \code{abc_smc()}.
#' @param model_func The model simulation function (same as used in ABC).
#' @param initial_state Named numeric vector of initial conditions.
#' @param tfinal Final time for simulation.
#' @param n_sim Number of posterior predictive simulations to run.
#'
#' @return A ggplot of observed vs simulated epidemic with 50% and 90% envelopes.
#' @export
plot_epidemic_envelope <- function(obs_data, result, model_func, initial_state, tfinal, n_sim = 500) {
  # Get particles and weights from final step
  final_step <- length(result$particles[, 1, 1])
  particles <- result$particles[final_step, , ]
  weights <- result$weights[final_step, ]
  
  # Resample particles with replacement
  sampled_indices <- sample(1:nrow(particles), size = n_sim, replace = TRUE, prob = weights)
  sampled_particles <- particles[sampled_indices, ]
  
  time_grid <- 0:tfinal
  
  sim_list <- lapply(1:n_sim, function(i) {
    theta <- sampled_particles[i, ]
    sim <- model_func(initial_state, theta, tfinal)
    I_interp <- approx(x = sim$time, y = sim$I, xout = time_grid, method = "linear", rule = 2)$y
    tibble(time = time_grid, I = I_interp)
  })
  
  # Bind and summarize trajectories
  all_sims <- bind_rows(sim_list, .id = "sim") %>%
    group_by(time) %>%
    summarise(
      median = median(I),
      q05 = quantile(I, 0.05),
      q25 = quantile(I, 0.25),
      q75 = quantile(I, 0.75),
      q95 = quantile(I, 0.95),
      .groups = "drop"
    )
  
  # Prepare observed data with color label
  obs_line <- obs_data %>%
      mutate(Source = "Observed")
  
  # Plot with legend
  ggplot(all_sims, aes(x = time)) +
    geom_ribbon(aes(ymin = q05, ymax = q95, fill = "90% interval"), alpha = 0.3) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = "50% interval"), alpha = 0.4) +
    geom_line(aes(y = median, color = "Posterior median"), linewidth = 1) +
    geom_line(data = obs_line, aes(x = time, y = I, color = "Observed", group = 1),
              linewidth = 1, linetype = "dashed") +
    scale_fill_manual(name = "Posterior interval", values = c("50% interval" = "blue", "90% interval" = "lightblue")) +
    scale_color_manual(name = "Trajectory", values = c("Observed" = "black", "Posterior median" = "blue")) +
    labs(title = "Posterior Predictive Epidemic Fit",
         x = "Time", y = "Infected") +
    theme_minimal()
  
}

