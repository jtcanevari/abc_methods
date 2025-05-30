# Load dependencies
source("code/sir_gillespie.R")
source("code/abc_helpers.R")
source('code/abc_smc.R')
source("code/plot_results.R")

set.seed(123)

#-------------------------------
# Step 0: Simulate observed data
true_params <- c(beta = 0.5, gamma = 0.1)
initial <- c(time = 0, S = 999, I = 1, R = 0)
tfinal <- 100

repeat {
  obs_data <- generate_observed_data(true_params, initial, tfinal)
  if (max(obs_data$I) > 10) break  # or some other reasonable threshold
}
obs_stats <- calculate_summary_stats(obs_data)

#-------------------------------
# Step 1: run abc smc routine
system.time(
  result <- abc_smc(
    obs_data = obs_data,
    obs_stats = obs_stats,
    model_func = SIR_Gillespie,
    summary_func = calculate_summary_stats,
    priors = define_priors(),
    initial_state = initial,
    tfinal = 100,
    n_particles = 200,
    n_steps = 7
  )
)

#-------------------------------
# Step 3: plot
plot_posteriors(result)
plot_epidemic_fit(obs_data, result)
plot_epsilon(result)
plot_summary_errors(obs_stats, result)
plot_epidemic_envelope(
  n_sim = 50,
  obs_data = obs_data,
  result = result,
  model_func = SIR_Gillespie,  
  initial_state = initial,
  tfinal = 100
)

median(result$particles[5,,1])
median(result$particles[5,,2])
