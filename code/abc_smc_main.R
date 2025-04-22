# Load dependencies
source("code/sir_gillespie.R")
source("code/simulate_observed.R")
source("code/summary_stats.R")
source("code/abc_helpers.R")

set.seed(123)

# Step 0: Simulate observed data
true_params <- c(beta = 0.5, gamma = 0.1)
initial <- c(time = 0, S = 999, I = 1, R = 0)
tfinal <- 100

repeat {
  obs_data <- generate_observed_data(true_params, initial, tfinal)
  if (max(obs_data$I) > 10) break  # or some other reasonable threshold
}
obs_stats <- calculate_summary_stats(obs_data)

# Step 1: run abc smc routine
source('code/abc_smc.R')

result <- abc_smc(
  obs_data = obs_data,
  obs_stats = obs_stats,
  model_func = SIR_Gillespie,
  summary_func = calculate_summary_stats,
  priors = define_priors(),
  initial_state = initial,
  tfinal = 100,
  n_particles = 200,
  n_steps = 5
)

# Step 3: plot
source("code/plot_results.R")

plot_posteriors(result)
plot_epidemic_fit(obs_data, result)
plot_epsilon(result)

#--------------------

