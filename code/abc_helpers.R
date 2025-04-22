# Define priors
define_priors <- function() {
  list(
    beta = c("unif", 0.2, 0.8),
    gamma = c("unif", 0.05, 0.2)
  )
}

sample_from_priors <- function(priors) {
  sapply(priors, function(p) runif(1, min = as.numeric(p[2]), max = as.numeric(p[3])))
}

prior_density <- function(theta, priors) {
  prod(mapply(function(val, p) {
    lower <- as.numeric(p[2])
    upper <- as.numeric(p[3])
    if (val >= lower && val <= upper) return(1 / (upper - lower))
    else return(0)
  }, theta, priors))
}

distance_function <- function(sim_stats, obs_stats) {
  sum(abs(sim_stats - obs_stats))  # simple L1 distance
}

# Converts particles and weights into tidy dataframe for plotting
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