#' Approximate Bayesian Computation via Sequential Monte Carlo (ABC-SMC)
#'
#' Performs ABC-SMC inference for a model with specified priors and summary statistics.
#' Can be used with both stochastic and deterministic models.
#'
#' @param obs_data A dataframe of observed model output (e.g. incidence over time).
#' @param obs_stats A numeric vector of summary statistics from the observed data.
#' @param model_func A function to simulate model output. Should take arguments:
#'   \code{(initial_state, theta, tfinal)}.
#' @param summary_func A function to calculate summary statistics from simulated data.
#' @param priors A named list of prior definitions. Each should be a character vector like
#'   \code{c("unif", lower, upper)}.
#' @param initial_state A named numeric vector of initial conditions for the model.
#' @param tfinal The final simulation time (passed to the model).
#' @param n_particles Number of particles (simulations) per SMC step.
#' @param n_steps Number of SMC steps (generations).
#'
#' @return A list with elements:
#'   \item{particles}{An array of accepted parameter values for each step}
#'   \item{weights}{Normalized weights for each particle}
#'   \item{x_particles}{List of model outputs per accepted particle}
#'   \item{epsilons}{Tolerance thresholds per step}
#' @export
abc_smc <- function(
    obs_data,
    obs_stats,
    model_func,
    summary_func,
    priors,
    initial_state,
    tfinal,
    n_particles = 200,
    n_steps = 5
) {
  # Initialize storage
  particles <- array(NA, dim = c(n_steps, n_particles, length(priors)))
  weights <- matrix(1 / n_particles, n_steps, n_particles)
  distances <- matrix(NA, n_steps, n_particles)
  epsilons <- rep(NA, n_steps)
  
  x_particles <- vector("list", n_steps)
  for (s in 1:n_steps) {
    x_particles[[s]] <- vector("list", n_particles)
  }
  
  # Tolerance for step 1
  epsilons[1] <- 1
  # epsilons[1] <- 1e5  # arbitrarily large to guarantee acceptance in step 1
  
  for (step in 1:n_steps) {
    cat("Running step", step, "\n")
    
    if (step > 1) {
      epsilons[step] <- quantile(distances[step - 1, ], probs = 0.75)
    }
    
    for (j in 1:n_particles) {
      repeat {
        # Propose parameters
        if (step == 1) {
          theta <- sample_from_priors(priors)
        } else {
          idx <- sample(1:n_particles, 1, prob = weights[step - 1, ])
          theta_star <- particles[step - 1, idx, ]
          
          theta <- numeric(length(priors))
          for (k in seq_along(priors)) {
            lower <- as.numeric(priors[[k]][2])
            upper <- as.numeric(priors[[k]][3])
            sd_k <- sd(particles[step - 1, , k])
            repeat {
              theta[k] <- rnorm(1, mean = theta_star[k], sd = sd_k)
              if (theta[k] >= lower && theta[k] <= upper) break
            }
          }
        }
        
        # Simulate model
        sim <- tryCatch({
          model_func(initial_state, theta, tfinal)
        }, error = function(e) return(NULL))
        
        if (is.null(sim)) next
        
        sim_stats <- summary_func(sim)
        dist <- distance_function(sim_stats, obs_stats)
        
        if (dist <= epsilons[step]) {
          particles[step, j, ] <- theta
          distances[step, j] <- dist
          x_particles[[step]][[j]] <- sim
          
          if (step == 1) {
            weights[step, j] <- 1 / n_particles
          } else {
            prior_prob <- prior_density(theta, priors)
            denom <- 0
            for (i in 1:n_particles) {
              prod_kernel <- prod(mapply(
                function(tk, tk_prev, sd_k) {
                  dnorm(tk, mean = tk_prev, sd = sd_k)
                },
                theta,
                particles[step - 1, i, ],
                apply(particles[step - 1, , ], 2, sd)
              ))
              denom <- denom + weights[step - 1, i] * prod_kernel
            }
            weights[step, j] <- prior_prob / denom
          }
          break
        }
      }
    }
    
    # Normalize weights
    weights[step, ] <- weights[step, ] / sum(weights[step, ])
  }
  
  return(list(
    particles = particles,
    weights = weights,
    distances = distances,
    epsilons = epsilons,
    x_particles = x_particles
  ))
}
