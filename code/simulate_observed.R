source('code/sir_Gillespie.R')

generate_observed_data <- function(params, init_state, tfinal) {
  SIR_Gillespie(xstart = init_state, params = params, tfinal = tfinal)
}
