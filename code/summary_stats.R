calculate_summary_stats <- function(sim_data) {
  peak_I <- max(sim_data$I)
  time_to_peak <- sim_data$time[which.max(sim_data$I)]
  final_size <- tail(sim_data$R, 1)
  return(c(peak_I, time_to_peak, final_size))
}