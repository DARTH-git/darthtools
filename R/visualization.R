#---------------------------------------------------------------------------#
####                    R functions for visualization                    ####
#---------------------------------------------------------------------------#

# plot density of total cost
plot_tc <- function(tc) {
  # Histogram showing variability in individual total costs
  plot(density(tc), main = paste("Total cost per person"), xlab = "Cost ($)")
}

# plot density of total QALYs
plot_te <- function(te) {
  # Histogram showing variability in individual total QALYs
  plot(density(te), main = paste("Total QALYs per person"), xlab = "QALYs")
}
