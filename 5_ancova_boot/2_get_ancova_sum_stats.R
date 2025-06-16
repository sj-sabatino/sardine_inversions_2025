library(dplyr)
library(ggplot2)

# Read in combined results from all bootstraps
results <- read.csv("/Users/sabo/Desktop/ancova_bootstrap_results.csv")

# Focus on the interaction term
interaction_term <- "lat:groupreference"

# Get list of unique inversions
inversions <- unique(results$inversion)

# Create a directory for plots
output_dir <- "interaction_effect_plots"
dir.create(output_dir, showWarnings = FALSE)

# Loop through each inversion
for (inv in inversions) {
  inv_results <- results %>%
    filter(term == interaction_term, inversion == inv)

  # Summary statistics
  n_total <- nrow(inv_results)
  n_significant <- sum(inv_results$p.value < 0.05)
  prop_significant <- n_significant / n_total

  mean_effect <- mean(inv_results$estimate)
  sd_effect <- sd(inv_results$estimate)
  ci_95 <- quantile(inv_results$estimate, probs = c(0.025, 0.975))

  # Console output
  cat("\nSummary for Inversion:", inv, "\n")
  cat("--------------------------------------\n")
  cat("Total bootstraps: ", n_total, "\n")
  cat("Significant p-values (p < 0.05): ", n_significant, "\n")
  cat("Proportion significant: ", round(prop_significant, 3), "\n\n")

  cat("Mean effect size: ", round(mean_effect, 3), "\n")
  cat("Standard deviation: ", round(sd_effect, 3), "\n")
  cat("95% CI of effect: [", round(ci_95[1], 3), ", ", round(ci_95[2], 3), "]\n")

  # Plot
  p <- ggplot(inv_results, aes(x = estimate)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = ci_95, linetype = "dotted", color = "darkred") +
    labs(title = paste("Interaction Effect Sizes -", inv),
         subtitle = paste("Term:", interaction_term),
         x = "Effect Size", y = "Count") +
    theme_minimal()

  ggsave(filename = file.path(output_dir, paste0("interaction_", inv, ".pdf")),
         plot = p, width = 8, height = 5)
}
