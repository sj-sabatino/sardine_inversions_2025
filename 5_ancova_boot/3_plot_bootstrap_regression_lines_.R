library(dplyr)
library(ggplot2)
library(tidyr)

# Load bootstrap ANCOVA results
results <- read.csv("ancova_bootstrap_results.csv")

# Create output directory for plots
output_dir <- "bootstrap_regression_plots"
dir.create(output_dir, showWarnings = FALSE)

# Get list of unique inversions
inversions <- unique(results$inversion)

# Create a set of latitude values (x-axis)
lat_vals <- seq(20, 55, by = 0.5)  # Adjust based on your actual latitude range

for (inv in inversions) {
  # Filter results for current inversion
  inv_results <- results %>% filter(inversion == inv)

  # Split into bootstrap replicates
  model_list <- split(inv_results, inv_results$bootstrap_file)

  # Extract model parameters for each replicate
  model_params <- lapply(model_list, function(df) {
    intercept <- df$estimate[df$term == "(Intercept)"]
    slope_ref <- df$estimate[df$term == "lat"]
    group_shift <- df$estimate[df$term == "groupreference"]
    slope_shift <- df$estimate[df$term == "lat:groupreference"]

    data.frame(
      intercept_ref = intercept,
      slope_ref = slope_ref,
      intercept_alt = intercept + group_shift,
      slope_alt = slope_ref + slope_shift
    )
  })

  params_df <- bind_rows(model_params)

  # Predict values for plotting
  lines_df <- do.call(rbind, lapply(1:nrow(params_df), function(i) {
    data.frame(
      lat = lat_vals,
      freq_ref = params_df$intercept_ref[i] + params_df$slope_ref[i] * lat_vals,
      freq_alt = params_df$intercept_alt[i] + params_df$slope_alt[i] * lat_vals,
      replicate = i
    )
  }))

  # Long format for ggplot
  lines_long <- pivot_longer(
    lines_df,
    cols = c(freq_ref, freq_alt),
    names_to = "group",
    values_to = "freq"
  )
  lines_long$group <- factor(lines_long$group, levels = c("freq_ref", "freq_alt"),
                             labels = c("Reference", "Alternative"))


  # Plot
  p <- ggplot(lines_long, aes(x = lat, y = freq, group = interaction(replicate, group), color = group)) +
    geom_line(alpha = 0.1) +
    stat_summary(aes(group = group), fun = mean, geom = "line", linewidth = 1) +
    scale_color_manual(values = c("black", "blue")) +
    labs(title = paste("Bootstrap Regression Lines -", inv)) +
    ylim(0,1) +
    theme(
      panel.background = element_rect(fill = "grey", color = NA),
      plot.background = element_rect(fill = "grey", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      legend.key = element_rect(fill = "white")
    )

  ggsave(filename = file.path(output_dir, paste0("regression_", inv, ".pdf")),
         plot = p, width = 10, height = 6)

}
