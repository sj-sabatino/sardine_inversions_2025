library(ggplot2)
library(car)
library(rstatix)
library(broom)

# List of bootstrap files
#bootstrap_files <- list.files("/Users/sabo/Desktop/bootstrap_ancova/bootstrap_files/boots1", pattern = "*.csv", full.names = TRUE)
bootstrap_files <- list.files("/Users/sabo/Desktop/bootstrap_ancova/boot_test/", pattern = "*.csv", full.names = TRUE)


results_list <- list()

inversions <-list(
"freq_inv1",
"freq_inv2",
"freq_inv3",
"freq_inv4",
"freq_inv5",
"freq_inv6",
"freq_inv7",
"freq_inv8",
"freq_inv9",
"freq_inv10"
)

for (inversion in inversions){
  for (f in bootstrap_files) {
    data <- read.csv(f)
    data_atl <- data[data$ptype == "ATL", ]

    ref <- data.frame(lat = data_atl$latitude, freq = data_atl$freq_ref, group = "reference")
    alt <- data.frame(lat = data_atl$latitude, freq = data_atl[[inversion]], group = "alternative")

    combined <- rbind(ref, alt)
    model <- lm(freq ~ lat * group, data = combined)
    result <- tidy(model)
    result$bootstrap_file <- basename(f)
    result$inversion <- inversion

    results_list[[paste0(basename(f), "_", inversion)]] <- result
  }
}


# Combine all results into a single data frame
all_results <- do.call(rbind, results_list)

# Save final results
write.csv(all_results, "ancova_bootstrap_results.csv", row.names = FALSE)
