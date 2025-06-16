# Load required libraries
library(ggplot2)
library(maps)
library(mapdata)
library(scatterpie)
library("readxl")

# Define geographic bounds (longitude and latitude)
lon1 <- -15.853060
lon2 <- 31.738642
lat1 <- 15.105538
lat2 <- 55.559711

# Define the x and y limits for the map
xlim <- c(lon1, lon2)
ylim <- c(lat1, lat2)

# Set the inversion name
inverson_name <- "INV_5"

# Load geographic data (longitude and latitude) from an Excel file
my_data <- read_excel("lats_longs.xlsx")

# Extract longitude and latitude columns from the loaded data
longitude <- c(my_data$longitude)
latitude <- c(my_data$latitude)

# Define data for pie charts (locations and associated data values)
locations <- data.frame(
  country = c("AEG_A", "AEG_J", "SPA_A", "CAD1_J", "CAD3_J", "CAD3_A", "NPO_A", "CPL_A",
              "SPA_J", "CAN_J", "CAN_A", "NSA_A", "LYG_A", "TAR_A", "TOR_A", "ION_A",
              "ADR_A", "NFR_J", "NFR_A", "CFR_J", "CFR_A", "ENC_A", "BRC_J", "ENC_J",
              "MMO_A", "SMO_J", "NMO_A", "NMO_J", "CMO_A", "CMO_J", "NSS_A", "NSC_A",
              "NSR_A", "NSB_J"),
  longitude = longitude,  # Longitude values loaded from the Excel file
  latitude = latitude,    # Latitude values loaded from the Excel file
  longitude2 = c(25.514, 23.928, -9.03447, -6.792, -7.007, -6.418, -8.759, -9.378,
                 -8.967, -15.349, -16.919, -6.105964, 3.483, 1.386, -0.523, 20.590671,
                 13.016, -2.056, -2.156, -1.31, -1.391, -4.585, -5.138, -3.419,
                 -5.072, -16.452, -6.59, -6.695, -13.568, -9.738, -9.06, -8.45,
                 -7.14, -2.590582),
  latitude2 = c(40.836, 40.618, 37.005169, 36.838, 36.955, 36.398, 40.859, 39.208,
                37.435301, 28.038, 28.521, 43.729712, 42.806, 41.074, 38.05, 38.901515,
                45.201, 46.054, 46.054, 45.466, 44.061, 50.175, 50.725, 50.115,
                35.642, 23.668, 35.676, 34.766, 27.05, 30.115, 42.73, 43.43, 43.75,
                43.561872),
  maj_AF = c(0.689, 0.73, 0.496, 0.486, 0.452, 0.475, 0.583, 0.41, 0.484, 0.863,
             0.917, 0.707, 0.848, 0.857, 0.833, 0.715, 0.855, 0.602, 0.575, 0.575,
             0.6, 0.696, 0.692, 0.686, 0.526, 0.073, 0.539, 0.478, 0.399, 0.526,
             0.482, 0.508, 0.609, 0.704),
  min_AF = c(0.311, 0.270, 0.504, 0.514, 0.548, 0.525, 0.417, 0.590, 0.516, 0.137,
             0.083, 0.293, 0.152, 0.143, 0.167, 0.285, 0.145, 0.398, 0.425, 0.425,
             0.400, 0.304, 0.308, 0.314, 0.474, 0.927, 0.461, 0.522, 0.601, 0.474,
             0.518, 0.492, 0.391, 0.296)
)

# Create two versions of the locations dataset for mapping
# This was done for making points and then pie charts that were subsequently moved manually.
locations2 <- locations  # Clone the dataset
locations <- locations[-c(4, 5, 27), ]  # Remove specific rows from the first dataset
locations2 <- locations2[-c(4, 5, 18, 19, 20, 21, 27), ]  # Remove additional rows from the second dataset

# Create a map of Europe with defined longitude and latitude bounds
europe_map <- map_data('world', xlim = c(lon1, lon2), ylim = c(lat1, lat2), interior = TRUE)

# Define colors for the pie charts
pal_color <- c("yellow", "navy")

# Create the map plot
mapplot1 <- ggplot(europe_map) +
  geom_map(data = europe_map, map = europe_map, aes(map_id = region),
           col = "grey", fill = "lightgrey") +  # Draw the map background
  geom_point(data = locations2[, c(4:5)], aes(x = longitude2, y = latitude2),
             size = 1, shape = 16, fill = "black") +  # Add points for locations
  geom_scatterpie(data = locations, aes(x = longitude, y = latitude, 
                                        group = country, r = 0.75),
                  cols = colnames(locations[, c(6:7)])) +  # Add pie charts
  scale_fill_manual(values = pal_color) +  # Set colors for pie charts
  xlim(-19.853060, 30.738642) +  # Define x-axis limits
  ylim(21.105538, 55.559711) +  # Define y-axis limits
  theme(
    panel.background = element_rect(fill = "white"),  # Set background color
    panel.ontop = FALSE,
    panel.grid = element_blank()  # Remove grid lines
  )

# Save the plot as a PDF file
ggsave(paste(inverson_name, "_piechart.pdf"), plot = mapplot1, width = 10, height = 6, dpi = 400)
