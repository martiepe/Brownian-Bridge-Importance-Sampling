# prep workspace ####
source(here::here("functions/utility_functions.R"))
sourceDir("functions")
load_lib(here, dplyr, mvnfast, parallel, terra, ggplot2, viridis, RColorBrewer,
         ctmm, sf)

# import data ####
data("jaguar")   # load jaguar data from ctmm package
# select 1 track
data <- jaguar[[1]]
tracks <- data.frame(time = data$timestamp,
                     x = data$longitude,
                     y = data$latitude) |> 
  vect(geom = c("x","y"), crs = "EPSG:4326")
