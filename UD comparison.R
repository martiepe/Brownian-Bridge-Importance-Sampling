library(Rhabit)
library(terra)
library(ggplot2)
library(viridis)
library(gridExtra)

source(here::here("functions/utility_functions.R"))
sourceDir("functions")
load_lib(here, dplyr, mvnfast, parallel, terra, ggplot2, viridis, RColorBrewer, tidyr)

# =========================================================
# Helper functions
# =========================================================

raster_to_rhabit <- function(cov) {
  lim <- as.vector(extent(cov))
  r   <- res(cov)
  
  xgrid <- seq(lim[1] + r[1] / 2, lim[2] - r[1] / 2, by = r[1])
  ygrid <- seq(lim[3] + r[2] / 2, lim[4] - r[2] / 2, by = r[2])
  
  z <- t(apply(as.matrix(cov), 2, rev))
  
  list(x = xgrid, y = ygrid, z = z)
}

make_ud_plot <- function(dat, fill_var, fill_limits, fill_name, theme_obj) {
  ggplot(dat, aes(x, y)) +
    geom_raster(aes(fill = .data[[fill_var]])) +
    coord_equal() +
    scale_fill_viridis(name = fill_name, limits = fill_limits) +
    xlab("Easting (km)") +
    ylab("Northing (km)") +
    theme_obj
}

# =========================================================
# Load fitted estimates
# =========================================================

objname <- load(
  "C:/Users/marti/Desktop/Brownian-Bridge-Importance-Sampling/case_study/fitted_estimates/sea_lion_deltamax_studyM=100.Rda"
)

df <- get(objname)

beta <- c(
  median(df$beta1[1:10]),
  median(df$beta2[1:10]),
  median(df$beta3[1:10])
)

# =========================================================
# Prepare tracks
# =========================================================

tracks <- read.csv(here("case_study/SSLpreddat.csv"))

ID <- as.integer(tracks$ID)

time <- as.POSIXct(tracks$time)
time <- as.numeric(time)
time <- (time - min(time)) / 3600

xy <- matrix(c(tracks$x, tracks$y) / 1000, ncol = 2)  # convert to km
xydf <- data.frame(x = xy[, 1], y = xy[, 2])

# =========================================================
# Load and prepare habitat covariates
# =========================================================

hbfull <- brick(here("case_study/aleut_habitat.grd"), values = TRUE)

covlist0 <- list(
  bathy  = hbfull$bathy,
  slope  = hbfull$slope,
  d2site = hbfull$d2site
)

# Convert raster extent/projection from m to km
for (i in seq_along(covlist0)) {
  extent(covlist0[[i]]) <- extent(
    c(
      xmin(covlist0[[i]]), xmax(covlist0[[i]]),
      ymin(covlist0[[i]]), ymax(covlist0[[i]])
    ) / 1000
  )
  
  projection(covlist0[[i]]) <- gsub(
    "units=m", "units=km", projection(covlist0[[i]])
  )
}

# Resample all covariates to the same grid
ncov <- length(covlist0)
for (i in 2:ncov) {
  covlist0[[i]] <- resample(covlist0[[i]], covlist0[[1]])
}

# Crop covariates to area of interest
border <- 30
lim <- c(min(xy[,1])-border,max(xy[,1])+border,min(xy[,2])-border,max(xy[,2])+border)
covlist0 <- lapply(covlist0, crop, y=extent(lim))

covlist <- lapply(covlist0, raster_to_rhabit)

# =========================================================
# UD from supplied beta
# =========================================================

estUDrhabit <- getUD(covariates = covlist, beta = beta)
UD <- rasterToGGplot(estUDrhabit)

# =========================================================
# Model fitting
# =========================================================

gradarray <- bilinearGradArray(locs = xy, cov_list = covlist)

keeptracks <- c(1, 2, 3)
indobs <- which(ID %in% unique(ID)[keeptracks])

keepcovs <- c(1, 2, 3)

fit <- langevinUD(
  locs = xy[indobs, ],
  times = time[indobs],
  ID = ID[indobs],
  grad_array = gradarray[indobs, , keepcovs]
)

estUDrhabit <- getUD(covariates = covlist, beta = fit$betaHat)
UD_michelot <- rasterToGGplot(estUDrhabit)

# =========================================================
# Plot settings
# =========================================================

UD$log_val <- log(UD$val)
UD_michelot$log_val <- log(UD_michelot$val)

ud_michelot_df <- as.data.frame(UD_michelot, xy = TRUE, na.rm = FALSE)
names(ud_michelot_df)[3] <- "val"

ggtheme <- theme(
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 12),
  legend.title = element_text(size = 15),
  legend.text = element_text(size = 12),
  legend.key.height = unit(2, "line"),
  title = element_text(size = 12)
)

lim_pi <- range(c(UD$val, UD_michelot$val), na.rm = TRUE)
lim_logpi <- range(c(UD$log_val, UD_michelot$log_val), na.rm = TRUE)

# =========================================================
# Plots
# =========================================================

p1 <- ggplot() +
  geom_raster(data = ud_df, aes(x = x, y = y, fill = val)) +
  geom_point(data = xydf, aes(x = x, y = y), size = 0.5, color = "grey10") +
  scale_fill_viridis_c(
    trans = "sqrt",
    limits = lim_pi,
    name = expression(pi),
    labels = label_scientific()
  ) +
  coord_equal() +
  labs(x = "Easting (km)", y = "Northing (km)") +
  ggtheme


p3 <- ggplot() +
  geom_raster(data = ud_michelot_df, aes(x = x, y = y, fill = val)) +
  geom_point(data = xydf, aes(x = x, y = y), size = 0.5, color = "grey10") +
  scale_fill_viridis_c(
    trans = "sqrt",
    limits = lim_pi,
    name = expression(pi),
    labels = label_scientific()
  ) +
  coord_equal() +
  labs(x = "Easting (km)", y = "Northing (km)") +
  ggtheme


p2 <- ggplot() +
  geom_raster(data = ud_df, aes(x = x, y = y, fill = log_val)) +
  geom_point(data = xydf, aes(x = x, y = y), size = 0.5, color = "grey10") +
  scale_fill_viridis_c(
    limits = lim_logpi,
    name = expression(log(pi)),
    labels = label_scientific()
  ) +
  coord_equal() +
  labs(x = "Easting (km)", y = "Northing (km)") +
  ggtheme



p4 <- ggplot() +
  geom_raster(data = ud_michelot_df, aes(x = x, y = y, fill = log_val)) +
  geom_point(data = xydf, aes(x = x, y = y), size = 0.5, color = "grey10") +
  scale_fill_viridis_c(
    limits = lim_logpi,
    name = expression(log(pi)),
    labels = label_scientific()
  ) +
  coord_equal() +
  labs(x = "Easting (km)", y = "Northing (km)") +
  ggtheme


grid.arrange(p1, p3, ncol = 2)
grid.arrange(p2, p4, ncol=2)





