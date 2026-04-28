# prep workspace ####
source(here::here("functions/utility_functions.R"))
sourceDir("functions")
load_lib(here, dplyr, mvnfast, parallel, terra, ggplot2, viridis, RColorBrewer, tidyr)

# import data ####
hbfull <- rast(here("case_study/aleut_habitat.grd"))
tracks <- read.csv(here("case_study/SSLpreddat.csv")) %>% 
  mutate(time = as.POSIXct(time)) %>% 
  vect(geom = c("x","y"), crs = crs(hbfull))

# Change the resolution and extent from m to km
crs <- gsub("units=m","units=km",crs(hbfull, T))
r <- rast(nrows = nrow(hbfull), ncols = ncol(hbfull),
          ext = as.vector(ext(hbfull))/1000,
          crs = crs) # define template raster

# standardise projections
hbfull <- project(hbfull, r)     # transform raster
tracks <- project(tracks, r) %>% # transform tracks
  as.data.frame(geom = "XY")

# drop distance to shelf due to colliniarity with distance to SSL sites
hbfull <- hbfull[[c("bathy", "slope", "d2site")]]

#### fit Langevin BBIS - single fit####
ncores <- 10
M <- 25
dt_max <- 0.01

# quick single fit
out <- fit_langevin_bbis(tracks, hbfull, 
                         M = M,
                         dt_max = dt_max, 
                         dt_units = "hours",
                         ncores = ncores, 
                         fixed_sampling = FALSE) 


#### fit Langevin BBIS - dt_max & M refits ####
# define fitting criteria
ncores <- 10
Ms <- c(25, 50, 100)
deltas <- exp(seq(log(0.01), log(25), length.out = 30))
nrefits <- 10

# number of pars 
npar <- nlyr(hbfull) + 1
# add 1 column to track dt_max
params <- matrix(NA, ncol = npar + 2, nrow = nrefits * length(deltas))


for (M in Ms){  
  for (k in 1:length(deltas)) {  
    for (i in 1:nrefits) {  
      delta_max <- deltas[k]
      
      print(sprintf("Fitting M = %s, delta_max = %.4f, refit = %s",M, delta_max, i))
      # fit
      out <- fit_langevin_bbis(tracks, hbfull, 
                               M = M,
                               dt_max = delta_max, 
                               dt_units = "hours",
                               ncores = ncores, 
                               fixed_sampling = FALSE) 
      # store outputs (par + delta_max)
      params[(k - 1L) * nrefits + i, ] <- c(out$par, delta_max, as.numeric(out$time, units = "secs"))


    }
  }
  # convert to data.frame 
  as.data.frame(params) %>% 
    # update names 
    setNames(c(paste0("beta", seq_len(npar - 1L)), "gammasq", "delta_max", "time")) %>% 
    # save
    save(file =
         sprintf("case_study/fitted_estimates/sea_lion_deltamax_studyM=%s.Rda",
                 M))
}


#### generate summary plots ####
# import estimates

objs <- load("case_study/fitted_estimates/sea_lion_deltamax_studyM=25.Rda")
df25 <- get(objs)

objs <- load("case_study/fitted_estimates/sea_lion_deltamax_studyM=50.Rda")
df50 <- get(objs)

objs <- load("case_study/fitted_estimates/sea_lion_deltamax_studyM=100.Rda")
df100 <- get(objs)



library(dplyr)
library(tidyr)
library(ggplot2)

# combine all data
df_all <- bind_rows(
  mutate(df25, M = "M=25"),
  mutate(df50, M = "M=50"),
  mutate(df100, M = "M=100")
) %>% 
  pivot_longer(
    cols = c(beta1, beta2, beta3, gammasq),
    names_to = "par",
    values_to = "mu"
  ) %>% 
  mutate(
    par = dplyr::recode_factor(
      par,
      beta1   = '"Bathymetry"',
      beta2   = '"Slope"',
      beta3   = '"Distance to SSL sites"',
      gammasq = 'gamma^2'
    )
  )

z <- qnorm(0.95)

# summarise estimates (median, sd, & confidence intervals)
sum_all <- df_all %>%
  group_by(par, delta_max, M) %>%
  summarise(
    sd = sd(mu),
    mu = median(mu),
    .groups = "drop"
  ) %>% 
  mutate(
    lo = mu - z * sd,
    hi = mu + z * sd
  )

# Michelot 2019 estimates
michelot_par_est <- data.frame(
  par = c('"Bathymetry"', '"Slope"', '"Distance to SSL sites"', "gamma^2"),
  mu = c(1.34e-4, 0.76e-3, -2.06e-5, 12.4),
  lo = c(0.004e-4, -1.74e-3, -3.07e-5, 11.9),
  hi = c(2.72e-4, 3.25e-3, -1.05e-5, 12.8)
)


michelot_par_est <- bind_rows(
  mutate(michelot_par_est, delta_max = min(sum_all$delta_max)),
  mutate(michelot_par_est, delta_max = max(sum_all$delta_max))
)

# levels for legend
M_levels <- c("M=25", "M=50", "M=100")

# plot
plot <- ggplot(sum_all, aes(x = delta_max, y = mu)) +
  # Michelot 2019 estimates
  geom_ribbon(
    data = michelot_par_est,
    aes(x = delta_max, ymin = lo, ymax = hi),
    inherit.aes = FALSE,
    fill = "grey20",
    alpha = 0.1,
    show.legend = FALSE
  ) +
  geom_line(
    data = michelot_par_est,
    aes(x = delta_max, y = mu),
    inherit.aes = FALSE,
    colour = "grey20",
    linetype = "dashed",
    alpha = 0.5,
    show.legend = FALSE
  ) +
  
  # BBIS estimates
  geom_point(
    data = df_all,
    aes(color = factor(M, levels = M_levels)),
    alpha = 0.15,
    stroke = NA,
    show.legend = FALSE
  ) +
  geom_line(
    aes(
      color = factor(M, levels = M_levels),
      linetype = factor(M, levels = M_levels)
    ),
    linewidth = 0.7
  ) +
  
  facet_wrap(
    ~ par,
    scales = "free",
    labeller = label_parsed
  ) +
  
  scale_x_log10() +
  
  # blue -> red
  scale_color_manual(
    values = c(
      "M=25" = "#2C7BB6",
      "M=50" = "#7B3294",
      "M=100" = "#D7191C"
    ),
    breaks = M_levels,
    name = NULL
  ) +
  scale_linetype_manual(
    values = c(
      "M=25" = "33",
      "M=50" = "55",
      "M=100" = "solid"
    ),
    breaks = M_levels,
    name = NULL
  ) +
  
  labs(
    x = expression(Delta[max]),
    y = expression(Estimate)
  ) +
  
  theme_bw()

plot
