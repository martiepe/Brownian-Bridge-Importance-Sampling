# prep workspace ####
source(here::here("functions/utility_functions.R"))
sourceDir("functions")
load_lib(here, dplyr, mvnfast, parallel, terra, ggplot2, viridis, RColorBrewer)

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
dt_max <- 1

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

for (M in Ms){  # for each number of bridges
  for (k in seq_along(deltas)) {  # for each delta_max
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
    setNames(c(paste0("beta", seq_len(npar - 1L)), "sigma", "delta_max", "time")) %>% 
    # save
    save(file =
         sprintf("case_study/fitted_estimates/sea_lion_deltamax_studyM=%s.Rda",
                 M))
}

#### generate summary plots ####
# import estimates
load("case_study/fitted_estimates/sea_lion_deltamax_studyM=25.Rda")
load("case_study/fitted_estimates/sea_lion_deltamax_studyM=50.Rda")
load("case_study/fitted_estimates/sea_lion_deltamax_studyM=100.Rda")

# combine all data
df_all <- bind_rows(mutate(df25, M = "M=25"),
                    mutate(df50, M = "M=50"),
                    mutate(df100, M = "M=100")) %>% 
  pivot_longer(cols = c(beta1, beta2, beta3, # beta4, 
                        gammasq), 
               names_to = "par", values_to = "mu") %>% 
  mutate(par = dplyr::recode_factor(par,
    beta1 = "beta[1]", beta2 = "beta[2]", 
    beta3 = "beta[3]", beta4 = "beta[4]",
    gammasq = "gamma^2"))


# summarise estimates (median, sd, & confidence intervals)
sum_all <- df_all %>%
  dplyr::group_by(par, delta_max, M) %>%
  dplyr::summarise(
    sd = sd(mu),
    mu = median(mu),
    .groups = "drop"
  ) %>% 
  dplyr::mutate(
    lo = mu - z * sd,
    hi = mu + z * sd
  ) 

# define michelot 2019 estimates
michelot_par_est <- data.frame(
  par = c("beta[1]", "beta[2]", "beta[3]", "gamma^2"),
  mu = c(1.34*10^-4,    0.76 *10^-3,   -2.06*10^-5, 12.4),
  lo = c(0.004*10^-4,   -1.74*10^-3,   -3.07*10^-5, 11.9),
  hi = c(2.72*10^-4,    3.25 *10^-3,   -1.05*10^-5, 12.8)
)
# add dumy delta_max for plotting
michelot_par_est <- bind_rows(mutate(michelot_par_est, delta_max = min(sum_all$delta_max)),
          mutate(michelot_par_est, delta_max = max(sum_all$delta_max)))

# generate plot
plot <- ggplot(sum_all, aes(x = delta_max, y = mu, 
                            color = factor(M, levels = c("M=25", "M=50", "M=100")))) +
  # michelot 2019 estimates
  geom_ribbon(data = michelot_par_est, aes(x = delta_max, ymin = lo, ymax = hi),
              fill = "grey20", alpha = 0.1) +
  geom_line(data = michelot_par_est, aes(delta_max, mu), 
            col = "grey20", linetype = "dashed", alpha = 0.5) +
  # BBIS estimates
  geom_point(data = df_all, alpha = 0.15, stroke = NA) +
  geom_line(aes(linetype = factor(M, levels = c("M=25", "M=50", "M=100"))),
            linewidth = 0.7) +
  # design
  facet_wrap(~ par, scales = "free", labeller = label_parsed) + 
  scale_x_log10() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = expression(Delta[max]), y = expression(Estimate),
       color = NULL, linetype = NULL) +
  theme_bw()

# save plot
ggsave(here("case_study", "SSL_par_est_dmax_M.png"), plot)