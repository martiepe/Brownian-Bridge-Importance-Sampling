library(ggplot2)
library(gridExtra)

#varying dt, fixed n_obs esimates plots
load("simulation_study/varying_thin_estimates.Rda")


p1 <- ggplot(data = df, aes(x = dt, y = beta1)) +
  geom_boxplot() +
  geom_hline(yintercept  = 4, color = "red", linetype = 2) +
  labs(title = expression(beta[1]), y = expression(beta[1]), x = expression(Delta[t])) +
  theme_bw()

p2 <- ggplot(data = df, aes(x = dt, y = beta2)) +
  geom_boxplot() +
  geom_hline(yintercept  = 2, color = "red", linetype = 2) +
  labs(title = expression(beta[2]), y = expression(beta[2]), x = expression(Delta[t])) +
  theme_bw()

p3 <- ggplot(data = df, aes(x = dt, y = beta3)) +
  geom_boxplot() +
  geom_hline(yintercept  = -0.1, color = "red", linetype = 2) +
  labs(title = expression(beta[3]), y = expression(beta[3]), x = expression(Delta[t])) +
  theme_bw()

p4 <- ggplot(data = df, aes(x = dt, y = gammasq)) +
  geom_boxplot() +
  geom_hline(yintercept  = 5, color = "red", linetype = 2) +
  labs(title = expression(gamma^2), y = expression(gamma^2), x = expression(Delta[t])) +
  theme_bw()

grid.arrange(p1,p2,p3,p4)






#varying thin, fixed Tmax estimates plots

load("simulation_study/varying_thin_estimates_fixed_Tmax.Rda")

p1 <- ggplot(data = df, aes(x = dt, y = beta1)) +
  geom_boxplot() +
  geom_hline(yintercept  = 4, color = "red", linetype = 2) +
  labs(title = expression(beta[1]), y = expression(beta[1]), x = expression(Delta[t])) +
  theme_bw()

p2 <- ggplot(data = df, aes(x = dt, y = beta2)) +
  geom_boxplot() +
  geom_hline(yintercept  = 2, color = "red", linetype = 2) +
  labs(title = expression(beta[2]), y = expression(beta[2]), x = expression(Delta[t])) +
  theme_bw()

p3 <- ggplot(data = df, aes(x = dt, y = beta3)) +
  geom_boxplot() +
  geom_hline(yintercept  = -0.1, color = "red", linetype = 2) +
  labs(title = expression(beta[3]), y = expression(beta[3]), x = expression(Delta[t])) +
  theme_bw()

p4 <- ggplot(data = df, aes(x = dt, y = gammasq)) +
  geom_boxplot() +
  geom_hline(yintercept  = 5, color = "red", linetype = 2) +
  labs(title = expression(gamma^2), y = expression(gamma^2), x = expression(Delta[t])) +
  theme_bw()

grid.arrange(p1,p2,p3,p4)


#varying M estimates plots

load("simulation_study/varying_M_estimates.Rda")



p1 <- ggplot(data = df, aes(x = M, y = beta1)) +
  geom_boxplot() +
  geom_hline(yintercept  = 4, color = "red", linetype = 2) +
  labs(title = expression(beta[1]), y = expression(beta[1])) +
  theme_bw()

p2 <- ggplot(data = df, aes(x = M, y = beta2)) +
  geom_boxplot() +
  geom_hline(yintercept  = 2, color = "red", linetype = 2) +
  labs(title = expression(beta[2]), y = expression(beta[2])) +
  theme_bw()

p3 <- ggplot(data = df, aes(x = M, y = beta3)) +
  geom_boxplot() +
  geom_hline(yintercept  = -0.1, color = "red", linetype = 2) +
  labs(title = expression(beta[3]), y = expression(beta[3])) +
  theme_bw()

p4 <- ggplot(data = df, aes(x = M, y = gammasq)) +
  geom_boxplot() +
  geom_hline(yintercept  = 5, color = "red", linetype = 2) +
  labs(title = expression(gamma^2), y = expression(gamma^2)) +
  theme_bw()

grid.arrange(p1,p2,p3,p4)







#varying N estimates plots
load("simulation_study/varying_N_estimates.Rda")



p1 <- ggplot(data = df, aes(x = N, y = beta1)) +
  geom_boxplot() +
  geom_hline(yintercept  = 4, color = "red", linetype = 2) +
  labs(title = expression(beta[1]), y = expression(beta[1])) +
  theme_bw()

p2 <- ggplot(data = df, aes(x = N, y = beta2)) +
  geom_boxplot() +
  geom_hline(yintercept  = 2, color = "red", linetype = 2) +
  labs(title = expression(beta[2]), y = expression(beta[2])) +
  theme_bw()

p3 <- ggplot(data = df, aes(x = N, y = beta3)) +
  geom_boxplot() +
  geom_hline(yintercept  = -0.1, color = "red", linetype = 2) +
  labs(title = expression(beta[3]), y = expression(beta[3])) +
  theme_bw()

p4 <- ggplot(data = df, aes(x = N, y = gammasq)) +
  geom_boxplot() +
  geom_hline(yintercept  = 5, color = "red", linetype = 2) +
  labs(title = expression(gamma^2), y = expression(gamma^2)) +
  theme_bw()


grid.arrange(p1,p2,p3,p4)
