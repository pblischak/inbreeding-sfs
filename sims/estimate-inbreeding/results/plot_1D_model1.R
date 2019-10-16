library(tidyverse)

# F = 0.1
f0.1  <- readr::read_csv("1D_model1_0.1.out")
cat("RMSD for estimated F (F=0.1): ", sqrt(mean((f0.1$F_true - f0.1$F_est)^2)), "\n")

# F = 0.25
f0.25 <- readr::read_csv("1D_model1_0.25.out")
cat("RMSD for estimated F (F=0.25): ", sqrt(mean((f0.25$F_true - f0.25$F_est)^2)), "\n")

# F = 0.5
f0.5  <- readr::read_csv("1D_model1_0.5.out")
cat("RMSD for estimated F (F=0.5): ", sqrt(mean((f0.5$F_true - f0.5$F_est)^2)), "\n")

# F = 0.75
f0.75 <- readr::read_csv("1D_model1_0.75.out")
cat("RMSD for estimated F (F=0.75): ", sqrt(mean((f0.75$F_true - f0.75$F_est)^2)), "\n")

# F = 0.9
f0.9  <- readr::read_csv("1D_model1_0.9.out")
cat("RMSD for estimated F (F=0.9): ", sqrt(mean((f0.9$F_true - f0.9$F_est)^2)), "\n")

df <- rbind(f0.1, f0.25, f0.5, f0.75, f0.9)
cat("Overall RMSD for F: ", sqrt(mean((df$F_true-df$F_est)^2)), "\n")

df %>% ggplot(aes(x=F_true, y=F_est, group=F_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  ggtitle("Inbreeding Coefficient (Model 1: Standard Neutral Model)")

ggsave("1D_model1_F.pdf", units = "in", width=12, height=10)
