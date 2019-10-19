library(tidyverse)

df <- readr::read_csv("2D_model1_0.5.out")
cat("Mean for F: ", mean(df$F_est), "\n")
cat("Overall RMSD for Nu2: ", sqrt(mean((df$Nu2_true-df$Nu2_est)^2)), "\n")
cat("Overall RMSD for T: ", sqrt(mean((df$T_true-df$T_est)^2)), "\n")
cat("Overall RMSD for M21: ", sqrt(mean((df$M21_true-df$M21_est)^2)), "\n")

df %>% ggplot(aes(F_est)) +
  geom_density(fill="gray40", color="white") +
  #geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color="red", size=0.5) +
  theme_bw() +
  facet_grid(Nu2_true ~ T_true) +
  xlim(c(0,1)) +
  ggtitle("Inbreeding Coefficient")

ggsave("2D_model1_F.pdf", units = "in", width=12, height=10)

df %>% ggplot(aes(x=Nu2_true, y=Nu2_est, group=Nu2_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw() +
  facet_grid(T_true ~ M21_true) +
  ggtitle("Pop 2 Size")

ggsave("2D_model1_Nu2.pdf", units = "in", width=12, height=10)

df %>% ggplot(aes(x=T_true, y=T_est, group=T_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw() +
  facet_grid(Nu2_true ~ M21_true) +
  ggtitle("Pop 2 Divergence")

ggsave("2D_model1_T.pdf", units = "in", width=12, height=10)

df %>% ggplot(aes(x=M21_true, y=M21_est, group=M21_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw() +
  facet_grid(Nu2_true ~ T_true) +
  ggtitle("Migration (M_21)")

ggsave("2D_model1_M21.pdf", units = "in", width=12, height=10)
