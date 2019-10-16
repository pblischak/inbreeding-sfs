library(tidyverse)

f0.1  <- readr::read_csv("2D_model1_noF_0.1.out")
f0.2 <- readr::read_csv("2D_model1_noF_0.2.out")
f0.3  <- readr::read_csv("2D_model1_noF_0.3.out")
f0.4 <- readr::read_csv("2D_model1_noF_0.4.out")
f0.5  <- readr::read_csv("2D_model1_noF_0.5.out")

df <- rbind(f0.1, f0.2, f0.3, f0.4, f0.5)
cat("Overall RMSD for Nu2: ", sqrt(mean((df$Nu2_true-df$Nu2_est)^2)), "\n")
cat("Overall RMSD for T: ", sqrt(mean((df$T_true-df$T_est)^2)), "\n")
cat("Overall RMSD for M21: ", sqrt(mean((df$M21_true-df$M21_est)^2)), "\n")
df %>% group_by(F_true,Nu2_true) %>% mutate(Nu2_sq=(Nu2_true-Nu2_est)^2) %>% summarise(rmsd=sqrt(mean(Nu2_sq)))
df %>% group_by(F_true,T_true) %>% mutate(T_sq=(T_true-T_est)^2) %>% summarise(rmsd=sqrt(mean(T_sq)))
df %>% group_by(F_true,M21_true) %>% mutate(M21_sq=(M21_true-M21_est)^2) %>% summarise(rmsd=sqrt(mean(M21_sq)))

df %>% ggplot(aes(x=Nu2_true, y=Nu2_est, group=Nu2_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(. ~ F_true) +
  ggtitle("Pop 2 Size")

ggsave("2D_model1_noF_Nu2.pdf", units = "in", width=12, height=10)

df %>% ggplot(aes(x=T_true, y=T_est, group=T_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(. ~ F_true) +
  ggtitle("Pop 2 Divergence")

ggsave("2D_model1_noF_T.pdf", units = "in", width=12, height=10)

df %>% ggplot(aes(x=M21_true, y=M21_est, group=M21_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(. ~ F_true) +
  ggtitle("Migration (M_21)")

ggsave("2D_model1_noF_M21.pdf", units = "in", width=12, height=10)
