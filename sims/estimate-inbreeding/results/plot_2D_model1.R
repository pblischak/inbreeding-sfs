library(tidyverse)

# F = 0.1
f0.1  <- readr::read_csv("2D_model1_0.1.out")
cat("RMSD for estimated F (F=0.1): ", sqrt(mean((f0.1$F_true - f0.1$F_est)^2)), "\n")

# F = 0.25
f0.25 <- readr::read_csv("2D_model1_0.25.out")
cat("RMSD for estimated F (F=0.25): ", sqrt(mean((f0.25$F_true - f0.25$F_est)^2)), "\n")

# F = 0.5
f0.5  <- readr::read_csv("2D_model1_0.5.out")
cat("RMSD for estimated F (F=0.5): ", sqrt(mean((f0.5$F_true - f0.5$F_est)^2)), "\n")

# F = 0.75
f0.75 <- readr::read_csv("2D_model1_0.75.out")
cat("RMSD for estimated F (F=0.75): ", sqrt(mean((f0.75$F_true - f0.75$F_est)^2)), "\n")

# F = 0.9
f0.9  <- readr::read_csv("2D_model1_0.9.out")
cat("RMSD for estimated F (F=0.9): ", sqrt(mean((f0.9$F_true - f0.9$F_est)^2)), "\n")

df <- rbind(f0.1, f0.25, f0.5, f0.75, f0.9)
df2 <- rbind(f0.1, f0.25, f0.5)
cat("Overall RMSD for F: ", sqrt(mean((df$F_true-df$F_est)^2)), "\n")
cat("Overall RMSD for F (F<=0.5): ", sqrt(mean((df2$F_true-df2$F_est)^2)), "\n")
cat("Overall RMSD for Nu2: ", sqrt(mean((df$Nu2_true-df$Nu2_est)^2)), "\n")
cat("Overall RMSD for Nu2 (F<=0.5): ", sqrt(mean((df2$Nu2_true-df2$Nu2_est)^2)), "\n")
cat("Overall RMSD for T: ", sqrt(mean((df$T_true-df$T_est)^2)), "\n")
cat("Overall RMSD for T (F<=0.5): ", sqrt(mean((df2$T_true-df2$T_est)^2)), "\n")
cat("Overall RMSD for M21: ", sqrt(mean((df$M21_true-df$M21_est)^2)), "\n")
cat("Overall RMSD for M21 (F<=0.5): ", sqrt(mean((df2$M21_true-df2$M21_est)^2)), "\n")
df %>% group_by(F_true,Nu2_true) %>% mutate(Nu2_sq=(Nu2_true-Nu2_est)^2) %>% summarise(rmsd=sqrt(mean(Nu2_sq)))
df %>% group_by(F_true,T_true) %>% mutate(T_sq=(T_true-T_est)^2) %>% summarise(rmsd=sqrt(mean(T_sq)))
df %>% group_by(F_true,M21_true) %>% mutate(M21_sq=(M21_true-M21_est)^2) %>% summarise(rmsd=sqrt(mean(M21_sq)))

df %>% ggplot(aes(x=F_true, y=F_est, group=F_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  ggtitle("Inbreeding Coefficient (Model 3: Divergence + One-Way Migration)")

ggsave("2D_model1_F.pdf", units = "in", width=12, height=10)

df %>% ggplot(aes(x=Nu2_true, y=Nu2_est, group=Nu2_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(. ~ F_true) +
  ggtitle("Pop 2 Size (Model 3: Divergence + One-Way Migration)")

ggsave("2D_model1_Nu2.pdf", units = "in", width=12, height=10)

df %>% ggplot(aes(x=T_true, y=T_est, group=T_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(. ~ F_true) +
  ggtitle("Pop 2 Divergence (Model 3: Divergence + One-Way Migration)")

ggsave("2D_model1_T.pdf", units = "in", width=12, height=10)

df %>% ggplot(aes(x=M21_true, y=M21_est, group=M21_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(. ~ F_true) +
  ggtitle("Migration (Model 3: Divergence + One-Way Migration)")

ggsave("2D_model1_M21.pdf", units = "in", width=12, height=10)
