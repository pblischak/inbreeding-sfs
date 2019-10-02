library(tidyverse)

f0.1  <- readr::read_csv("2D_model1_0.1.out")
f0.25 <- readr::read_csv("2D_model1_0.25.out")
f0.5  <- readr::read_csv("2D_model1_0.5.out")
f0.75 <- readr::read_csv("2D_model1_0.75.out")
f0.9  <- readr::read_csv("2D_model1_0.9.out")

df <- rbind(f0.1, f0.25, f0.5, f0.75, f0.9)

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
