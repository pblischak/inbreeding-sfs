library(tidyverse)

f0.1  <- readr::read_csv("2D_model1_0.1.out")
f0.2 <- readr::read_csv("2D_model1_0.2.out")
f0.3  <- readr::read_csv("2D_model1_0.3.out")
f0.4 <- readr::read_csv("2D_model1_0.4.out")
f0.5  <- readr::read_csv("2D_model1_0.5.out")

df <- rbind(f0.1, f0.2, f0.3, f0.4, f0.5)

df %>% ggplot(aes(x=Nu2_true, y=Nu2_est, group=Nu2_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(. ~ F_true) +
  ggtitle("Pop 2 Size")

ggsave("2D_model1_Nu2.pdf", units = "in", width=12, height=10)

df %>% ggplot(aes(x=T_true, y=T_est, group=T_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(. ~ F_true) +
  ggtitle("Pop 2 Divergence")

ggsave("2D_model1_T.pdf", units = "in", width=12, height=10)

df %>% ggplot(aes(x=M21_true, y=M21_est, group=M21_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(. ~ F_true) +
  ggtitle("Migration (M_21)")

ggsave("2D_model1_M21.pdf", units = "in", width=12, height=10)
