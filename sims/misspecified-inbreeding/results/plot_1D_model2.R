library(tidyverse)

df <- readr::read_csv("1D_model2_0.5.out")

df %>% ggplot(aes(F_est)) +
  geom_density(fill="gray40", color="white") +
  #geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(Nu0_true ~ T_true) +
  xlim(c(0,1)) +
  ggtitle("Inbreeding Coefficient (Model 2: Bottle-Growth Model)")

ggsave("1D_model2.pdf", units = "in", width=12, height=10)


df %>% ggplot(aes(x=Nu0_true, y=Nu0_est, group=Nu0_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(. ~ T_true) +
  ggtitle("Bottleneck Size (Model 2: Bottle-Growth Model)")

ggsave("1D_model2_Nu0.pdf", units = "in", width=12, height=10)

df %>% ggplot(aes(x=T_true, y=T_est, group=T_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(. ~ Nu0_true) +
  ggtitle("Bottleneck Duration (Model 2: Bottle-Growth Model)")

ggsave("1D_model2_T.pdf", units = "in", width=12, height=10)
