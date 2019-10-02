library(tidyverse)

f0.1  <- readr::read_csv("1D_model2_0.1.out")
f0.25 <- readr::read_csv("1D_model2_0.25.out")
f0.5  <- readr::read_csv("1D_model2_0.5.out")
f0.75 <- readr::read_csv("1D_model2_0.75.out")
f0.9  <- readr::read_csv("1D_model2_0.9.out")

df <- rbind(f0.1, f0.25, f0.5, f0.75, f0.9)

df %>% ggplot(aes(x=F_true, y=F_est, group=F_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  #facet_grid(Nu0_true ~ T_true) +
  ggtitle("Inbreeding Coefficient (Model 2: Bottle-Growth Model)")

ggsave("1D_model2_F.pdf", units = "in", width=12, height=10)


df %>% ggplot(aes(x=Nu0_true, y=Nu0_est, group=Nu0_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(F_true ~ T_true) +
  ggtitle("Bottleneck Size (Model 2: Bottle-Growth Model)")

ggsave("1D_model2_Nu0.pdf", units = "in", width=12, height=10)

df %>% ggplot(aes(x=T_true, y=T_est, group=T_true)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.02), alpha=0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size=0.5) +
  theme_bw(base_size = 20) +
  facet_grid(F_true ~ Nu0_true) +
  ggtitle("Bottleneck Duration (Model 2: Bottle-Growth Model)")

ggsave("1D_model2_T.pdf", units = "in", width=12, height=10)
