library(tidyverse)

f0.1  <- readr::read_csv("1D_model2_0.1.out")
f0.2 <- readr::read_csv("1D_model2_0.2.out")
f0.3  <- readr::read_csv("1D_model2_0.3.out")
f0.4 <- readr::read_csv("1D_model2_0.4.out")
f0.5  <- readr::read_csv("1D_model2_0.5.out")

df <- rbind(f0.1, f0.2, f0.3, f0.4, f0.5)
cat("Overall RMSD for Nu0: ", sqrt(mean((df$Nu0_true-df$Nu0_est)^2)), "\n")
cat("Overall RMSD for T: ", sqrt(mean((df$T_true-df$T_est)^2)), "\n")
df %>% group_by(F_true,Nu0_true) %>% mutate(Nu0_sq=(Nu0_true-Nu0_est)^2) %>% summarise(rmsd=sqrt(mean(Nu0_sq)))
df %>% group_by(F_true,T_true) %>% mutate(T_sq=(T_true-T_est)^2) %>% summarise(rmsd=sqrt(mean(T_sq)))

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
