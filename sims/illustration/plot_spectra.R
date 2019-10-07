library(tidyverse)

df <- read_csv("spectra.csv")
df2 <- data.frame(F=as.factor(df$Fis),Frequency=df$X,Expectation=df$Y)

p <- df2 %>% ggplot(aes(Frequency,Expectation,group=F,color=F)) + 
  geom_line(size=0.8,alpha=0.6) +
  geom_point(size=3) +
  theme_bw(base_size = 24) +
  ggtitle("Example Spectra with Inbreeding")

ggsave("inbreeding-spectra.pdf", height=8, width=10)
