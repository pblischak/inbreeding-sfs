# Check if tidyverse is installed
if("tidyverse" %in% rownames(installed.packages()) == FALSE){
  print("Package viridis is not installed. Installing now...")
  install.packages("viridis")
}

# Check if viridis is installed
if("viridis" %in% rownames(installed.packages()) == FALSE){
  print("Package viridis is not installed. Installing now...")
  install.packages("viridis")
}

# Load libraries
library(tidyverse)
library(viridis)

# Read in spectra, reorder, and convert Fis to factor
df <- read_csv("spectra.csv")
df <- df[order(df$Fis),]
df2 <- data.frame(F=as.factor(df$Fis),Frequency=df$X,Expectation=df$Y)

# Make plot and save
p <- df2 %>% ggplot(aes(Frequency,Expectation,group=F,color=F)) + 
  geom_line(size=0.8,alpha=0.6) +
  geom_point(size=5) +
  theme_bw(base_size = 24) +
  scale_color_viridis(end=0.95,discrete=TRUE) +
  ggtitle("Example Spectra with Inbreeding")

ggsave("inbreeding-spectra.pdf", height=8, width=12)
