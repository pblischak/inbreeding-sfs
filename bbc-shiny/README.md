# N-Fold Convolution of Beta-Binomials

This folder contains the code for a small R Shiny application to visualize the probability
mass function for the n-fold convolution of beta-binomial distributions used to derived the
expected site frequency spectrum with inbreeding. The easiest way to run the app is to clone
the GitHub repository and to navigate to this folder within RStudio. Then, open up the `app.R`
script and click the **Run App** button in the top right corner of the text editor where the
script is displayed.

```
# 1. Clone the repo if you don't already have it
git clone https://github.com/pblischak/inbreeding-sfs.git

# 2. Navigate to the bbc-shiny folder
cd inbreeding-sfs/bbc-shiny/

# 3. Launch RStudio
rstudio bbc-shiny.Rproj  # if this doesn't work, launch RStudio and navigate to bbc-shiny/
```

Running the Shiny app requires three R packages: `shiny`, `tidyverse`, and `partitions`.
These can be installed from CRAN in the usual way using the `install.packages()` function
in R.

```r
# Install these packages if you don't have them
install.packages("shiny")
install.packages("tidyverse")
install.packages("partitions")
```

![](fig.png)

```r
> sessionInfo()
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] partitions_1.9-19 forcats_0.4.0     stringr_1.4.0     dplyr_0.8.3      
 [5] purrr_0.3.2       readr_1.3.1       tidyr_0.8.3       tibble_2.1.3     
 [9] ggplot2_3.2.1     tidyverse_1.2.1   shiny_1.4.0      

loaded via a namespace (and not attached):
 [1] tidyselect_0.2.5  haven_2.1.1       lattice_0.20-38   colorspace_1.4-1 
 [5] generics_0.0.2    vctrs_0.2.0       htmltools_0.4.0   gmp_0.5-13.5     
 [9] rlang_0.4.0       later_1.0.0       pillar_1.4.2      glue_1.3.1       
[13] withr_2.1.2       modelr_0.1.5      readxl_1.3.1      munsell_0.5.0    
[17] gtable_0.3.0      cellranger_1.1.0  rvest_0.3.4       fastmap_1.0.1    
[21] httpuv_1.5.2      broom_0.5.2       Rcpp_1.0.2        polynom_1.4-0    
[25] xtable_1.8-4      promises_1.1.0    scales_1.0.0      backports_1.1.4  
[29] jsonlite_1.6      mime_0.7          hms_0.5.1         digest_0.6.20    
[33] stringi_1.4.3     grid_3.6.1        cli_1.1.0         tools_3.6.1      
[37] magrittr_1.5      shinythemes_1.1.2 lazyeval_0.2.2    crayon_1.3.4     
[41] pkgconfig_2.0.2   zeallot_0.1.0     xml2_1.2.2        lubridate_1.7.4  
[45] assertthat_0.2.1  httr_1.4.1        rstudioapi_0.10   R6_2.4.0         
[49] nlme_3.1-140      compiler_3.6.1 
```