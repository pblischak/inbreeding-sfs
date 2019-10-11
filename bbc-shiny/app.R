library(shiny)
library(tidyverse)
library(partitions)

lmultichoose <- function(x){
  # Returns the log of the multinomial
  # coefficient for a vector of integers.
  n <- sum(x)
  denom <- 0
  for(i in x){
    denom <- denom + lfactorial(i)
  }
  return(lfactorial(n) - denom)
}

lbetaBinom <- function(i,p,Fis){
  # Returns the log of the pmf of the beta-binomial
  # distribution.
  if(p == 0){
    p <- 1e-8
  } else if(p == 1){
    p <- 1 - 1e-8
  }
  if(Fis == 0){
    Fis <- 1e-8
  } else if(Fis == 1){
    Fis <- 1 - 1e-8
  }
  alpha <- p * ((1-Fis)/Fis)
  beta  <- (1-p) * ((1-Fis)/Fis)
  res <- lchoose(2,i) + lbeta(i+alpha,2-i+beta) - lbeta(alpha,beta)
  return(res)
}

calc_prt <- function(prt,p,Fis){
  # Calculates the joint probability of sampling the genotypes
  # in a single partition.
  n0 <- sum(prt==0)
  n1 <- sum(prt==1)
  n2 <- sum(prt==2)
  C <- lmultichoose(c(n0,n1,n2))
  res <- C + lbetaBinom(0,p,Fis)*n0 +
    lbetaBinom(1,p,Fis)*n1 +
    lbetaBinom(2,p,Fis)*n2
  return(exp(res))
}

betaBinomConv_pmf <- function(n,p,Fis){
  # This is the main function for the n-fold convolution
  # of beta-binomials. It calculates the full pmf for a sample
  # size of n, allele frequency p, and inbreeding coefficient Fis.
  pmf <- rep(NA,2*n + 1)
  z0 <- rep(0,n)
  z1 <- z0
  z1[1] <- 1
  pmf[1] <- calc_prt(z0,p,Fis)
  pmf[2] <- calc_prt(z1,p,Fis)
  for(d in 2:(2*n)){
    prts <- parts(d)
    if(d <= n){
      good_prts <- prts[,!colSums(prts>2)]
      extra_zeros <- matrix(0, nrow = n-nrow(good_prts),ncol=ncol(good_prts))
      good_prts <- rbind(good_prts,extra_zeros)
    } else {
      good_prts <- prts[1:n, !colSums(prts>2)]
      good_prts <- good_prts[,colSums(good_prts)==d]
    }
    pmf[d+1] <- sum(apply(as.matrix(good_prts),2,calc_prt, p=p, Fis=Fis))
  }
  return(pmf)
}


# Everything after this is just setup for the Shiny app interface.

ui <- fluidPage(
  theme = shinythemes::shinytheme("darkly"),
  titlePanel("N-Fold Convolution of Beta-Binomials"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("p", "Allele frequency", 
                   value = 0.1, min = 0.0, max = 1.0, step=0.01),
      sliderInput("fis", "Inbreeding coefficient", 
                   value = 0.5, min = 0.0, max = 1.0, step=0.01),
      sliderInput("n", "Number of individuals",
                  value = 10, min = 5, max = 20)
    ),
    mainPanel(
      plotOutput("hist")
    )
  )
)

server <- function(input, output, session){
  output$hist <- renderPlot({
    pmf <- betaBinomConv_pmf(input$n,input$p,input$fis)
    plot(0:(2*input$n), pmf, type='h', lwd=5, ylim=c(0,max(pmf)+0.05),
         xlab = "Derived alleles", ylab="Probability",
         main="Beta-Binomial Convolution PMF")
  })
}

shinyApp(ui, server)
