library(dplyr)
library(skewt)
library(moments)
#source("R/DefineFmodels.R")

musigma_fs <- function(df, gamma=1) {
  if (is.infinite(df)) {
    M1 <- 2*dnorm(0)
    M2 <- 1
  } else {
    M1 <- dt(0,df=df)*df/(df-1)*2 
    M2 <- df/(df-2)
  }
  dist.mn <- M1*(gamma-1/gamma)
  dist.var <- (M2-M1^2)*((gamma^2)+1/(gamma^2)) + 2*(M1^2) - M2
  list(mu=dist.mn, sigma=sqrt(dist.var))
}


dflist <- c(Inf, 10, 5, 3)
gma <- c(3, 2, 3/2, 21/20, 1)
n0 <- 100000

dat <- expand.grid(zeta=dflist, gamma=gma)

rPIT <-function(n, df, gamma=1) {
  mm <- musigma_fs(df,gamma)
  L <- (rskt(n,df,gamma)-mm$mu)/mm$sigma
  pnorm(L)
}


dat <- dat |> rowwise() |>
  mutate(skw = skewness(rPIT(n0, zeta, gamma)), 
          krt = kurtosis(rPIT(n0, zeta, gamma)))

  