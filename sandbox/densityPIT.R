library(ggplot2)
library(dplyr)
library(skewt)
source("R/DefineFmodels.R")

df <- Inf
gma <- 3/2

du <- 0.01
u <- seq(0.01,0.99,by=du)
dPIT <- function(u, df, gamma=1) {
  z <- qnorm(u)
  mm <- musigma_fs(df,gamma)
#  mm$sigma*dsst_fs(z*mm$sigma+mm$mu, df, gamma)/dnorm(z)
  dsst_fs(z,df=df,gamma=gamma)/dnorm(z)
}

dPIT2 <- function(u, df, gamma=1) {
  z <- qnorm(u)
  mm <- musigma_fs(df,gamma)
  mm$sigma*dskt(mm$sigma*z+mm$mu,df=df,gamma=gamma)/dnorm(z)
}

# Check dsst by simulation
gma <- 1.5
df <- Inf
Fsim <- rsst_fs(df,gma)(50000) |> ecdf()
mm <- musigma_fs(df, gma)
xx <- seq(-3,4,by=0.01)
fa <- dsst_fs(xx, df, gamma=gma)
z <- (xx-mm$mu)/mm$sigma
fb <- dskt(z, df, gma)/mm$sigma
plot(xx,fa,type='l')
lines(xx,fb,col=2)


# Check PIT density by simulation
Fsim <- rsst_fs(df,gma)(50000) |> pnorm() |> ecdf()
# uu<-seq(0.0001,0.9999,by=0.0001)
# Fa <- cumsum(dPIT(uu,df=df,gamma=gma))*0.0001
# plot(uu,Fsim(uu),type='l')
# lines(uu,Fa,col=2)

Fsim2 <- ((rskt(50000,df,gma)-mm$mu)/mm$sigma) |> pnorm() |> ecdf()
uu<-seq(0.0001,0.9999,by=0.0001)
Fa <- cumsum(dPIT2(uu,df=df,gamma=gma))*0.0001
plot(uu,Fsim2(uu),type='l')
lines(uu,Fa,col=2)
