# Plot cdf of PIT under FS-Normal misspecification
library(ggplot2)
library(patchwork)
library(dplyr)
library(skewt)
source("R/DefineFmodels.R")

df <- Inf
gma <- 3/2
legfont <- 11

du <- 0.005
u <- seq(0.005,0.995,by=du)
# dPIT <- function(u, df, gamma=1) {
#   z <- qnorm(u)
#   dsst_fs(z,df=df,gamma=gamma)/dnorm(z)
# }
dPIT <- function(u, df, gamma=1) {
  z <- qnorm(u)
  mm <- musigma_fs(df,gamma)
  mm$sigma*dskt(mm$sigma*z+mm$mu,df=df,gamma=gamma)/dnorm(z)
}

pPIT <- function(u, df, gamma=1) {
  z <- qnorm(u)
  mm <- musigma_fs(df,gamma)
  pskt(mm$sigma*z+mm$mu,df=df,gamma=gamma)
}

rPIT <-function(n, df, gamma=1) {
  mm <- musigma_fs(df,gamma)
  L <- (rskt(n,df,gamma)-mm$mu)/mm$sigma
  pnorm(L)
}

# To verify by simulation
#pSim <- rPIT(20000,5,gma) |> ecdf()
#plot(u, pPIT(u,5,gma), type='l')
#lines(u,pSim(u),col=2)

plotpPIT <- function(degfr) {
  dat <- bind_rows(
    tibble(u=u,cdf=pPIT(u,df=degfr,gamma=gma), gamma='3/2'),
    tibble(u=u,cdf=pPIT(u,df=degfr,gamma=1/gma), gamma='2/3')
  )
  if (is.infinite(degfr)) {
    dflab <- "df==infinity"
  } else {
    dflab <- paste0("df==",as.character(degfr))
  }
  p1 <- ggplot(dat, aes(x=u,y=cdf,color=gamma, group=gamma)) + 
    geom_line() + xlim(c(0,1)) + ylim(c(0,1)) +
    geom_abline(intercept=0, slope=1, linetype=3) + 
    theme_bw() +
    labs(x='PIT', y='cdf', color=expression(gamma))  + 
    theme(legend.position = "bottom", legend.direction="horizontal",
          legend.title=element_text(size=legfont), 
          legend.text=element_text(size=legfont)) +
    annotate('text', x=0.82,y=0.1, label=dflab, parse=TRUE)
             #parse=is.infinite(degfr))
  p1
}

pN <- plotpPIT(Inf)
p5 <- plotpPIT(5)

# ggsave mishandles the gamma in legend
grDevices::cairo_pdf("output/PITcdf.pdf", width=7, height=3.5)
pN+p5
dev.off()
