# Plot cdf of PIT and transformed PIT under FS-t misspecification
library(ggplot2)
library(patchwork)
library(dplyr)
library(skewt)
source("R/DefineFmodels.R")

df <- 5
gma <- 6/5
legfont <- 11

du <- 0.0005
u <- seq(0,0.9995,by=du)

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

pPITdelta <- function(u, df, gamma=1, dlta=1/2) {
  lu <- dlta*(1-u)
  pPIT(lu+u,df,gamma)-pPIT(lu,df,gamma)
}

rPIT <-function(n, df, gamma=1) {
  mm <- musigma_fs(df,gamma)
  L <- (rskt(n,df,gamma)-mm$mu)/mm$sigma
  pnorm(L)
}

rPITdelta <-function(n, df, gamma=1, delta=1/2) {
  PIT <- rPIT(n, df, gamma)
  ifelse(PIT<delta,
         1-PIT/delta,
         (PIT-delta)/(1-delta))
}

# To verify by simulation
# dlta <- 0.81
# pSim <- rPITdelta(50000,5,gma,dlta) |> ecdf()
# plot(u, pPITdelta(u,5,gma,dlta), type='l')
# lines(u,pSim(u),col=2)


plotpPIT <- function(degfr, xlims=c(0,1), ylims=c(0,1), xanno=0.82, yanno=0.1) {
  u <- u[u>=xlims[1] & u<=xlims[2]]
  dat <- bind_rows(
    tibble(u=u,cdf=pPIT(u,degfr,gamma=gma), delta="identity"),
    tibble(u=u,cdf=pPITdelta(u,degfr,gamma=gma,1/2), delta="1/2"),
    tibble(u=u,cdf=pPITdelta(u,degfr,gamma=gma,2/3), delta="2/3")
  ) |>
    mutate(delta=as.factor(delta))
  if (is.infinite(degfr)) {
    dflab <- "df==infinity"
  } else {
    dflab <- paste0("df==",as.character(degfr))
  }
  p1 <- ggplot(dat, aes(x=u,y=cdf,color=delta, group=delta)) +
    geom_line() + xlim(xlims) + ylim(ylims) +
    geom_abline(intercept=0, slope=1, linetype=3) +
    theme_bw() +
    labs(x='v-transformed PIT', y='cdf', color=expression(delta))  +
    scale_colour_manual(values=c("blue","red","green")) +
    theme(legend.position = "bottom", legend.direction="horizontal",
          legend.title=element_text(size=legfont),
          legend.text=element_text(size=legfont)) +
    annotate('text', x=xanno,y=yanno, label=dflab, parse=TRUE)
             #parse=is.infinite(degfr))
  p1
}

pfull <- plotpPIT(df)
ptail <- plotpPIT(df, c(0.95,1),c(0.92,1),0.99,0.925)
# pN <- plotpPIT(Inf)
# p5 <- plotpPIT(5)
# 
# # ggsave mishandles the gamma in legend
# grDevices::cairo_pdf("output/PITcdf.pdf", width=7, height=3.5)
# pN+p5
# dev.off()
