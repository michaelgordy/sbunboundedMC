# Tails of distribution functions for v-transformed PIT-values.

library(ggplot2)
library(grDevices)
library(patchwork)
library(latex2exp)
library(dplyr)
library(skewt)
source("R/DefineFmodels.R")

df <- 5
gmak <- 25
gma <- (gmak+1)/gmak
legfont <- 11

du <- 0.00025
u <- seq(0,1-du/2,by=du)

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


plotpPIT <- function(degfr, xlims=c(0,1), ylims=c(0,1), xanno=0.82, yanno=0.1, drophalf=FALSE) {
  u <- u[u>=xlims[1] & u<=xlims[2]]
  dat <- bind_rows(
    tibble(u=u,cdf=pPIT(u,degfr,gamma=gma), delta="identity"),
    tibble(u=u,cdf=pPITdelta(u,degfr,gamma=gma,1/3), delta="1/3"),
    tibble(u=u,cdf=pPITdelta(u,degfr,gamma=gma,1/2), delta="1/2"),
    tibble(u=u,cdf=pPITdelta(u,degfr,gamma=gma,2/3), delta="2/3")
  ) |>
    mutate(delta=factor(delta, levels=c('identity','1/3','1/2','2/3')))
  linecol <- c("blue","purple","red","green")
  if (drophalf) {
    linecol <- linecol[-3]
    dat <- filter(dat, delta!="1/2")
  }
  if (is.infinite(degfr)) {
    dflab <- as.character(expression(paste0(zeta,"=",infinity))) #\U03B6\U003D\U221E"))
  } else {
    dflab <- TeX(paste0("$\\zeta=",as.character(degfr),
                        ",\ \\gamma=", glue::glue("{gmak+1}/{gmak}"),"$"))
  }
#  
  p1 <- ggplot(dat, aes(x=u,y=cdf,color=delta, group=delta)) +
    geom_line() + xlim(xlims) + ylim(ylims) +
    geom_abline(intercept=0, slope=1, linetype=3) +
    theme_bw() +
    labs(y=expression(Pr(T(P[t]) <= u)), x=expression(u), color=expression(delta))  +
    scale_colour_manual(values=linecol) +
    theme(legend.position = "bottom", legend.direction="horizontal",
          legend.title=element_text(size=legfont),
          legend.text=element_text(size=legfont)) +
    annotate('text', x=xanno,y=yanno, label=dflab, parse=TRUE)
             #parse=is.infinite(degfr))
  p1
}

#pfull <- plotpPIT(df)
ptail <- plotpPIT(df, c(0.95,1),c(0.945,1),0.9928,0.948, drophalf=TRUE)
# pN <- plotpPIT(Inf)
# p5 <- plotpPIT(5)
# 
# # ggsave mishandles the gamma in legend
grDevices::cairo_pdf("output/figure6.pdf", width=7, height=3.5)
ptail
dev.off()
