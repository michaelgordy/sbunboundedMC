# Demonstrate insensitivity to F0

source("R/DefineVtransforms.R")
library(ggplot2)
library(tidyr)

delta <- 0.3
kappa <- 3/2
xi <- 2/3

pst3 <- purrr::partial(tscopula::pst, mu=0, df=3, sigma=sqrt(1/3))
qst3 <- purrr::partial(tscopula::qst, mu=0, df=3, sigma=sqrt(1/3))
pst8 <- purrr::partial(tscopula::pst, mu=0, df=8, sigma=sqrt(3/4))
qst8 <- purrr::partial(tscopula::qst, mu=0, df=8, sigma=sqrt(3/4))

df <- data.frame(
     u=seq(from = 0, to = 1, length.out = 500),
     Normal = V_F0(u,pnorm,qnorm,delta,kappa,xi),
     t8 = V_F0(u,pst8,qst8,delta,kappa,xi),
     t3 = V_F0(u,pst3,qst3,delta,kappa,xi),
     Laplace = V_F0(u,plaplace0,qlaplace0,delta,kappa,xi) ) |>
  tidyr::pivot_longer(cols=-u, names_to="F0", values_to = "V")

ggplot(df, aes(x=u,y=V,color=F0)) + geom_line() + theme_bw()
ggsave("robustF0.pdf")