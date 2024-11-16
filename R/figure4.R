# Plot examples of v-transforms
library(dplyr)
library(ggplot2)
library(patchwork)

V_exponential <- function(u, delta=1/2,kappa=1,xi=1) {
  Psi <- function(x) exp(-kappa*(-log(x))^xi)
  Psiinv <- function(x) exp(-(log(x)/(-kappa))^(1/xi))
  ifelse(u<delta,
         (1-u)-(1-delta)*Psi(u/delta),
         u-delta*Psiinv((1-u)/(1-delta)))
}

x <- seq(0,1,by=0.001)
dflinear <- list(
  tibble(x=x, v=V_exponential(x,1/3), delta='1/3'),
  tibble(x=x, v=V_exponential(x,1/2), delta='1/2'),
  tibble(x=x, v=V_exponential(x,2/3), delta='2/3')) |>
  bind_rows() |>
  mutate(delta=as.factor(delta))

dfnonlinear <- list(
  tibble(x=x, v=V_exponential(x,1/2,2), linear=FALSE, delta='1/2', kappa='2'),
  tibble(x=x, v=V_exponential(x,1/2,1/2), linear=FALSE, delta='1/2', kappa='1/2')) |>
  bind_rows() |>
  mutate(kappa=as.factor(kappa))

legfont <- 11  
p1 <- ggplot(dflinear, aes(x=x,y=v,group=delta,color=delta)) +
  geom_line() + xlim(c(0,1)) + ylim(c(0,1)) +
  theme_bw() + 
  labs(x='v', y="T(v)", color=expression(delta)) + 
  theme(legend.position = "bottom", legend.direction="horizontal",
        legend.title=element_text(size=legfont), legend.text=element_text(size=legfont)) 

#  annotate("text", x=0.95,y=0.05,label="theta==1",size=12,parse=TRUE)

p2 <- ggplot(dfnonlinear, aes(x=x,y=v,group=kappa,color=kappa)) +
  geom_line() + xlim(c(0,1)) + ylim(c(0,1)) +
  theme_bw() + 
  labs(x='v', y="T(v)", color=expression(kappa)) + 
  theme(legend.position = "bottom", legend.direction="horizontal",
        legend.title=element_text(size=legfont), legend.text=element_text(size=legfont)) 

p1+p2
ggsave("output/figure4.pdf", width=6.5, height=3.5, units="in")


