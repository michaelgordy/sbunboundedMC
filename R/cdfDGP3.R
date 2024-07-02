# Script to plot tails under three DGP

library(tidyr)
library(dplyr)
library(forcats)
library(ggplot2)
library(RColorBrewer)

Fscaledt <- function(x, nu)
  pt(x*sqrt(nu/(nu-2)), df=nu)


PIT <- seq(0.9459, 0.9999, by=0.0005)
Fnormal <- PIT
Ft10 <- Fscaledt(qnorm(PIT), 10)
Ft5 <- Fscaledt(qnorm(PIT), 5)
Ft3 <- Fscaledt(qnorm(PIT), 3)
df <- data.frame(PIT=PIT, Normal=Fnormal, t10=Ft10, t5=Ft5, t3=Ft3) |>
  pivot_longer(cols=c(Normal,t10, t5, t3), names_to = 'DGP', values_to = "CDF") |>
  transform(DGP=factor(DGP,levels=c('Normal','t10','t5','t3'))) |>
  mutate(DGP = fct_recode(DGP, "Scaled t3"="t3", "Scaled t5"="t5", "Scaled t10"="t10"))

#
# "Pr(P(t)<=PIT)"
# p <- ggplot(df, aes(x=PIT, y=CDF, colour=DGP)) + geom_line() +
#   labs(y=expression(Pr(P[t] <= "PIT")), x="PIT", colour="True model") +
#   theme(legend.position = c(0.95, 0.05),
#         legend.justification = c("right", "bottom")) +
#   scale_x_continuous(breaks=c(0.95, 0.985, 0.995, 1.0),1)

p <- ggplot(df, aes(x=PIT, y=CDF, colour=DGP, linetype=DGP)) + geom_line() +
  scale_x_continuous(breaks=c(0.95, 0.975, 0.99, 1.0)) +
  scale_color_brewer(palette="Dark2") +
#  scale_linetype_manual(name = NULL, values = c("solid", "longdash", "dotdash")) +
  theme_bw() +
  theme(legend.title=element_text(size=13),
        legend.text=element_text(size=11),
        legend.position = c(0.95, 0.05),
        legend.justification = c("right", "bottom"),
        legend.key.width = unit(2.2,"line"),
        panel.border = element_rect(colour = 'black', size=1),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13)) +
  labs(y=expression(Pr(P[t] <= u)), x=expression(u), colour="True model", linetype="True model")

ggsave(filename='cdfDGP3.pdf', plot=p)