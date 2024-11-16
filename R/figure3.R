# Tail behavior of the G2 kernel

source("R/simSetup.R")
library(ggplot2)
library(spectralBacktest)
library(purrr)
library(dplyr)
library(tidyr)

support<- c(0.95,1)

PIT <- 1-10^seq(log10(1-1.01*support[1]),log10(1-0.999999*support[2]),length.out=250)
kern_vec2 <- c("PNS","LLS", "GS","GcS")
kernel_list <- define_kernels_tlsf(support[1], support[2])[kern_vec2]
kernnames <- data.frame(mnemonic=kern_vec2) |>
  mutate(kernel=purrr::map_chr(mnemonic, ~kernel_list[[.x]]$name))

G2 <- function(tlsfkernel, PIT) {
  W_list <- tlsfkernel$nu(tlsfkernel$support, tlsfkernel$param)(PIT)
  W_list[[2]]
}
df3 <- purrr::map_dfc(kernel_list, ~G2(.x,PIT)) |>
  mutate(PIT=PIT) |>
  tidyr::pivot_longer(cols=-PIT, names_to = "mnemonic", values_to = "value") |>
  dplyr::left_join(kernnames, by="mnemonic") |>
  mutate(x=-log(1-PIT), 
         kernel=forcats::fct_recode(kernel, Probitnormal="Probitnormal score", 
                            Logistic="Logit-Logistic score",
                            Gumbel="Gumbel score", `Comp Gumbel`="Comp Gumbel score"))
lw <- 0.4  # linewidth
p <- ggplot(df3, aes(x, y=value, color=kernel, group=kernel)) + 
  geom_line() + ylim(0,35) +
  labs(x=expression(-ln(1-PIT)),y=expression(G[2](PIT)-mu[2])) +
  theme_bw() +
  geom_vline(xintercept=-log(c(0.01, 0.001,0.0001)), linetype=3, linewidth=lw, show.legend = FALSE) +
  annotate("text", x=-log(0.006), y=32, label="99%") + 
  annotate("text", x=-log(0.00054), y=32, label="99.9%") +
  annotate("text", x=-log(0.000043), y=32, label="99.99%") +
  geom_abline(intercept=-1,slope=1, linetype=4, linewidth=lw) +
  geom_abline(intercept=-6,slope=2, linetype=4, linewidth=lw) 
ggsave("output/figure3.pdf", plot=p, width=6.5, height=3.5, units="in")
