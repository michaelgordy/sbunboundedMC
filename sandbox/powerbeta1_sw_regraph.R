# Size and Power of Unbounded monovariate continuous Kernels
# Michael Gordy

library(simsalapar)
library(spectralBacktest)
#library(abind)
library(future)
#library(magrittr)
library(dplyr)
library(tidyr)
library(moments)
library(purrr)
library(stringi)
library(ggplot2)
library(viridis)
library(gt)
#library(xtable)
num_cores <- as.integer(parallelly::availableCores(omit=8))
table_location <- "output/"
sim_location <- "simdata/"

load(file=paste0(sim_location, "power1grid.RData"))

gdf <- select(resdf, rejectrate, a, b) |> 
  mutate(bfct = as.factor(b))
ggplot(gdf, aes(x=a,y=rejectrate,group=bfct, color=bfct)) + geom_line() +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10() +
  theme_bw() +
  labs(y="Rejection Rate (log scale)", color="b")
ggsave(paste0(table_location, "powerbeta1.pdf"), width =6.5, height=3, units = "in")




