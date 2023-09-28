library(simsalapar)
library(spectralBacktest)
#library(abind)
library(future)
#library(magrittr)
library(dplyr)
library(tidyr)
library(gt)
options(gt.html_tag_check = FALSE)

### Ideally this works
source("R/DefineUtilities.R")
source("R/DefineKernels.R")

num_cores <- as.integer(parallelly::availableCores(omit=8))
table_location <- "output/"
sim_location <- "simdata/"
alpha_narrow <- c(0.98,1)
alpha_wide <-  c(0.95,1)
#alpha_star <- 0.99
