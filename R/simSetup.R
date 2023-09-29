# Shared libraries and fixed parameter choices
library(future)
library(simsalapar)
library(spectralBacktest)
#library(abind)
#library(magrittr)
#library(rlang)
library(glue)
library(tibble)
library(dplyr)
library(tidyr)
library(gt)
options(gt.html_tag_check = FALSE)

stopifnot((packageVersion("spectralBacktest")>="0.5.4"))

source("R/DefineUtilities.R")
source("R/DefineKernels.R")

# Parameters shared across simulations
num_cores <- as.integer(parallelly::availableCores(omit=8))
table_location <- "output/"
sim_location <- "simdata/"
alpha_narrow <- c(0.98,1)
alpha_wide <-  c(0.95,1)
