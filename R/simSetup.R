library(simsalapar)
library(spectralBacktest)
library(abind)
#library(here)
library(future)
library(magrittr)
library(xtable)

### Ideally this works
source("helperFiles/DefineUtilities.R")
source("helperFiles/DefineKernels.R")
source("helperFiles/DefineCVTs.R")
source("helperFiles/DefineLRtests.R")


table_location <- "tables/"
sim_location <- "simresults/"
alpha_narrow <- c(0.985,1)
alpha_wide <- c(0.95,1)
alpha_star <- 0.99

narrow_name <- paste("[",paste(alpha_narrow,collapse=", "),"]",sep="")
wide_name <- paste("[",paste(alpha_wide,collapse=", "),"]",sep="")

n_days <- 750
nsims <- 2^12
blk_size <- nsims/2^5

# Setting core number and blocksize

num_cores <- as.integer(parallelly::availableCores(omit=8))
