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

source("R/DefineFmodels.R")
source("R/DefineKernels.R")

# Parameters shared across simulations
num_cores <- as.integer(parallelly::availableCores(omit=8))
table_location <- "output/"
sim_location <- "simdata/"
alpha_narrow <- c(0.98,1)
alpha_wide <-  c(0.95,1)

# Standard "doOne" functions for simsalapar
# When there is no v-transform
doOne <- function(n,Fmodel,kernel){
  PIT <- Fmodel_list[[Fmodel]](n) |> pnorm() |> 
    pmin(1-.Machine$double.eps)
  purrr::map_dbl(kernel, ~spectral_Ztest(.x, PIT))
}
# With v-transform
doOneV <- function(n,Fmodel,Vmodel,kernel){
  PIT <- Fmodel_list[[Fmodel]](n) |> pnorm() |>
    vtransform_list[[Vmodel]]() |> pmin(1-.Machine$double.eps)
  purrr::map_dbl(kernel, ~spectral_Ztest(.x, PIT))
}

