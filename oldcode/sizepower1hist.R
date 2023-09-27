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
table_location <- "tables/"
sim_location <- "simresults/"
alpha_narrow <- c(0.98,1)
alpha_wide <-  c(0.95,1)
#alpha_star <- 0.99
n_days <- 500
nsims <- 2^16
savedata <- TRUE  # TRUE to have simsalapar save simulation data
gtsavename <- 'sizepower_beta1new'

### Ideally this works
source("helperFiles/DefineUtilities.R")
source("helperFiles/DefineKernels.R")
# source("helperFiles/DefineCVTs.R")
blk_size <- nsims/2^5

#kernelgrid <- expand.grid(a=c(0.25,0.5,1,2,4), b=c(0, 0.5, 1, 2, 4, 8)) |>
#  mutate(name=stri_rand_strings(n(), 8, '[A-Za-z]'), .before=1)

kernelgrid <- tibble::tribble(
  ~name, ~a, ~b,
  "(1,1)", 1, 1,
  "(2, 1/2)", 2, 0.5,
  "(4, 0)", 4, 0 
)
  
F_names <- c("Normal")

# doOne <- function(n,Fmodel,kernel){
#   PITs <- choose_dist(Fmodel,n) |> pnorm()
#   purrr::map_dbl(as.list(kernel), ~spectral_Ztest(eval(parse(text=.x)), PITs))
# }
doOne <- function(n,Fmodel,kernel){
  PIT <- choose_dist(Fmodel,n) |> pnorm() |> 
    pmin(1-.Machine$double.eps)
  purrr::map_dbl(kernel, ~spectral_Ztest(.x, PIT))
}

makeBetaKernel <- function(name, a, b, support) {
  list(
    name=name,
    type='mono',
    nu=nu_beta,
    support=support,
    param=c(a,b) )
}

# Run the simulation for a given window
runSim <- function(bk_list,support,windowname) {
  varList <-
    varlist( # constructor for an object of class 'varlist'
      n.sim = list(type="N", expr = quote(m), value = nsims),
      n = list(type="frozen", value = c(n_days)),
      Fmodel=list(type="grid", expr=quote(F), value = F_names),
      kernel = list(type="inner",  expr = quote(Kernel),
                    value=bk_list))
  if (savedata) {
    savefilename <- paste0(sim_location, "continuous_unconditional_",
                           windowname, "_", n_days, ".rds")
  } else { savefilename <- NULL }
  doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
            doOne=doOne, block.size = blk_size, cores = num_cores,
            exports = ls(), extraPkgs = (.packages())) 
}

res <- purrr::pmap(kernelgrid, 
                        function(name,a,b) 
                          makeBetaKernel(name,a,b,support=alpha_narrow)) |>
     set_names(kernelgrid$name) |>
     runSim(alpha_narrow,"narrow") 

narrowdf <- purrr::map_dfc(1:nrow(kernelgrid),
                         function(k)
                           purrr::map_dbl(res,~.x$value[k])) |>
  set_names(kernelgrid$name) |>
  tidyr::pivot_longer(cols=everything(), 
                      names_to = "Kernel", values_to = "pvalue") |>
  arrange(Kernel)

res <- purrr::pmap(kernelgrid, 
                   function(name,a,b) 
                     makeBetaKernel(name,a,b,support=alpha_wide)) |>
  set_names(kernelgrid$name) |>
  runSim(alpha_wide,"wide") 

widedf <- purrr::map_dfc(1:nrow(kernelgrid),
                         function(k)
                           purrr::map_dbl(res,~.x$value[k])) |>
  set_names(kernelgrid$name) |>
  tidyr::pivot_longer(cols=everything(), 
                      names_to = "Kernel", values_to = "pvalue") |>
  arrange(Kernel)

brks <- seq(0,1,by=0.05)
ggplot(widedf, aes(x=pvalue)) +
  geom_histogram(mapping=aes(y=after_stat(density)), breaks=brks)+
  facet_wrap(.~Kernel,nrow=nrow(kernelgrid)) +
  labs(x="p-value")

