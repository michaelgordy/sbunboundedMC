# Size and Power of Unbounded monovariate continuous Kernels
# Michael Gordy

library(simsalapar)
library(spectralBacktest)
#library(abind)
library(future)
#library(magrittr)
library(dplyr)
library(tidyr)
library(gt)
#library(xtable)
num_cores <- as.integer(parallelly::availableCores(omit=8))
table_location <- "tables/"
sim_location <- "simresults/"
alpha_narrow <- c(0.98,1)
alpha_wide <-  c(0.95,1)
#alpha_star <- 0.99
n_days <- 500
nsims <- 2^14
savedata <- FALSE  # TRUE to have simsalapar save simulation data
gtsavename <- 'sizepower_beta1'

### Ideally this works
source("helperFiles/DefineUtilities.R")
source("helperFiles/DefineKernels.R")
# source("helperFiles/DefineCVTs.R")
blk_size <- nsims/2^5

kern_vec2 <- c("ZU", "ZA", "ZE", "ZLp", "Z1Q", "Z1E", "ZZ1", "ZZ5") 
F_names <- c("Normal", "Scaled t5", "Scaled t3")

# doOne <- function(n,Fmodel,kernel){
#   PITs <- choose_dist(Fmodel,n) |> pnorm()
#   purrr::map_dbl(as.list(kernel), ~spectral_Ztest(eval(parse(text=.x)), PITs))
# }
doOne <- function(n,Fmodel,kernel){
  PIT <- choose_dist(Fmodel,n) |> pnorm() |> 
    pmin(1-.Machine$double.eps)
  purrr::map_dbl(kernel, ~spectral_Ztest(.x, PIT))
}


# Run the simulation for a given window
rejectionrate <- function(support,windowname,level=0.05) {
  kernel_list <- define_kernels_beta1(support[1], support[2])[kern_vec2]
  kernelnames <- purrr::map_chr(kernel_list,"name") 
  # for (s in names(kernel_list)) assign(s, kernel_list[[s]])
  varList <-
    varlist( # constructor for an object of class 'varlist'
      n.sim = list(type="N", expr = quote(m), value = nsims),
      n = list(type="frozen", value = c(n_days)),
      Fmodel=list(type="grid", expr = quote(F), value = F_names),
      kernel = list(type="inner",  expr = quote(Kernel),
                    value=kernel_list))
  if (savedata) {
    savefilename <- paste0(sim_location, "continuous_unconditional_",
                    windowname, "_", n_days, ".rds")
  } else { savefilename <- NULL }
  windowlabel <- paste("Window: ", windowname, " [",
                        paste0(support,collapse=", "),"]")
  doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
                 doOne=doOne, block.size = blk_size, cores = num_cores,
                 exports = ls(), extraPkgs = (.packages())) |>
        getArray() |>
        apply(c(1,2),function(p) mean(p<=level)) |>
        as_tibble() |>
        mutate(kernel=kernelnames, window=windowlabel, .before=1) |>
        pivot_longer(cols=all_of(F_names), names_to="F", values_to="rejectrate")
}

res <- dplyr::bind_rows(rejectionrate(alpha_narrow, "narrow"),
                        rejectionrate(alpha_wide,"wide")) |>
  pivot_wider(names_from = kernel, values_from = rejectrate) 

tabL <- gt(res, groupname_col = "window") |>
  fmt_percent(columns=where(is.numeric), decimals=1) |>
  tab_header(
    title='Size and power of unconditional tests',
    subtitle='Best-performing monokernels from JBF2020 plus Beta(a,0)'
  ) |>
  tab_footnote(paste0('2^',log(nsims,2), ' trials with ', n_days,
               " observations per trial."))
gtfile<- paste0(table_location,gtsavename) 
gt::gtsave(tabL,filename = paste0(gtfile,".html"), inline_css=TRUE)
gt::gtsave(tabL,filename = paste0(gtfile,".tex"))


