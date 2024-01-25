# Size and Power of Unbounded Bivariate Kernels for single window
# Assumes current directory is the project folder
# Michael Gordy

source("R/simSetup.R")
n_days <- 500
nsims <- 2^16
blk_size <- nsims/2^5
kernwindow <- c(0.975,1)
savedata <- FALSE  # TRUE to have simsalapar save simulation data
gtsavename <- 'spsw_beta2'

betakerns <- list(
  ZLL = list(name="(2,1)<br>(1,2)", param=list(c(2,1),c(1,2))),
  ZPP = list(name="(25,1)<br>(1,25)", param=list(c(25,1),c(1,25))),
  ZP4h = list(name="(2,0)<br>(1,3)", param=list(c(2,0),c(1,3))),
  ZP5h = list(name="(5/2,0)<br>(1/2,3)", param=list(c(2.5,0),c(0.5,3))),
  ZP9h = list(name="(9/2,0)<br>(1/2,6)", param=list(c(4.5,0),c(0.5,6)))
)

kern_vec2 <- c("ZLL", "ZPP", "ZP4h", "ZP5h", "ZP9h")
F_names <- c("Normal", "Scaled t10", "Scaled t5", "Scaled t3")

makeBetaKernel <- function(betaspec, type="multi") {
  param <- betaspec$param
  #type <- c('mono', 'bi', 'multi')[min(3,length(param))]
  if (length(betaspec$param)==1) 
    param <- unlist(param)
  list(
    name=betaspec$name,
    type=type,
    nu=nu_beta,
    correlation=rho_beta_beta,
    support=kernwindow,
    param=param )
}

# Run the simulation for a given window
rejectionrate <- function(bk_list,level=0.05) {
  kernelnames <- purrr::map_chr(bk_list,"name") 
  varList <-
    varlist( # constructor for an object of class 'varlist'
      n.sim = list(type="N", expr = quote(m), value = nsims),
      n = list(type="frozen", value = c(n_days)),
      Fmodel=list(type="grid", expr = quote(F), value = F_names),
      kernel = list(type="inner",  expr = quote(Kernel),
                    value=bk_list))
  if (savedata) {
    savefilename <- paste0(sim_location, "continuous_unconditional_sw_",
                           n_days, ".rds")
  } else { savefilename <- NULL }
  doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
            doOne=doOne, block.size = blk_size, cores = num_cores,
            exports = ls(), extraPkgs = (.packages())) |>
    getArray() |>
    apply(c(1,2),function(p) mean(p<=level)) |>
    as_tibble() |>
    mutate(kernel=kernelnames, .before=1) |>
    pivot_longer(cols=all_of(F_names), names_to="Parameters:", values_to="rejectrate")
}

res <- purrr::map(betakerns[kern_vec2], makeBetaKernel) |>
  rejectionrate() |>
  pivot_wider(names_from = kernel, values_from = rejectrate) 

save(res,kernwindow,nsims,n_days,
     file=paste0(sim_location,gtsavename,".RData"))

gttabl <- gt(res) |>
  fmt_percent(columns=where(is.numeric), decimals=1) |>
  cols_label_with(fn=gt::md) |>
  tab_header(
    title='Size and power of tests based on Beta bikernels'
  ) |>
  tab_footnote(paste0('2^',log(nsims,2), ' trials with ', n_days,
                      " observations per trial.  Kernel window is [",
                      paste0(kernwindow,collapse=", "),"].") )
gtfile<- paste0(table_location,gtsavename) 
gt::gtsave(gttabl,filename = paste0(gtfile,".html"), inline_css=TRUE)
gt::gtsave(gttabl,filename = paste0(gtfile,".tex"))

# Report on correlations
getcorrelation <- function(betakernel) {
  betakernel$correlation(betakernel$support,betakernel$param) 
}
correlationBeta <- 
  define_kernels_beta2(kernwindow[1], kernwindow[2])[kern_vec2] |>
  purrr::map_dbl(getcorrelation)


