# Size and Power of Unbounded Bivariate Kernels
# Assumes current directory is the project folder
# Michael Gordy

source("R/simSetup.R")
n_days <- 500
kernwindow <- c(0.95,1)
nsims <- 2^16
blk_size <- nsims/2^5
savedata <- FALSE  # TRUE to have simsalapar save simulation data
gtsavename <- 'spsw_tlsf'


kern_vec2 <- c("PNS", "GcS","GS", "LLS", "LB2","LB3c")
F_names <- c("Normal", "Scaled t10", "Scaled t5", "Scaled t3")

# Run the simulation for a given window
rejectionrate <- function(level=0.05) {
  kernel_list <- define_kernels_tlsf(kernwindow[1], kernwindow[2])[kern_vec2]
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
    savefilename <- paste0(sim_location, "tlsf_sw_",n_days, ".rds")
  } else { savefilename <- NULL }
  doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
                 doOne=doOne, block.size = blk_size, cores = num_cores,
                 exports = ls(), extraPkgs = (.packages())) |>
        getArray() |>
        apply(c(1,2),function(p) mean(p<=level)) |>
        as_tibble() |>
        mutate(kernel=kernelnames, .before=1) |>
        pivot_longer(cols=all_of(F_names), names_to="F", values_to="rejectrate")
}

# Correlation table
getcorrelation <- function(tlsfkernel) {
  corkern <- tlsfkernel$VCV(tlsfkernel$support,tlsfkernel$param) |>
    cov2cor() |> (\(x) x[1,2])()
}
correlationTLSF <-  define_kernels_tlsf(kernwindow[1], kernwindow[2])[kern_vec2] |>
  purrr::map_dbl(getcorrelation)

res <- rejectionrate() |>
  pivot_wider(names_from = kernel, values_from = rejectrate) 

save(res,kernwindow,nsims,n_days, correlationTLSF,
     file=paste0(sim_location,gtsavename,".RData"))

gttabl <- gt(res) |>
  fmt_percent(columns=where(is.numeric), decimals=1) |>
  tab_header(
    title='Size and power of tests based on TLSF kernels'
  ) |>
  tab_footnote(paste0('2^',log(nsims,2), ' trials with ', n_days,
                      " observations per trial.  Kernel window is [",
                      paste0(kernwindow,collapse=", "),"].") )
gtfile<- paste0(table_location,gtsavename) 
gt::gtsave(gttabl,filename = paste0(gtfile,".html"), inline_css=TRUE)
gt::gtsave(gttabl,filename = paste0(gtfile,".tex"))


