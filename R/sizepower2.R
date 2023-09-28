# Size and Power of Unbounded Bivariate Kernels
# Assumes current directory is the project folder
# Michael Gordy

source("R/simSetup.R")
n_days <- 500
nsims <- 2^8
blk_size <- nsims/2^5
savedata <- FALSE  # TRUE to have simsalapar save simulation data
gtsavename <- 'sizepower_beta2'

betakerns <- list(
  ZLL = list(name="(2,1)<br>(1,2)", param=list(c(2,1),c(1,2))),
  ZPP = list(name="(25,1)<br>(1,25)", param=list(c(25,1),c(1,25))),
  ZP4h = list(name="(2,0)<br>(1,3)", param=list(c(2,0),c(1,3))),
  ZP5h = list(name="(5/2,0)<br>(1/2,3)", param=list(c(2.5,0),c(0.5,3))),
  ZP9h = list(name="(9/2,0)<br>(1/2,6)", param=list(c(4.5,0),c(0.5,6)))
)

kern_vec2 <- c("ZLL", "ZPP", "ZP4h", "ZP5h", "ZP9h")
# kern_vec2 <- c("ZLL","ZPP","ZP0","PNS","LLS","GS","GcS")
F_names <- c("Normal", "Scaled t10", "Scaled t5", "Scaled t3")

doOne <- function(n,Fmodel,kernel){
  PIT <- choose_dist(Fmodel,n) |> pnorm() |> 
    pmin(1-.Machine$double.eps)
  purrr::map_dbl(kernel, ~spectral_Ztest(.x, PIT))
}

makeBetaKernel <- function(betaspec,support, type="multi") {
  param <- betaspec$param
  #type <- c('mono', 'bi', 'multi')[min(3,length(param))]
  if (length(betaspec$param)==1) 
    param <- unlist(param)
  list(
    name=betaspec$name,
    type=type,
    nu=nu_beta,
    correlation=rho_beta_beta,
    support=support,
    param=param )
}

# Run the simulation for a given window
rejectionrate <- function(bk_list,support,windowname,level=0.05) {
  kernelnames <- purrr::map_chr(bk_list,"name") 
  # kernelpstr <- purrr::map_chr(bk_list,"pstr") 
  # for (s in names(kernel_list)) assign(s, kernel_list[[s]])
  varList <-
    varlist( # constructor for an object of class 'varlist'
      n.sim = list(type="N", expr = quote(m), value = nsims),
      n = list(type="frozen", value = c(n_days)),
      Fmodel=list(type="grid", expr = quote(F), value = F_names),
      kernel = list(type="inner",  expr = quote(Kernel),
                    value=bk_list))
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
    pivot_longer(cols=all_of(F_names), names_to="Parameters:", values_to="rejectrate")
}

# rejectionrate <- function(support,windowname,level=0.05) {
#   kernel_list <- define_kernels(support[1], support[2], alpha_star=0.99)[kern_vec2]
#   kernelnames <- purrr::map_chr(kernel_list,"name") 
#   # for (s in names(kernel_list)) assign(s, kernel_list[[s]])
#   varList <-
#     varlist( # constructor for an object of class 'varlist'
#       n.sim = list(type="N", expr = quote(m), value = nsims),
#       n = list(type="frozen", value = c(n_days)),
#       Fmodel=list(type="grid", expr = quote(F), value = F_names),
#       kernel = list(type="inner",  expr = quote(Kernel),
#                     value=kernel_list))
#   if (savedata) {
#     savefilename <- paste0(sim_location, "continuous_unconditional_",
#                     windowname, "_", n_days, ".rds")
#   } else { savefilename <- NULL }
#   windowlabel <- paste("Window: ", windowname, " [",
#                         paste0(support,collapse=", "),"]")
#   doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
#                  doOne=doOne, block.size = blk_size, cores = num_cores,
#                  exports = ls(), extraPkgs = (.packages())) |>
#         getArray() |>
#         apply(c(1,2),function(p) mean(p<=level)) |>
#         as_tibble() |>
#         mutate(kernel=kernelnames, window=windowlabel, .before=1) |>
#         pivot_longer(cols=all_of(F_names), names_to="F", values_to="rejectrate")
# }

narrow <- purrr::map(betakerns[kern_vec2], 
                     ~makeBetaKernel(.,alpha_narrow)) |>
  rejectionrate(alpha_narrow,"narrow")
wide <- purrr::map(betakerns[kern_vec2], 
                   ~makeBetaKernel(.,alpha_wide)) |>
  rejectionrate(alpha_wide,"wide")
res <- dplyr::bind_rows(narrow,wide) |>
  pivot_wider(names_from = kernel, values_from = rejectrate) 

save(narrow,wide,res,alpha_wide,alpha_narrow,nsims,n_days,
     file=paste0(sim_location,gtsavename,".RData"))

tabL <- gt(res, groupname_col = "window") |>
  fmt_percent(columns=where(is.numeric), decimals=1) |>
  cols_label_with(fn=gt::md) |>
  tab_header(
    title='Size and power of tests based on Beta bikernels'
  ) |>
  tab_footnote(paste0('2^',log(nsims,2), ' trials with ', n_days,
               " observations per trial."))
gtfile<- paste0(table_location,gtsavename) 
gt::gtsave(tabL,filename = paste0(gtfile,".html"), inline_css=TRUE)
gt::gtsave(tabL,filename = paste0(gtfile,".tex"))

# Report on correlations
getcorrelation <- function(betakernel) {
  betakernel$correlation(betakernel$support,betakernel$param) 
}
correlationBeta <- function(support) {
  kernel_list <- define_kernels_beta2(support[1], support[2])[kern_vec2]
  purrr::map_dbl(kernel_list,getcorrelation)
}

