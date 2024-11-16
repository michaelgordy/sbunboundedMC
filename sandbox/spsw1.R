# Size and Power of Unbounded monovariate beta kernels for single window
# Assumes current directory is the project folder
# Michael Gordy

source("R/simSetup.R")
n_days <- 500
nsims <- 2^16
blk_size <- nsims/2^5
kernwindow <- c(0.975,1)
savedata <- FALSE  # TRUE to have simsalapar save simulation data
gtsavename <- 'spsw_beta1'

betakerns <- list(
  ZU = list(name='Uniform', param=list(c(1,1)), 
            pstr="(1,1)"),
  ZA = list(name='Arcsin', param=list(c(1/2,1/2)), pstr="(1/2,1/2)"),
  ZE = list(name='Epinechnikov', param=list(c(2,2)), pstr="(2,2)"),
  ZLp = list(name = 'LinearUp', param=list(c(2,1)), pstr="(2,1)"),
  ZQp = list(name = 'QuadraticUp', param=list(c(3,1)), pstr="(3,1)"),
  Z1Q = list(name = 'Beta(1,1/4)', param=list(c(1,1/4)), pstr="(1,1/4)"),
  Z1E = list(name = 'Beta(1,1/8)', param=list(c(1,1/8)), pstr="(1,1/8)"),
  Z1Z = list(name = 'Beta(1,0)', param=list(c(1,0)), pstr="(1,0)"),
  Z2Z = list(name = 'Beta(2,0)', param=list(c(2,0)), pstr="(2,0)"),
  Z5Z = list(name = 'Beta(5,0)', param=list(c(5,0)), pstr="(5,0)")
)

kern_vec2 <- c("ZU", "ZA", "ZE", "ZLp", "Z1Q", "Z1E", "Z1Z", "Z2Z", "Z5Z") 
F_names <- c("Normal", "Scaled t10", "Scaled t5", "Scaled t3")

makeBetaKernel <- function(betaspec, type=NULL) {
  param <- betaspec$param
  if (is.null(type))
    type <- c('mono', 'bi', 'multi')[min(3,length(param))]
  if (length(betaspec$param)==1) 
    param <- unlist(param)
  list(
    name=betaspec$name,
    pstr=betaspec$pstr,
    type=type,
    nu=nu_beta,
    correlation=rho_beta_beta,
    support=kernwindow,
    param=param )
}

# Run the simulation for a given window
rejectionrate <- function(bk_list,level=0.05) {
  kernelnames <- purrr::map_chr(bk_list,"name") 
  kernelpstr <- purrr::map_chr(bk_list,"pstr") 
  # for (s in names(kernel_list)) assign(s, kernel_list[[s]])
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
  # windowlabel <- paste("Window: ", windowname, " [",
  #                      paste0(support,collapse=", "),"]")
  doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
            doOne=doOne, block.size = blk_size, cores = num_cores,
            exports = ls(), extraPkgs = (.packages())) |>
    getArray() |>
    apply(c(1,2),function(p) mean(p<=level)) |>
    as_tibble() |>
    mutate(kernel=kernelnames, pstr=kernelpstr, .before=1) |>
    pivot_longer(cols=all_of(F_names), names_to="Parameters", values_to="rejectrate")
}

res <- purrr::map(betakerns[kern_vec2], makeBetaKernel) |>
  rejectionrate() |>
  mutate(notBeta=if_else(stringr::str_starts(kernel,"Beta")," ",kernel),
         collab=paste(notBeta,pstr, sep="<br>")) |>
  dplyr::select(-kernel,-pstr,-notBeta) |>
  pivot_wider(names_from = collab, values_from = rejectrate) 

save(res,kernwindow,nsims,n_days,
     file=paste0(sim_location,gtsavename,".RData"))

gttabl <- gt(res) |>
  fmt_percent(columns=where(is.numeric), decimals=1) |>
  cols_label_with(fn=gt::md) |>
  tab_header(
    title='Size and power of tests based on Beta monokernels'
  ) |>
  tab_footnote(paste0('2^',log(nsims,2), ' trials with ', n_days,
               " observations per trial.  Kernel window is [",
               paste0(kernwindow,collapse=", "),"].") )
gtfile<- paste0(table_location,gtsavename) 
gt::gtsave(gttabl,filename = paste0(gtfile,".html"), inline_css=TRUE)
gt::gtsave(gttabl,filename = paste0(gtfile,".tex"))


