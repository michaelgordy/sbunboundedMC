# Size and Power of Unbounded monovariate beta kernels
# Assumes current directory is the project folder
# Michael Gordy

source("R/simSetup.R")
source("R/DefineVtransforms.R")
n_days <- 500
nsims <- 2^16
level <- 0.05
blk_size <- nsims/2^5
savedata <- FALSE  # TRUE to have simsalapar save simulation data
gtsavename <- 'sizepower_beta1v_varywindow'

# Define two kernels, each with same beta parameters.
# First has standard window, second has wider window
# 
alpha1 <- 0.98
Zwide <- Znarrow <- list( name = 'Narrow',
                type = 'mono',
                nu = nu_uniform,
                support = c(alpha1, 1),
                param = NULL )
Zwide$support[1] <- 2*alpha1-1   # double the width
Zwide$name <- "Wide"
bk_list <- list(Znarrow, Zwide)
F_names <- c("Normal", "Scaled t10", "Scaled t5", "Scaled t3")
vtransform_list <- vlaplace_list[1:2]
V_names <- names(vtransform_list)

kernelnames <- purrr::map_chr(bk_list,"name") 
varList <-
    varlist( # constructor for an object of class 'varlist'
      n.sim = list(type="N", expr = quote(m), value = nsims),
      n = list(type="frozen", value = c(n_days)),
      Fmodel=list(type="grid", expr = quote(F), value = F_names),
      Vmodel=list(type="grid", expr = quote(V), 
                  value = names(vtransform_list)),
      kernel = list(type="inner",  expr = quote(Kernel),
                    value=bk_list))
if (savedata) {
    savefilename <- paste0(sim_location, "beta1v_",
                           windowname, "_", n_days, ".rds")
} else { savefilename <- NULL }
df <- doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
            doOne=doOneV, block.size = blk_size, cores = num_cores,
            exports = ls(), extraPkgs = (.packages())) |>
    getArray() |>
    apply(c(1,2,3),function(p) mean(p<=level)) |>
    as_tibble() |>
    mutate(kernel=kernelnames, .before=1) |> 
    pivot_longer(cols=c(-kernel), names_to="gridvar", values_to="rejectrate") |>
    separate_wider_delim(gridvar, delim='.', names=c("F","Vtransform"))

gttabl <-  df |>  
    pivot_wider(names_from = Vtransform, values_from = rejectrate) |> 
    mutate(kernel=paste0("Kernel: ", kernel)) |>
    gt(groupname_col = "kernel") |>
    fmt_percent(columns=where(is.numeric), decimals=1) |>
    cols_label_with(fn=gt::md) |>
    tab_header(
      title=glue::glue('Size and power of v-transformed tests')
    ) |>
    tab_footnote(paste0('2^',log(nsims,2), ' trials with ', n_days,
                        " observations per trial."))
gtfile<- paste0(table_location,gtsavename)
gt::gtsave(gttabl,filename = paste0(gtfile,".html"), inline_css=TRUE)
gt::gtsave(gttabl,filename = paste0(gtfile,".tex"))
  
# Make separate tables by window
# Each table is grouped by kernel
# makePowerTable <- function(support, windowname) {
#  gttabl <- purrr::map(betakerns[kern_vec2], 
#                     ~makeBetaKernel(.,support)) |>
#     rejectionrate(support,windowname) |>
#     pivot_wider(names_from = Vtransform, values_from = rejectrate) |> 
#     mutate(kernel=paste0("Kernel: ", kernel)) |> 
#     gt(groupname_col = "kernel") |>
#     fmt_percent(columns=where(is.numeric), decimals=1) |>
#     cols_label_with(fn=gt::md) |>
#     tab_header(
#       title=glue::glue('Size and power of v-transformed tests ({windowname} window)')
#     ) |>
#     tab_footnote(paste0('2^',log(nsims,2), ' trials with ', n_days,
#                         " observations per trial. Kernel window is [",
#                         paste0(support,collapse=", "),"]"))
#  gtfile<- paste0(table_location,gtsavename,"_",windowname) 
#  gt::gtsave(gttabl,filename = paste0(gtfile,".html"), inline_css=TRUE)
#  gt::gtsave(gttabl,filename = paste0(gtfile,".tex"))
#  return(gttabl)
# }

