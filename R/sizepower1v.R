# Size and Power of Unbounded monovariate beta kernels
# Assumes current directory is the project folder
# Michael Gordy

source("R/simSetup.R")
source("R/DefineVtransforms.R")
n_days <- 500
nsims <- 2^12
blk_size <- nsims/2^5
savedata <- FALSE  # TRUE to have simsalapar save simulation data
gtsavename <- 'sizepower_beta1v'

betakerns <- list(
  ZU = list(name='Uniform', param=list(c(1,1)), pstr="(1,1)"),
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

kern_vec2 <- c("ZU", "ZLp", "Z1Z", "Z2Z") 
F_names <- c("Normal", "Scaled t10", "Scaled t5", "Scaled t3")
V_names <- names(vtransform_list)

doOne <- function(n,Fmodel,Vmodel,kernel){
  PIT <- choose_dist(Fmodel,n) |> pnorm() |>
    vtransform_list[[Vmodel]]() |> pmin(1-.Machine$double.eps)
  purrr::map_dbl(kernel, ~spectral_Ztest(.x, PIT))
}

makeBetaKernel <- function(betaspec,support, type=NULL) {
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
      Vmodel=list(type="grid", expr = quote(V), 
                  value = names(vtransform_list)),
      kernel = list(type="inner",  expr = quote(Kernel),
                    value=bk_list))
  if (savedata) {
    savefilename <- paste0(sim_location, "beta1v_",
                           windowname, "_", n_days, ".rds")
  } else { savefilename <- NULL }
  doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
            doOne=doOne, block.size = blk_size, cores = num_cores,
            exports = ls(), extraPkgs = (.packages())) |>
    getArray() |>
    apply(c(1,2,3),function(p) mean(p<=level)) |>
    as_tibble() |>
    mutate(kernel=kernelnames, .before=1) |> 
    pivot_longer(cols=c(-kernel), names_to="gridvar", values_to="rejectrate") |>
    separate_wider_delim(gridvar, delim='.', names=c("F","Vtransform"))
}

# Make separate tables by window
# Each table is grouped by kernel
makePowerTable <- function(support, windowname) {
 gttabl <- purrr::map(betakerns[kern_vec2], 
                    ~makeBetaKernel(.,support)) |>
    rejectionrate(support,windowname) |>
    pivot_wider(names_from = Vtransform, values_from = rejectrate) |> 
    mutate(kernel=paste0("Kernel: ", kernel)) |> 
    gt(groupname_col = "kernel") |>
    fmt_percent(columns=where(is.numeric), decimals=1) |>
    cols_label_with(fn=gt::md) |>
    tab_header(
      title=glue::glue('Size and power of v-transformed tests ({windowname} window)')
    ) |>
    tab_footnote(paste0('2^',log(nsims,2), ' trials with ', n_days,
                        " observations per trial. Kernel window is [",
                        paste0(support,collapse=", "),"]"))
 gtfile<- paste0(table_location,gtsavename,"_",windowname) 
 gt::gtsave(gttabl,filename = paste0(gtfile,".html"), inline_css=TRUE)
 gt::gtsave(gttabl,filename = paste0(gtfile,".tex"))
 return(gttabl)
}

gtnarrow <- makePowerTable(alpha_narrow,"narrow")
gtwide <- makePowerTable(alpha_wide,"wide")
