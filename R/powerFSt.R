# Power when varying delta and kappa of v-transform
# Michael Gordy

library(stringi)
source("R/simSetup.R")
source("R/DefineVtransforms.R")
n_days <- 500
nsims <- 2^14
blk_size <- nsims/2^5
num_cores <- as.integer(parallelly::availableCores(omit=8))
savedata <- TRUE  # TRUE to have simsalapar save simulation data
gtsavename <- "power_deltakappa"

betakerns <- list(
  ZU = list(name='Beta(1,1) [Uniform]', param=list(c(1,1)), pstr="(1,1)"),
  ZA = list(name='Arcsin', param=list(c(1/2,1/2)), pstr="(1/2,1/2)"),
  ZLp = list(name = 'Beta(2,1) [LinearUp]', param=list(c(2,1)), pstr="(2,1)"),
  Z1Z = list(name = 'Beta(1,0)', param=list(c(1,0)), pstr="(1,0)"),
  Z2Z = list(name = 'Beta(2,0)', param=list(c(2,0)), pstr="(2,0)"),
  Z5Z = list(name = 'Beta(5,0)', param=list(c(5,0)), pstr="(5,0)")
)

F_names <- "FS-t(10,6/5)"
v_F0_param <- expand.grid(delta=seq(0.025,0.976,by=0.05), kappa=seq(0.5,2.5,by=0.5)) |>
  mutate(name=stri_rand_strings(n(), 8, '[A-Za-z]'), .before=1) |>
  mutate(xi=1)

# Convert to namd list of functions
Vlaplace <- purrr::partial(V_F0, pF0=plaplace0, qF0=qlaplace0)
vlaplace_list <- purrr::pmap(v_F0_param,
                             function(name, delta, kappa, xi)
                               (\(u) Vlaplace(u, delta, kappa, xi))) |>
  setNames(v_F0_param$name)
vtransform_list <- vlaplace_list
V_names <- names(vtransform_list)

kern_vec2 <- c("ZU", "ZLp", "Z1Z") 
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
      Fmodel=list(type="frozen", expr = quote(F), value = F_names),
      Vmodel=list(type="grid", expr = quote(V), 
                  value = names(vlaplace_list)),
      kernel = list(type="inner",  expr = quote(Kernel),
                    value=bk_list))
  if (savedata) {
    savefilename <- paste0(sim_location, gtsavename,
                           windowname, "_", n_days, ".rds")
  } else { savefilename <- NULL }
  doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
            doOne=doOneV, block.size = blk_size, cores = num_cores,
            exports = ls(), extraPkgs = (.packages())) |>
    getArray() |>
    apply(c(1,2),function(p) mean(p<=level)) |>
    as_tibble() |>
    mutate(kernel=kernelnames, .before=1) |> 
    pivot_longer(cols=c(-kernel), names_to="name", values_to="rejectrate") |>
    left_join(v_F0_param, by="name") |>
    rename(Vtransform=name)
}

# Make separate tables by window
# Each table is grouped by kernel
makePowerTable <- function(support, windowname) {
  gttabl <- purrr::map(betakerns[kern_vec2], 
                       ~makeBetaKernel(.,support)) |>
    rejectionrate(support,windowname)  |>
    mutate(kernel=paste0("Kernel: ", kernel)) |> 
    select(kernel,delta, kappa, rejectrate) |> 
    gt(groupname_col = "kernel") |>
    fmt_percent(columns=c(rejectrate), decimals=1) |>
    cols_label_with(fn=gt::md) |>
    tab_header(
      title=glue::glue('Power of v-transformed tests for DGP {F_names} ({windowname} window)')
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



# save(resdf, kernelgrid,n_days, nsims, alpha_narrow, alpha_wide,
#      file="data/power1grid.RData")

# gdf <- filter(resdf,window=="wide", F=="Scaled t5") |>
#   mutate(bfct = as.factor(b))
# ggplot(gdf, aes(x=a,y=rejectrate,group=bfct, color=bfct)) + geom_line() +
#   scale_color_viridis(discrete = TRUE) +
#   scale_y_log10() +
#   labs(y="Rejection Rate (log scale)", color="b")
# ggsave("tables/powerbeta1.pdf")


