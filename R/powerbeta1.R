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
savedata <- FALSE  # TRUE to have simsalapar save simulation data
gtsavename <- 'sizepower_beta1new'

### Ideally this works
source("helperFiles/DefineUtilities.R")
source("helperFiles/DefineKernels.R")
# source("helperFiles/DefineCVTs.R")
blk_size <- nsims/2^5

kernelgrid <- expand.grid(a=seq(0.25,8,by=0.25), b=c(0, 0.5, 1, 2, 4, 8)) |>
  mutate(name=stri_rand_strings(n(), 8, '[A-Za-z]'), .before=1)

F_names <- c("Scaled t5")

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
rejectionrate <- function(bk_list,support,windowname,level=0.05) {
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
  doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
            doOne=doOne, block.size = blk_size, cores = num_cores,
            exports = ls(), extraPkgs = (.packages())) |>
    getArray() |>
    apply(c(1,2),function(p) mean(p<=level)) |>
    as_tibble() |>
    mutate(kernel=kernelgrid$name, window=windowname, .before=1) |>
    pivot_longer(cols=all_of(F_names), names_to="F", values_to="rejectrate")
}

# narrowdf <- purrr::pmap(kernelgrid, 
#                         function(name,a,b) 
#                           makeBetaKernel(name,a,b,support=alpha_narrow)) |>
#   set_names(kernelgrid$name) |> 
#   rejectionrate(alpha_narrow,"narrow")
widedf <-  purrr::pmap(kernelgrid, 
                        function(name,a,b) 
                          makeBetaKernel(name,a,b,support=alpha_wide))  |>
  set_names(kernelgrid$name) |> 
  rejectionrate(alpha_wide,"wide")
#resdf <- dplyr::bind_rows(narrowdf, widedf) |>
#  left_join(kernelgrid,by=join_by(kernel==name))
resdf <-  left_join(widedf,kernelgrid,by=join_by(kernel==name))
save(resdf, kernelgrid,n_days, nsims, alpha_narrow, alpha_wide,
     file="data/power1grid.RData")

gdf <- filter(resdf,window=="wide", F=="Scaled t5") |>
  mutate(bfct = as.factor(b))
ggplot(gdf, aes(x=a,y=rejectrate,group=bfct, color=bfct)) + geom_line() +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10() +
  labs(y="Rejection Rate (log scale)", color="b")
ggsave("tables/powerbeta1.pdf")

U <- runif(100000)
betakurt <- kernelgrid |> 
  mutate(v=purrr::map2_dbl(a,b,
                           ~mean(nu_beta(param=c(.x,.y))(U)^2)),
    skw=purrr::map2_dbl(a,b,
                             ~mean(nu_beta(param=c(.x,.y))(U)^3)),
      krt=purrr::map2_dbl(a,b,
                             ~mean(nu_beta(param=c(.x,.y))(U)^4)),
      bfct = as.factor(b))
ggplot(betakurt, aes(x=a,y=krt,group=bfct, color=bfct)) + geom_line() +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10() + 
  labs(y="Kurtosis (log scale)", color="b")
ggsave("tables/betakurt1.pdf")



