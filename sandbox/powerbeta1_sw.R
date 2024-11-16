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
table_location <- "output/"
sim_location <- "simdata/"
kernwindow <- c(0.975,1)
n_days <- 500
nsims <- 2^16
savedata <- FALSE  # TRUE to have simsalapar save simulation data
#gtsavename <- 'powerbeta1_sw'
blk_size <- nsims/2^5

kernelgrid <- expand.grid(a=seq(0.25,8,by=0.25), b=c(0, 0.5, 1, 2, 4, 8)) |>
  mutate(name=stri_rand_strings(n(), 8, '[A-Za-z]'), .before=1)

F_names <- c("Scaled t5")

makeBetaKernel <- function(name, a, b, support=kernwindow) {
  list(
    name=name,
    type='mono',
    nu=nu_beta,
    support=support,
    param=c(a,b) )
}

# Run the simulation for a given window
# Key function doOne() defined in simSetup.R
rejectionrate <- function(bk_list,level=0.05) {
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
    mutate(kernel=kernelgrid$name, .before=1) |>
    pivot_longer(cols=all_of(F_names), names_to="F", values_to="rejectrate")
}

resdf <-  purrr::pmap(kernelgrid, 
                        function(name,a,b) 
                          makeBetaKernel(name,a,b))  |>
  set_names(kernelgrid$name) |> 
  rejectionrate() |> 
  left_join(kernelgrid,by=join_by(kernel==name))
save(resdf, kernelgrid,n_days, nsims, kernwindow,
     file=paste0(sim_location, "power1grid.RData"))

gdf <- select(resdf, rejectrate, a, b) |> 
  mutate(bfct = as.factor(b))
ggplot(gdf, aes(x=a,y=rejectrate,group=bfct, color=bfct)) + geom_line() +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10() +
  theme_bw() +
  labs(y="Rejection Rate (log scale)", color="b")
ggsave(paste0(table_location, "powerbeta1.pdf"), width =6.5, height=3, units = "in")

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
ggsave(paste0(table_location, "betakurt1.pdf"))




