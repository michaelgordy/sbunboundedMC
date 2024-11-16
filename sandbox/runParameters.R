# This script sets and saves the run parameters that should be held constant across analyses.

minPlausible <- 1e-5  # lower threshold for Plausibility
sampleyears <- 3L
strtsamp <- c("2013-01-01", "2014-01-01", "2015-01-01") %>% 
              as.Date(origin="1970-01-01")

getEnddate <- function(dtstr,n) {
  yy <- format(dtstr,"%Y") %>% as.integer()
  endyr <- yy+n
  paste(endyr, format(dtstr,"-%m-%d"), sep='') %>%
    as.Date(origin="1970-01-01")
}

# Specify sample period for second exercise
daterangedf <- data_frame(
  startdate = strtsamp,
  enddate = getEnddate(strtsamp, sampleyears)
)
daterange <- daterangedf[2,] %>% unlist() %>% as.Date(origin="1970-01-01")

windowlist <- list(
  narrow=list(kernel_window="narrow", alpha1=0.985, alpha_star=0.99, alpha2=0.995),
  wide=list(kernel_window="wide",   alpha1=0.95,  alpha_star=0.99, alpha2=0.995)
)

createKernels <- function(wp) {
  
  BIN <- list( name = 'Binomial score at 99%',
               type = 'mono',
               nu = nu_discrete,
               support = wp$alpha_star,
               param = 1 )
  
  ZU3 <- list( name = 'Discrete Uniform 3',
               type = 'mono',
               nu = nu_discrete,
               support = c(wp$alpha1, wp$alpha_star, wp$alpha2),
               param = c(1, 1, 1) )
  
  ZE <- list( name = 'Epanechnikov',
              type = 'mono',
              nu = nu_epanechnikov,
              support = c(wp$alpha1, wp$alpha2),
              param = NULL )
  
  ZA <- list( name = 'Arcsin',
              type = 'mono',
              nu = nu_arcsin,
              support = c(wp$alpha1, wp$alpha2),
              param = NULL )
  
  ZU <- list( name = 'Uniform',
              type = 'mono',
              nu = nu_uniform,
              support = c(wp$alpha1, wp$alpha2),
              param = NULL )
  
  ZLp <- list( name = 'Linear increasing',
               type = 'mono',
               nu = nu_linear,
               support = c(wp$alpha1, wp$alpha2),
               param = 1 )
  
  ZLn <- list( name = 'Linear decreasing',
               type = 'mono',
               nu = nu_linear,
               support = c(wp$alpha1, wp$alpha2),
               param = -1 )
  
  ZLL <- list( name = 'linear/Linear',
               type = 'bi',
               nu = list(nu_linear, nu_linear),
               correlation = rho_linear_linear,
               support = c(wp$alpha1, wp$alpha2),
               param = list(-1,1) )
  
  ZAE <- list( name = 'Arcsin/Epanechnikov',
               type = 'bi',
               nu = list(nu_arcsin, nu_epanechnikov),
               correlation = rho_arcsin_epanechnikov,
               support = c(wp$alpha1, wp$alpha2),
               param = list(NULL, NULL) )
  
  ZPP <- list( name = 'bipower25',
               type = 'bi',
               nu = list(nu_beta, nu_beta),
               correlation = rho_beta_beta,
               support = c(wp$alpha1, wp$alpha2),
               param = list(c(25,1), c(1,25)) )
  
  PNS <- list( name = 'Probitnormal score',
               type = 'bi',
               nu = list(nu_probitnormal, nu_probitnormal),
               support = c(wp$alpha1, wp$alpha2),
               correlation = rho_probitnormal,
               param = list(1L, 2L) )

  Pearson2 <- list(name = 'Pearson2',
                   type = 'multi',
                   nu = nu_pearson,
                   correlation = rho_pearson,
                   support=NULL,
                   param=list(wp$alpha1, wp$alpha2))
  
  Pearson3 <- list(name = 'Pearson3',
                   type = 'multi',
                   nu = nu_pearson,
                   correlation = rho_pearson,
                   support=NULL,
                   param=list(wp$alpha1, wp$alpha_star, wp$alpha2))
  
  # gather the tests into a list and return
  list(BIN=BIN, ZU3=ZU3, ZU=ZU, ZLp=ZLp, ZLn=ZLn, 
       Pearson2=Pearson2, ZLL=ZLL, ZPP=ZPP,  Pearson3=Pearson3)
}

save(minPlausible, daterange, windowlist, createKernels,
     file='./data/runParameters.Rdata')
