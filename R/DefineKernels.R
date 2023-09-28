define_kernels_discrete <- function(alpha1, alpha2, alpha_star) {
  alpha1_mid <- (alpha1+alpha_star)/2
  alpha2_mid <- (alpha2+alpha_star)/2

  list(
    BIN = list( name = paste0("Binomial score at ", alpha_star*100, "%"),
                type = 'mono',
                nu = nu_discrete,
                support = alpha_star,
                param = 1 ),
    
    ZU2 = list( name = 'Discrete Uniform 2',
                type = 'mono',
                nu = nu_discrete,
                support = c(alpha1, alpha2),
                param = c(1, 1) ),
    
    ZU3 = list( name = 'Discrete Uniform 3',
                type = 'mono',
                nu = nu_discrete,
                support = c(alpha1, alpha_star, alpha2),
                param = c(1, 1, 1) ),
    
    ZU5 = list( name = 'Discrete Uniform 5',
                type = 'mono',
                nu = nu_discrete,
                support = c(alpha1, alpha1_mid, alpha_star, alpha2_mid, alpha2),
                param = c(1, 1, 1, 1, 1) ),
    
    ZL3 = list( name = 'Discrete Linear 3',
                type = 'mono',
                nu = nu_discrete,
                support = c(alpha1, alpha_star, alpha2),
                param = (1:3)/3 ),
    
    ZL5 = list( name = 'Discrete Linear 5',
                type = 'mono',
                nu = nu_discrete,
                support = c(alpha1, alpha1_mid, alpha_star, alpha2_mid, alpha2),
                param = (1:5)/5 ),
    
    
    PE2 = list(name = 'Pearson 2',
               type = 'multi',
               nu = nu_pearson,
               correlation = rho_pearson,
               support=NULL,
               param=list(alpha1, alpha2)),
    
    PE3 = list(name = 'Pearson 3',
               type = 'multi',
               nu = nu_pearson,
               correlation = rho_pearson,
               support=NULL,
               param=list(alpha1, alpha_star, alpha2)),
    
    PE5 = list(name = 'Pearson 5',
               type = 'multi',
               nu = nu_pearson,
               correlation = rho_pearson,
               support=NULL,
               param=list(alpha1, alpha1_mid, alpha_star, alpha2_mid, alpha2)) )
}

define_kernels_beta1 <- function(alpha1, alpha2) {

  list(
  ZU = list( name = 'Uniform',
               type = 'mono',
               nu = nu_uniform,
               support = c(alpha1, alpha2),
               param = NULL ),
    
  ZE = list( name = 'Epanechnikov',
              type = 'mono',
              nu = nu_epanechnikov,
              support = c(alpha1, alpha2),
              param = NULL ),

  ZA = list( name = 'Arcsin',
              type = 'mono',
              nu = nu_arcsin,
              support = c(alpha1, alpha2),
              param = NULL ),

  ZQ = list( name = 'Quadratic Increasing',
             type = 'mono',
             nu = nu_beta,
             support = c(alpha1, alpha2),
             param = c(3,1) ),

  Z1Q = list( name = 'Beta(1,1/4)',
             type = 'mono',
             nu = nu_beta,
             support = c(alpha1, alpha2),
             param = c(1,1/4) ),
  Z1E = list( name = 'Beta(1,1/8)',
              type = 'mono',
              nu = nu_beta,
              support = c(alpha1, alpha2),
              param = c(1,1/8) ),
  
  ZZ1 = list( name = 'Beta(1,0)',
             type = 'mono',
             nu = nu_beta,
             support = c(alpha1, alpha2),
             param = c(1,0) ),
  ZZ5 = list( name = 'Beta(5,0)',
              type = 'mono',
              nu = nu_beta,
              support = c(alpha1, alpha2),
              param = c(5,0) ),
  
  ZLp = list( name = 'LinearUp',
               type = 'mono',
               nu = nu_linear,
               support = c(alpha1, alpha2),
               param = 1 ),

  ZLn = list( name = 'LinearDown',
               type = 'mono',
               nu = nu_linear,
               support = c(alpha1, alpha2),
               param = -1 ) )
}

define_kernels_tlsf <- function(alpha1, alpha2) {
  
  list(
  PNS = list( name = 'Probitnormal score',
              type = 'tlsf',
              nu = nu_probitnormal,
              VCV = vcv_probitnormal,
              support = c(alpha1, alpha2),
              param = NULL ),
  
  LLS = list( name = 'Logit-Logistic score',
               type = 'tlsf',
               nu = nu_logitlogistic,
               VCV = vcv_logitlogistic,
               support = c(alpha1, alpha2),
               param = NULL ),
  
  GS = list( name = 'Gumbel score',
              type = 'tlsf',
              nu = nu_gumbel,
              VCV = vcv_gumbel,
              support = c(alpha1, alpha2),
              param = FALSE ),
  GcS = list( name = 'Comp Gumbel score',
               type = 'tlsf',
               nu = nu_gumbel,
               VCV = vcv_gumbel,
               support = c(alpha1, alpha2),
               param = TRUE),
  
  LB1 = list( name = 'Logistic Beta(1,1/2)',
              type = 'tlsf',
              nu = nu_logisticbeta,
              VCV = vcv_logisticbeta,
              support = c(alpha1, alpha2),
              param = c(1,1/2)),
  LB1c = list( name = 'Logistic Beta(1/2,1)',
              type = 'tlsf',
              nu = nu_logisticbeta,
              VCV = vcv_logisticbeta,
              support = c(alpha1, alpha2),
              param = c(1/2,1)),
  LB2 = list( name = 'Logistic Beta(3/2,1/2)',
              type = 'tlsf',
              nu = nu_logisticbeta,
              VCV = vcv_logisticbeta,
              support = c(alpha1, alpha2),
              param = c(3/2,1/2)),
  LB2c = list( name = 'Logistic Beta(1/2,3/2)',
              type = 'tlsf',
              nu = nu_logisticbeta,
              VCV = vcv_logisticbeta,
              support = c(alpha1, alpha2),
              param = c(1/2,3/2)),
  LB3 = list( name = 'Logistic Beta(2/3,1/3)',
              type = 'tlsf',
              nu = nu_logisticbeta,
              VCV = vcv_logisticbeta,
              support = c(alpha1, alpha2),
              param = c(2/3,1/3)),
  LB3c = list( name = 'Logistic Beta(1/3,2/3)',
               type = 'tlsf',
               nu = nu_logisticbeta,
               VCV = vcv_logisticbeta,
               support = c(alpha1, alpha2),
               param = c(1/3,2/3)),
  LLtest = list( name = 'Logistic Beta(1,1)',
              type = 'tlsf',
              nu = nu_logisticbeta,
              VCV = vcv_logisticbeta,
              support = c(alpha1, alpha2),
              param = c(1,1)),
  LBtest = list( name = 'Logistic Beta(1.02,0.98)',
                 type = 'tlsf',
                 nu = nu_logisticbeta,
                 VCV = vcv_logisticbeta,
                 support = c(alpha1, alpha2),
                 param = c(1.02,0.98)),
  LBtestc = list( name = 'Logistic Beta(0.98,1.02)',
               type = 'tlsf',
               nu = nu_logisticbeta,
               VCV = vcv_logisticbeta,
               support = c(alpha1, alpha2),
               param = c(0.98,1.02)),
  LBtest1 = list( name = 'Logistic Beta(1.02,1.02)',
                 type = 'tlsf',
                 nu = nu_logisticbeta,
                 VCV = vcv_logisticbeta,
                 support = c(alpha1, alpha2),
                 param = c(1.02,1.02)),
  LBtest2 = list( name = 'Logistic Beta(0.98,0.98)',
                  type = 'tlsf',
                  nu = nu_logisticbeta,
                  VCV = vcv_logisticbeta,
                  support = c(alpha1, alpha2),
                  param = c(0.98,0.98)),
  LBup = list( name = 'Logistic Beta(3/2,3/2)',
                 type = 'tlsf',
                 nu = nu_logisticbeta,
                 VCV = vcv_logisticbeta,
                 support = c(alpha1, alpha2),
                 param = c(3/2,3/2)),
  LBdown = list( name = 'Logistic Beta(1/2,1/2)',
                  type = 'tlsf',
                  nu = nu_logisticbeta,
                  VCV = vcv_logisticbeta,
                  support = c(alpha1, alpha2),
                  param = c(1/2,1/2))
  )
}

define_kernels_beta2 <- function(alpha1, alpha2) {
 list(
  ZUL = list( name = 'Uniform/Linear',
              type = 'bi',
              nu = list(nu_beta, nu_beta),
              correlation = rho_beta_beta,
              support = c(alpha1, alpha2),
              param = list(c(1,1),c(2,1)) ),
  
  ZLL = list( name = 'linear/Linear',
              type = 'bi',
              nu = list(nu_linear, nu_linear),
              correlation = rho_linear_linear,
              support = c(alpha1, alpha2),
              param = list(-1,1) ),
  
  ZAE = list( name = 'Arcsin/Epanechnikov',
              type = 'bi',
              nu = list(nu_arcsin, nu_epanechnikov),
              correlation = rho_arcsin_epanechnikov,
              support = c(alpha1, alpha2),
              param = list(NULL, NULL) ),
  
  ZQQ = list( name = 'quadratic/Quadratic',
              type = 'bi',
              nu = list(nu_beta, nu_beta),
              correlation = rho_beta_beta,
              support = c(alpha1, alpha2),
              param = list(c(3,1), c(1,3)) ),
  
  ZQUQ = list( name = 'quadratic/Uniform/Quadratic',
               type = 'multi',
               nu = nu_beta,
               correlation = rho_beta_beta,
               support = c(alpha1, alpha2),
               param = list(c(3,1), c(1,1), c(1,3)) ), 
  
  ZPP = list( name = 'beta(25,1)/(1,25)',
              type = 'bi',
              nu = list(nu_beta, nu_beta),
              correlation = rho_beta_beta,
              support = c(alpha1, alpha2),
              param = list(c(25,1), c(1,25)) ),

  ZP5h = list( name = 'beta (5/2,0)/(1/2,3)',
               type = 'bi',
               nu = list(nu_beta, nu_beta),
               correlation = rho_beta_beta,
               support = c(alpha1, alpha2),
               param = list(c(5/2,0), c(1/2,3)) ),
  
  ZP4h = list( name = 'beta (2,0)/(1,3)',
               type = 'bi',
               nu = list(nu_beta, nu_beta),
               correlation = rho_beta_beta,
               support = c(alpha1, alpha2),
               param = list(c(2,0), c(1,3)) ),
  
  ZP9h = list( name = 'beta (9/2,0)/(1/2,5)',
              type = 'bi',
              nu = list(nu_beta, nu_beta),
              correlation = rho_beta_beta,
              support = c(alpha1, alpha2),
              param = list(c(9/2,0), c(1/2,5)) ),
  
  ZPUP = list( name = 'quadratic/Uniform/Quadratic',
               type = 'multi',
               nu = nu_beta,
               correlation = rho_beta_beta,
               support = c(alpha1, alpha2),
               param = list(c(25,1), c(1,1), c(1,25)) ) )
}

define_kernels <- function(alpha1, alpha2, alpha_star) {
  c(define_kernels_discrete(alpha1,alpha2,alpha_star),
    define_kernels_beta1(alpha1,alpha2),
    define_kernels_beta2(alpha1,alpha2),
    define_kernels_tlsf(alpha1,alpha2)) 
}

newBetaKernel <- function(support=c(0,1), param=c(1,1), kernelname=NULL) {
  if (is.list(param)) {
    if (length(param)==2) {
      type <- 'bi'
    } else {
      type <- 'multi'
    }
  }
  list(
      name=kernelname,
      type=type,
      nu=nu_beta,
      correlation=rho_beta_beta,
      support=support,
      param=param )
}

# PP10 = list( name = 'quadratic/Quadratic',
#             type = 'bi',
#             nu = list(nu_beta, nu_beta),
#             correlation = rho_beta_beta,
#             support = c(alpha1, alpha2),
#             param = list(c(10,1), c(1,10)) ),
# 
# PP25 = list( name = 'quadratic/Quadratic',
#              type = 'bi',
#              nu = list(nu_beta, nu_beta),
#              correlation = rho_beta_beta,
#              support = c(alpha1, alpha2),
#              param = list(c(25,1), c(1,25)) ),
# 
# PP50 = list( name = 'quadratic/Quadratic',
#              type = 'bi',
#              nu = list(nu_beta, nu_beta),
#              correlation = rho_beta_beta,
#              support = c(alpha1, alpha2),
#              param = list(c(50,1), c(1,50)) ),
# 
# PUP25 = list( name = 'quadratic/Uniform/Quadratic',
#              type = 'multi',
#              nu = nu_beta,
#              correlation = rho_beta_beta,
#              support = c(alpha1, alpha2),
#              param = list(c(25,1), c(1,1), c(1,25)) ), 
# 
# PUP50 = list( name = 'quadratic/Uniform/Quadratic',
#               type = 'multi',
#               nu = nu_beta,
#               correlation = rho_beta_beta,
#               support = c(alpha1, alpha2),
#               param = list(c(50,1), c(1,1), c(1,50)) ), 
# 
# PUP100 = list( name = 'quadratic/Uniform/Quadratic',
#               type = 'multi',
#               nu = nu_beta,
#               correlation = rho_beta_beta,
#               support = c(alpha1, alpha2),
#               param = list(c(100,1), c(1,1), c(1,100)) ),

# Implement Psi(x) = exp(-kappa*(-log(x))^xi)
Vexplog <- function(u, delta=1/2,kappa=1,xi=1) {
  Psi <- function(x) exp(-kappa*(-log(x))^xi)
  Psiinv <- function(x) exp(-(log(x)/(-kappa))^(1/xi))
  ifelse(u<=delta,
         (1-u)-(1-delta)*Psi(u/delta),
         u-delta*Psiinv((1-u)/(1-delta)))
}

vtransform_list <- list(
  identity,
  (\(u) abs(1-2*u)),
  (\(u) Vexplog(u, delta=0.4)),
  (\(u) Vexplog(u, xi=2)) 
)
