library(tscopula)

# distributions F0 all constrained to have variance 1

vt_t <- function(u,delta = 0.45, kappa = 1.5, xi = 0.5, df =4){
  output <- rep(NA,length(u))
  uless <- u[u <= delta]
  umore <- u[u > delta]
  inner1 <- -qst(uless/(2*delta), df= df, mu = 0, sigma = sqrt((df-2)/df))
  inner2 <- -qst((1-umore)/(2*(1-delta)), df=df, mu = 0, sigma = sqrt((df-2)/df))
  output[u <= delta] <- (1-uless) - 
    (1-delta) *2 * pst(-kappa * inner1^xi, df=df, mu = 0, sigma = sqrt((df-2)/df))
  output[u > delta] <- umore - 
    delta * 2 * pst(-(inner2/kappa)^(1/xi), df=df, mu = 0, sigma = sqrt((df-2)/df))
  output
}


vt_laplace <- function(u,delta = 0.45, kappa = 1.5, xi = 0.5){
  output <- rep(NA,length(u))
  uless <- u[u <= delta]
  umore <- u[u > delta]
  inner1 <- -qlaplace0(uless/(2*delta))
  inner2 <- -qlaplace0((1-umore)/(2*(1-delta)))
  output[u <= delta] <- (1-uless) - 
    (1-delta) *2 * plaplace0(-kappa * inner1^xi)
  output[u > delta] <- umore - 
    delta * 2 * plaplace0(-(inner2/kappa)^(1/xi))
  output
}

vt_norm <- function(u,delta = 0.45, kappa = 1.5, xi = 0.5){
  output <- rep(NA,length(u))
  uless <- u[u <= delta]
  umore <- u[u > delta]
  inner1 <- -qnorm(uless/(2*delta))
  inner2 <- -qnorm((1-umore)/(2*(1-delta)))
  output[u <= delta] <- (1-uless) - 
    (1-delta) *2 * pnorm(-kappa * inner1^xi)
  output[u > delta] <- umore - 
    delta * 2 * pnorm(-(inner2/kappa)^(1/xi))
  output
}

delta <- 0.45
kappa <- 3/2
xi <- 2/3

u <- seq(from = 0, to = 1, length.out = 500)
plot(u, vt_t(u, delta=delta, kappa = kappa, xi = xi, df = 3), type="l", ylab = "T(u)")
lines(u, vt_laplace(u, delta = delta, kappa = kappa, xi = xi), col=2)
lines(u, vt_norm(u, delta = delta, kappa = kappa, xi =xi), col=3)
lines(u, vt_t(u, delta = delta, kappa = kappa, xi =xi, df = 8), col=4)


# Quick check on linear case

delta <- 0.45
kappa <- 1
xi <- 1

u <- seq(from = 0, to = 1, length = 500)
plot(u, vt_t(u, delta=delta, kappa = kappa, xi = xi, df = 3), type="l", ylab = "T(u)")
lines(u, vt_laplace(u, delta = delta, kappa = kappa, xi = xi), col=2)
lines(u, vt_norm(u, delta = delta, kappa = kappa, xi =xi), col=3)
lines(u, vt_t(u, delta = delta, kappa = kappa, xi =xi, df = 8), col=4)
