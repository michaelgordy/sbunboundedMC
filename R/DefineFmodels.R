
# Closures provide RNG for Scaled t and Scaled Skew-t 
rscaled_t <- function(nu) {
  \(n) sqrt((nu-2)/nu)*rt(n, df=nu)
}
rscaledskewed_t <- function(nu, alpha=0) {
  delta <- alpha/sqrt(1+alpha^2)
  bnu <- sqrt(nu/pi)*exp(lgamma(nu/2-1/2)-lgamma(nu/2))
  omega <- 1/sqrt(nu/(nu-2)-(delta*bnu)^2)
  xi <- -omega*delta*bnu
  \(n) sn::rst(n,xi=xi,omega=omega,alpha=alpha,nu=nu)
}

# Named list of available F models
Fmodel_list <- list(
  Normal = rnorm,
  "Scaled t10" = rscaled_t(10),
  "Scaled t5" = rscaled_t(5),
  "Scaled t3" = rscaled_t(3),
  "SS-t(10,1)" = rscaledskewed_t(10,1),
  "SS-t(10,-1)" = rscaledskewed_t(10,-1),
  "SS-t(5,1)" = rscaledskewed_t(5,1),
  "SS-t(5,-1)" = rscaledskewed_t(5,-1)
)
