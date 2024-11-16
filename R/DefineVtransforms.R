# Implement Psi(x) = exp(-kappa*(-log(x))^xi)
# V_exponential <- function(u, delta=1/2,kappa=1,xi=1) {
#   Psi <- function(x) exp(-kappa*(-log(x))^xi)
#   Psiinv <- function(x) exp(-(log(x)/(-kappa))^(1/xi))
#   ifelse(u<delta,
#          (1-u)-(1-delta)*Psi(u/delta),
#          u-delta*Psiinv((1-u)/(1-delta)))
# }
# 
# v_exponential_param <- tibble::tribble(
#   ~name,~delta, ~kappa, ~xi,
#   "identity", 0, 1, 1,
#   "|1-2u|", 1/2, 1, 1,
#   "(1/3, 1)", 1/3, 1, 1,
#   "(2/3, 1)", 2/3, 1, 1,
#   "(1/2, 1/2)", 1/2, 1/2, 1, 
#   "(1/2, 2)", 1/2, 2, 1 )
# 
# # Convert to namd list of functions
# vexponential_list <- purrr::pmap(v_exponential_param,
#                    function(name, delta, kappa, xi)
#                            (\(u) V_exponential(u, delta, kappa, xi))) |>
#   setNames(v_exponential_param$name)

plaplace0 <- function(q) {
  ifelse(q<=0, 
         0.5*exp(sqrt(2)*q),
         1-0.5*exp(-sqrt(2)*q))
} 
qlaplace0 <- function(p) {
  ifelse(p<=1/2,
         sqrt(0.5)*log(2*p),
         -sqrt(0.5)*log(2*(1-p)))
}

V_F0 <- function(u, pF0, qF0, delta = 0.5, kappa = 1, xi = 1) {
  suppressWarnings(
  nq <- ifelse(u<=delta,
             -qF0(0.5*u/delta),
             -qF0(0.5*(1-u)/(1-delta))))
  suppressWarnings(
  ifelse(u<=delta,
         (1-u) - (1-delta)*2*pF0(-kappa*nq^xi),
         u - delta*2*pF0(-(nq/kappa)^(1/xi))))
}

v_F0_param <- tibble::tribble(
  ~name,~delta, ~kappa, ~xi,
  "identity", 0, 1, 1,
  "|1-2u|", 1/2, 1, 1,
  "(1/3, 1)", 1/3, 1, 1,
  "(2/3, 1)", 2/3, 1, 1,
  "(1/2, 1/2)", 1/2, 1/2, 1, 
  "(1/2, 2)", 1/2, 2, 1 )

# Convert to namd list of functions
Vlaplace <- purrr::partial(V_F0, pF0=plaplace0, qF0=qlaplace0)
vlaplace_list <- purrr::pmap(v_F0_param,
                                 function(name, delta, kappa, xi)
                                   (\(u) Vlaplace(u, delta, kappa, xi))) |>
  setNames(v_F0_param$name)

#pF0_laplace <- purrr::partial(tscopula::plaplace,scale = sqrt(0.5))
#qF0_laplace <- purrr::partial(tscopula::qlaplace,scale = sqrt(0.5))
