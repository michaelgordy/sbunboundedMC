# Implement Psi(x) = exp(-kappa*(-log(x))^xi)
V_exponential <- function(u, delta=1/2,kappa=1,xi=1) {
  Psi <- function(x) exp(-kappa*(-log(x))^xi)
  Psiinv <- function(x) exp(-(log(x)/(-kappa))^(1/xi))
  ifelse(u<delta,
         (1-u)-(1-delta)*Psi(u/delta),
         u-delta*Psiinv((1-u)/(1-delta)))
}

vtransform_param <- tibble::tribble(
  ~name,~delta, ~kappa, ~xi,
  "identity", 0, 1, 1,
  "|1-2u|", 1/2, 1, 1,
  "(1/3, 1)", 1/3, 1, 1,
  "(2/3, 1)", 2/3, 1, 1,
  "(1/2, 1/2)", 1/2, 1/2, 1, 
  "(1/2, 2)", 1/2, 2, 1 )

# Convert to namd list of functions
vtransform_list <- purrr::pmap(vtransform_param,
                   function(name, delta, kappa, xi)
                           (\(u) V_exponential(u, delta, kappa, xi))) |>
  setNames(vtransform_param$name)
