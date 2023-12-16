
# Closures provide RNG for Scaled t and Scaled Skew-t 
rscaled_t <- function(nu) {
  \(n) sqrt((nu-2)/nu)*rt(n, df=nu)
}
rscaledskewed_t_sn <- function(nu, alpha=0) {
  delta <- alpha/sqrt(1+alpha^2)
  bnu <- sqrt(nu/pi)*exp(lgamma(nu/2-1/2)-lgamma(nu/2))
  omega <- 1/sqrt(nu/(nu-2)-(delta*bnu)^2)
  xi <- -omega*delta*bnu
  \(n) sn::rst(n,xi=xi,omega=omega,alpha=alpha,nu=nu)
}

# Standardized t in Fernandez-Steel class
qsst_fs <- function (p, df, gamma=1) 
{
  probzero <- 1/(gamma^2 + 1)
  suppressWarnings(
    q <- ifelse(p<probzero,
            (1/gamma) * qt(0.5*p/probzero, df), 
            gamma * qt(1/2+0.5*(p-probzero)/(gamma^2*probzero), df))
  )
  M1 <- dt(0,df=df)*df/(df-1)*2 
  dist.mn <- M1*(gamma-1/gamma)
  M2 <- df/(df-2)
  dist.var <- (M2-M1^2)*((gamma^2)+1/(gamma^2)) + 2*(M1^2) - M2
  (q-dist.mn)/sqrt(dist.var)  # standardized q
}

rsst_fs <- function (df, gamma=1) 
{
  \(n) qsst_fs(runif(n), df, gamma)
}

# Return named list with a pair of FS-t(m,k/(k+1)) and FS-t(m,(k+1)/k)
rsst_fs_list <- function(m,k) {
  list(rsst_fs(m,(k+1)/k), rsst_fs(m,k/(k+1))) |>
     setNames(c(glue::glue("FS-t({m},{k}/{k+1})"),
                glue::glue("FS-t({m},{k+1}/{k})")))
}
  
# Named list of available F models
Fmodel_list2 <- list(
  Normal = rnorm,
  "Scaled t10" = rscaled_t(10),
  "Scaled t5" = rscaled_t(5),
  "Scaled t3" = rscaled_t(3),
  "SS-t(10,1)" = rscaledskewed_t_sn(10,1),
  "SS-t(10,-1)" = rscaledskewed_t_sn(10,-1),
  "SS-t(5,1)" = rscaledskewed_t_sn(5,1),
  "SS-t(5,-1)" = rscaledskewed_t_sn(5,-1),
  "FS-t(10,3/2)" = rsst_fs(10,3/2),
  "FS-t(10,2/3)" = rsst_fs(10,2/3),
  "FS-t(5,3/2)" = rsst_fs(5,3/2),
  "FS-t(5,2/3)" = rsst_fs(5,2/3),
  "FS-t(3,3/2)" = rsst_fs(3,3/2),
  "FS-t(3,2/3)" = rsst_fs(3,2/3),
  "FS-t(10,6/5)" = rsst_fs(10,6/5),
  "FS-t(10,5/6)" = rsst_fs(10,5/6),
  "FS-t(5,6/5)" = rsst_fs(5,6/5),
  "FS-t(5,5/6)" = rsst_fs(5,5/6),
  "FS-t(3,6/5)" = rsst_fs(3,6/5),
  "FS-t(3,5/6)" = rsst_fs(3,5/6)
)

Fmodel_list <-  list(
  list(
  Normal = rnorm,
  "Scaled t10" = rscaled_t(10),
  "Scaled t5" = rscaled_t(5),
  "Scaled t3" = rscaled_t(3),
  "SS-t(10,1)" = rscaledskewed_t_sn(10,1),
  "SS-t(10,-1)" = rscaledskewed_t_sn(10,-1),
  "SS-t(5,1)" = rscaledskewed_t_sn(5,1),
  "SS-t(5,-1)" = rscaledskewed_t_sn(5,-1)),
 purrr::map(c(10,5,3), ~rsst_fs_list(.x,2)),
 purrr::map(c(10,5,3), ~rsst_fs_list(.x,3)),
 purrr::map(c(10,5,3), ~rsst_fs_list(.x,5))
 ) |> purrr::list_flatten() |> purrr::list_flatten() 

