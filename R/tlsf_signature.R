library(ggplot2)

tlsfsig <- function(tlsfkernel, PIT) {
  W_list <- tlsfkernel$nu(tlsfkernel$support, tlsfkernel$param)(PIT)
  Wbar <- purrr::list_transpose(W_list, simplify=TRUE)
  Sigma<-tlsfkernel$VCV(tlsfkernel$support, tlsfkernel$param)
  purrr::map_dbl(Wbar, ~mahalanobis(.x, center=0, cov=Sigma))
}

support<-alpha_wide
PIT <- seq(1.0001*support[1],0.9999*support[2],length.out=100)
#kern_vec2 <- c("LLS", "LB1", "LB1c", "LB2", "LB2c", "LB3", "LB3c",
#               "LBtest", "LBtestc")
kern_vec2 <- c("LLS", 
               "LBtest", "LBtestc")
kernel_list <- define_kernels_tlsf(support[1], support[2])[kern_vec2]
zzsig <- purrr::map_dfc(kernel_list, ~tlsfsig(.x,PIT)) |>
  dplyr::mutate(alpha=PIT) |> 
  tidyr::pivot_longer(cols=-alpha, names_to = "kernel", values_to = "value")

p<- ggplot(zzsig, aes(x=alpha, y=value, color=kernel, group=kernel)) + 
  geom_line() +
  scale_y_log10() +
  labs(x="alpha",y="Mahalanobis")

