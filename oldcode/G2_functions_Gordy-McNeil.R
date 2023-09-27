library(qrmtools)

alpha1 <- 0.95 
k <- 1:5
u <- 1-10^(-k)
x <- -log(1-u)

G2_norm <- qnorm(alpha1)*dnorm(qnorm(alpha1))/alpha1 + qnorm(u)^2
G2_Gumbel <- qGEV(alpha1, 0)*dGEV(qGEV(alpha1, 0), 0)/alpha1 + qGEV(u, 0)*(1+log(u))
G2_cGumbel <- -qGEV(1-alpha1, 0)*dGEV(qGEV(1-alpha1, 0), 0)/alpha1 + qGEV(1-u, 0)*(1+log(1-u))
G2_logistic <- qlogis(alpha1)*dlogis(qlogis(alpha1))/alpha1 + qlogis(u)*(2*u-1)

plot(x,G2_norm, ylim = range(G2_norm, G2_Gumbel, G2_cGumbel, G2_logistic), type = "l", 
     ylab = expression(G[2](u)), xlab = "-ln(1-u)", col = 6)
lines(x, G2_logistic, col=3)
lines(x, G2_Gumbel, col=4)
lines(x, G2_cGumbel, col=5)

# Linear growth functions
abline(0,1, lty =2, col =1)
abline(0,2, lty =2, col =1)

# 99th and 99.5th percentiles
abline(v = -log(0.01), lty=3, col=2)
abline(v = -log(0.001), lty =3, col =2)

