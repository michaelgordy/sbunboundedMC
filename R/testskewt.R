# Standardizing the skew-t distribution
# Page 12 http://www.eief.it/files/2008/11/monti-8-maggio.pdf
# Page 17 http://azzalini.stat.unipd.it/SN/se-ext.pdf

library(sn)

# Note sign error in Monti-8-maggio. Azzalini correct
rsst <- function(n=1,alpha=0,nu=Inf) {
  delta <- alpha/sqrt(1+alpha^2)
  bnu <- sqrt(nu/pi)*exp(lgamma(nu/2-1/2)-lgamma(nu/2))
  omega <- 1/sqrt(nu/(nu-2)-(delta*bnu)^2)
  xi <- -omega*delta*bnu
  sn::rst(n=n,xi=xi,omega=omega,alpha=alpha,nu=nu)
}
