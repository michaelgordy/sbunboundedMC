# Recurrence rules based on:
# 2F1 = http://functions.wolfram.com/07.23.03.0003.01
# 3F2 = http://functions.wolfram.com/07.27.17.0001.01
# 4F3 = http://functions.wolfram.com/07.28.17.0001.01

# Generalized hypergeometric _{q+1}F_q with parameter vectors
# c(-k,rep(a,q)) and rep(a+1,q) with argument 1.
hyperFqk <- function(a,q,k,Fqklag=c()) {
  if (length(Fqklag)<q) {
    Fk <- sum((-1)^(0:k)* 
                exp(lchoose(k,0:k)+q*(log(a)-log(a+(0:k)))))
  } else {
    m <- length(Fqklag)+1
    if (q==1) {
      Fk <- Fqklag[m-1]*k/(k+a)
    } else if (q==2) {
      Fk <- ((2*(k+a)-1)*Fqklag[m-1]-(k-1)*Fqklag[m-2])*k/(k+a)^2
    } else if (q==3) {
      Fk <- -((3*(k+a)-3*(k+a)^2-1)*Fqklag[m-1]+
               3*(k-1)*(a+k-1)*Fqklag[m-2]-
               (k-1)*(k-2)*Fqklag[m-3])*k/(k+a)^3
    } else {
      stop('q must be leq 3')
    }
  }
  c(Fqklag, Fk) |> tail(q)
}

a<-2.7
q<-1
kmax<-25

# Compare direct and indirect calculation
yd <-purrr::map_dbl(0:kmax,~hyperFqk(a,q,.x)) # direct
yr <- rep(0,kmax+1)
Fqklag<-c()
for (k in 0:kmax) {
  Fqklag <- hyperFqk(a,q,k,Fqklag)
  yr[k+1] <- tail(Fqklag,1)
}

# test an idea
library(hypergeo)
theta <- 0.64
a <- 1/2+theta/pi
x <- 1.9
Rx <- pbeta(plogis(pi*x),a,1-a)
rho<-(pi/beta(a,1-a))*plogis(pi*x)^a*plogis(-pi*x)^(1-a)
lambda <- pi*(plogis(pi*x)-a)
lambdaprime <- pi^2*plogis(pi*x)*plogis(-pi*x)
B0 <- (pi^2/2)*a*(1-a)
B0x <- B0*Rx+pi^2/(2*beta(a,1-a))*plogis(pi*x)^a*plogis(-pi*x)^(2-a)


# Direct calculation of y[k]
# y<-purrr::map_dbl(0:kmax,~feefie(a,q,.x))
# r<-log(1-exp(diff(log(y))))

# k<-7
# a<-1.3
# a1<--k
# a2<-a3<-a
# b1<-b2<-a+1
# B1C1 <- (b1*b2+(a1+1)*(3*a1-2*b1-2*b2+4)+
#           (-a1+a2-1)*(a1-a3+1))/((a1-b1+1)*(a1-b2+1))

# Case q=2
# B1C1s <- k*(2*(k+a)-1)/(k+a)^2

# Indirect, starting with first two terms
# ynew <- rep(0,kmax+1)
# ynew[1:2] <- y[1:2]
# for (k in 2:kmax) {
#   ynew[k+1] <- (k*(2*(k+a)-1)*ynew[k]-k*(k-1)*ynew[k-1])/(k+a)^2
# }
