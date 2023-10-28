# skew-normal

dsnorm <- function (x, gamma, mu = 0, sigma = 1, log = FALSE) 
{
  result <- rep(NA, length(x))
  x <- (x - mu)/sigma
  result[x < 0] <- dnorm(gamma * x[x < 0], log = log)
  result[x >= 0] <- dnorm(x[x >= 0]/gamma, log = log)
  if (log) {
    return(result + log(2/(gamma + 1/gamma)) - log(sigma))
  }
  else {
    return(result * (2/(gamma + 1/gamma))/sigma)
  }
}

psnorm <- function (q, gamma, mu = 0, sigma = 1) 
{
  result <- rep(NA, length(q))
  x <- (q - mu)/sigma
  result[x < 0] <- 2/(gamma^2 + 1) * pnorm(gamma * x[x < 0])
  result[x >= 0] <- 1/(gamma^2 + 1) + 2/(1 + (1/gamma^2)) * 
    (pnorm(x[x >= 0]/gamma) - 1/2)
  result
}

qsnorm <- function (p, gamma, mu = 0, sigma = 1) 
{
  result <- rep(NA, length(p))
  probzero <- 1/(gamma^2 + 1)
  result[p < probzero] <- 1/gamma * qnorm(((gamma^2 + 1) * p[p < 
                                                            probzero])/2)
  result[p >= probzero] <- gamma * qnorm((1 + 1/gamma^2)/2 * (p[p >= 
                                                               probzero] - probzero) + 1/2)
  result * sigma + mu
}

rsnorm <- function (n, gamma, mu = 0, sigma = 1) 
{
  qsnorm(runif(n), gamma, mu, sigma)
}

momentsnorm <- function(gamma, mu = 0, sigma = 1){
  M1 <- 2*dnorm(0)
  dist.mn <- M1*(gamma-1/gamma)
  dist.var <- (1-M1^2)*((gamma^2)+1/(gamma^2)) + 2*(M1^2) - 1
  c(dist.mn*sigma + mu, dist.var*sigma^2)			 
}

# Some checks

integrate(dsnorm, lower = -Inf, upper = Inf, gamma = 2, mu = -0.5, sigma = 1.1)
integrate(dsnorm, lower = -Inf, upper = 3, gamma = 2, mu = -0.5, sigma = 1.1)
psnorm(3, gamma = 2, mu = -0.5, sigma = 1.1)
qsnorm(0.910696, gamma = 2, mu = -0.5, sigma = 1.1)
integrate(dsnorm, lower = -Inf, upper = -1, gamma = 2, mu = -0.5, sigma = 1.1)
psnorm(-1, gamma = 2, mu = -0.5, sigma = 1.1)
qsnorm(0.07266043, gamma = 2, mu = -0.5, sigma = 1.1)
integrate(function(x,gamma,mu,sigma){x*dsnorm(x,gamma,mu,sigma)}, 
          lower = -Inf, upper = Inf, gamma = 2, mu = -0.5, sigma = 1.1)
integrate(function(x,gamma,mu,sigma){(x^2)*dsnorm(x,gamma,mu,sigma)}, 
          lower = -Inf, upper = Inf, gamma = 2, mu = -0.5, sigma = 1.1)
mn <- momentsnorm(gamma = 2, mu = -0.5, sigma = 1.1)[1]
vr <- momentsnorm(gamma = 2, mu = -0.5, sigma = 1.1)[2]
c(mn, vr + mn^2)
mean(rsnorm(10000, gamma = 2, mu = -0.5, sigma = 1.1))

# skew t

dsst <- function (x, df, gamma, mu = 0, sigma = 1, log = FALSE) 
{
  result <- rep(NA, length(x))
  x <- (x - mu)/sigma
  result[x < 0] <- dt(gamma * x[x < 0], df, log = log)
  result[x >= 0] <- dt(x[x >= 0]/gamma, df, log = log)
  if (log) {
    return(result + log(2/(gamma + 1/gamma)) - log(sigma))
  }
  else {
    return(result * (2/(gamma + 1/gamma))/sigma)
  }
}

psst <- function (q, df, gamma, mu = 0, sigma = 1) 
{
  result <- rep(NA, length(q))
  x <- (q - mu)/sigma
  result[x < 0] <- 2/(gamma^2 + 1) * pt(gamma * x[x < 0], df)
  result[x >= 0] <- 1/(gamma^2 + 1) + 2/(1 + (1/gamma^2)) * 
    (pt(x[x >= 0]/gamma, df) - 1/2)
  result
}

qsst <- function (p, df, gamma, mu = 0, sigma = 1) 
{
  result <- rep(NA, length(p))
  probzero <- 1/(gamma^2 + 1)
  result[p < probzero] <- 1/gamma * qt(((gamma^2 + 1) * p[p < 
                                                            probzero])/2, df)
  result[p >= probzero] <- gamma * qt((1 + 1/gamma^2)/2 * (p[p >= 
                                                               probzero] - probzero) + 1/2, df)
  result * sigma + mu
}

rsst <- function (n, df, gamma, mu = 0, sigma = 1) 
{
  qsst(runif(n), df, gamma, mu, sigma)
}

momentsst <- function(df, gamma, mu = 0, sigma = 1){
  M1 <- dt(0,df=df)*df/(df-1)*2 
  dist.mn <- M1*(gamma-1/gamma)
  M2 <- df/(df-2)
  dist.var <- (M2-M1^2)*((gamma^2)+1/(gamma^2)) + 2*(M1^2) - M2
  c(dist.mn*sigma + mu, dist.var*sigma^2)			 
}

# Some checks

integrate(dsst, lower = -Inf, upper = Inf, df = 4, gamma = 2, mu = -0.5, sigma = 1.1)
integrate(dsst, lower = -Inf, upper = 3, df = 4, gamma = 2, mu = -0.5, sigma = 1.1)
psst(3, df = 4, gamma = 2, mu = -0.5, sigma = 1.1)
qsst(0.8505283, df = 4, gamma = 2, mu = -0.5, sigma = 1.1)
integrate(dsst, lower = -Inf, upper = -1, df = 4, gamma = 2, mu = -0.5, sigma = 1.1)
psst(-1, df = 4, gamma = 2, mu = -0.5, sigma = 1.1)
qsst(0.08294478, df = 4, gamma = 2, mu = -0.5, sigma = 1.1)
integrate(function(x,df,gamma,mu,sigma){x*dsst(x,df,gamma,mu,sigma)}, 
          lower = -Inf, upper = Inf, df = 4, gamma = 2, mu = -0.5, sigma = 1.1)
integrate(function(x,df,gamma,mu,sigma){(x^2)*dsst(x,df,gamma,mu,sigma)}, 
          lower = -Inf, upper = Inf, df = 4, gamma = 2, mu = -0.5, sigma = 1.1)
mn <- momentsst(df = 4, gamma = 2, mu = -0.5, sigma = 1.1)[1]
vr <- momentsst(df = 4, gamma = 2, mu = -0.5, sigma = 1.1)[2]
c(mn, vr + mn^2)
mean(rsst(10000, df = 4, gamma = 2, mu = -0.5, sigma = 1.1))
