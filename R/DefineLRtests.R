define_lr_tests <- function(alpha1, alpha2, alpha_star) {
        lr_test_list = list(
        LR1 = list( name = 'Binomial LRT at 99%',
                     pars = 0.99),
        
        LR2 = list( name = "Binomial LRT at 99%",
                    pars = c(alpha1, alpha2)),
        

        LR3 = list( name = 'Binomial LRT at 99%',
                     pars = c(alpha1, alpha_star, alpha2)),

        LRB = list( name = "Berkowitz",
                     pars = c(alpha1, alpha2))
        )
        return(lr_test_list)
}

spectral_LRtest <- function(test,PITs){
  PITs <- as.numeric(PITs)
  if (test$name=="Berkowitz")
    output <- Berkowitz_test(PITs,test$pars)
  else
    output <- mn_test(PITs,test$pars)
  output
}

# Multinomial LR-test

mn_test<- function(PITs,levels)
{
  mlevels <- matrix(levels,
                    ncol=length(levels),
                    nrow=length(PITs),byrow=TRUE)
  breaches <- apply(PITs > mlevels,1,sum)
  obs <- rep(NA,length(levels)+1)
  for (j in 0:length(levels))
    obs[j+1] <- sum(breaches==j)
  probs <- diff(c(0,levels,1))
  n <- sum(obs)
  stat  <- -2*(dmultinom(obs,n,probs,log=TRUE)-
                 dmultinom(obs,n,obs/n,log=TRUE))
  df <- length(obs)-1
  1-pchisq(stat,df)
}

# Berkowitz LR Test

Berkowitz_test <- function(PITs,interval){
  alpha1 <- interval[1]
  alpha2 <- interval[2]
  PITs <- as.numeric(PITs)
  PITs <- pmax(pmin(PITs,1-0.1^10),0.1^10) # truncate away from [0,1] regardless of alpha values
  zval <- qnorm(PITs)
  lTP <- qnorm(alpha1) # lower truncation point
  uTP <- qnorm(alpha2) # upper truncation point
  par.init <- c(mean(zval),sd(zval))
  par.null <- c(0,1)
  ll.null <- 	BK_loglik(par.null,zval,lTP,uTP)
  optim.init <- suppressWarnings(try(optim(par.init,BK_loglik,control=list(fnscale=-1),zval=zval,lTP=lTP,uTP=uTP),silent=T))
  optim.null <- suppressWarnings(try(optim(par.null,BK_loglik,control=list(fnscale=-1),zval=zval,lTP=lTP,uTP=uTP),silent=T))
  if (is(optim.init,"try-error"))
    ll.opt.init <- -Inf
  else{
    ll.opt.init <- optim.init$value
    par.opt <- optim.init$par
  }
  if (is(optim.null,"try-error"))
    ll.opt.null <- -Inf
  else{
    ll.opt.null <- optim.null$value
    if(ll.opt.null>ll.opt.init)
      par.opt <- optim.null$par
  }
  ll.opt <- max(ll.opt.init,ll.opt.null)
  test.statistic <- 2*(ll.opt-ll.null)
  df <- length(par.null)
  if (test.statistic >= 0){
    p.value <- 1-pchisq(test.statistic,df=df)
  }else{
    p.value <- NA
  }
  return(p.value)
}

BK_loglik <- function(par,zval,lTP,uTP){
  mu <- par[1]
  sigma <- par[2]
  ll.u <- which(zval >= uTP)
  ll.m <- which((zval < uTP) & (zval >lTP))
  ll.l <- which(zval <= lTP)
  ll <- sum(dnorm(zval[ll.m], mean = mu, sd = sigma, log=TRUE))
  ll <- ll + sum(pnorm(pmax(zval[ll.l],lTP), mean = mu, sd = sigma, lower.tail=TRUE,log.p=TRUE))
  ll <- ll + sum(pnorm(pmin(zval[ll.u],uTP), mean = mu, sd = sigma, lower.tail=FALSE,log.p=TRUE))
  return(ll)
}

