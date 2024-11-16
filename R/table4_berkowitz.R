# Berkowitz column in Table 4

source("R/simSetup.R")
n_days <- 500
nsims <- 2^16
blk_size <- nsims/2^5
kernwindow <- c(0.975,1)
level <- 0.05
savedata <- FALSE  # TRUE to have simsalapar save simulation data
gtsavename <- 'table4_berkowitz'

Berkowitz_test <- function(PITs,alpha1){
  PITs <- as.numeric(PITs)
  PITs <- pmax(pmin(PITs,1-0.1^10),0.1^10) # truncate away from [0,1] regardless of alpha values
  zval <- qnorm(PITs)
  lTP <- qnorm(alpha1) # lower truncation point
  par.init <- c(mean(zval),sd(zval))
  par.null <- c(0,1)
  ll.null <- 	BK_loglik(par.null,zval,lTP)
  optim.init <- suppressWarnings(try(optim(par.init,BK_loglik,control=list(fnscale=-1),zval=zval,lTP=lTP),silent=T))
  optim.null <- suppressWarnings(try(optim(par.null,BK_loglik,control=list(fnscale=-1),zval=zval,lTP=lTP),silent=T))
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

BK_loglik <- function(par,zval,lTP,uTP=Inf){
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

doBerkowitz <- function(n,Fmodel, alpha1){
  PIT <- Fmodel_list[[Fmodel]](n) |> pnorm() |> 
    pmin(1-.Machine$double.eps)
  Berkowitz_test(PIT, alpha1)
}

F_names <- c("Normal", "Scaled t10", "Scaled t5", "Scaled t3")

varList <- varlist( # constructor for an object of class 'varlist'
    n.sim = list(type="N", expr = quote(m), value = nsims),
    n = list(type="frozen", value = c(n_days)),
    Fmodel=list(type="grid", expr = quote(F), value = F_names),
    alpha1=list(type="grid", expr = quote(alf1), value = kernwindow[1]))
  
  # windowlabel <- paste("Window: ", windowname, " [",
  #                      paste0(support,collapse=", "),"]")
  # doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
  #           doOne=doOne, block.size = blk_size, cores = num_cores,
  #           exports = ls(), extraPkgs = (.packages())) |>
  #   getArray() |>
  #   apply(c(1,2),function(p) mean(p<=level)) |>
  #   as_tibble() |>
  #   mutate(kernel=kernelnames, pstr=kernelpstr, window=windowlabel, .before=1) |>
  #   pivot_longer(cols=all_of(F_names), names_to="Parameters", values_to="rejectrate")

res <- doForeach(varList, seed="seq", sfile=NULL, cluster = NULL,
                 doOne=doBerkowitz, block.size = blk_size, cores = num_cores,
                 exports = ls(), extraPkgs = (.packages()))
val <- getArray(res) |>
     apply(c(1,2),function(p) mean(p<=level)) |>
     as_tibble() |>
     rename('Rejection Rate'=1) |>
     mutate(F=F_names, .before=1) 

gttabl <- gt(val) |>
  fmt_percent(columns=where(is.numeric), decimals=1) |>
  tab_header(
    title='Size and power of Berkowitz LR-tests'
  ) |>
  tab_footnote(paste0('2^',log(nsims,2), ' trials with ', n_days,
                      " observations per trial.  Kernel window is ",
                      paste0(kernwindow,collapse=", "),"].") )
gtfile<- paste0(table_location,gtsavename) 
gt::gtsave(gttabl,filename = paste0(gtfile,".html"), inline_css=TRUE)
#gt::gtsave(gttabl,filename = paste0(gtfile,".tex"))
