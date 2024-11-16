source("R/simSetup.R")
n_days <- 500
nsims <- 2^16
blk_size <- nsims/2^5
level <- 0.05
savedata <- FALSE  # TRUE to have simsalapar save simulation data
gtsavename <- 'berkowitz'

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

# doOne <- function(n,Fmodel,kernel){
#   PIT <- Fmodel_list[[Fmodel]](n) |> pnorm() |> 
#     pmin(1-.Machine$double.eps)
#   purrr::map_dbl(kernel, ~spectral_Ztest(.x, PIT))
# }

doBerkowitz <- function(n,Fmodel, alpha1){
  PIT <- Fmodel_list[[Fmodel]](n) |> pnorm() |> 
    pmin(1-.Machine$double.eps)
  Berkowitz_test(PIT, alpha1)
}

F_names <- c("Normal", "Scaled t10", "Scaled t5", "Scaled t3")
alf1_vec <- c(alpha_narrow[1],alpha_wide[1])

varList <- varlist( # constructor for an object of class 'varlist'
    n.sim = list(type="N", expr = quote(m), value = nsims),
    n = list(type="frozen", value = c(n_days)),
    Fmodel=list(type="grid", expr = quote(F), value = F_names),
    alpha1=list(type="grid", expr = quote(alf1), value = alf1_vec))
  
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
     as_tibble()



# narrow <- purrr::map(betakerns[kern_vec2], 
#                      ~makeBetaKernel(.,alpha_narrow)) |>
#   rejectionrate(alpha_narrow,"narrow")
# wide <- purrr::map(betakerns[kern_vec2], 
#                    ~makeBetaKernel(.,alpha_wide)) |>
#   rejectionrate(alpha_wide,"wide")
# res <- dplyr::bind_rows(narrow,wide) |>
#   mutate(notBeta=if_else(stringr::str_starts(kernel,"Beta")," ",kernel),
#          collab=paste(notBeta,pstr, sep="<br>")) |>
#   dplyr::select(-kernel,-pstr,-notBeta) |>
#   pivot_wider(names_from = collab, values_from = rejectrate) 
# 
# save(narrow,wide,res,alpha_wide,alpha_narrow,nsims,n_days,
#      file=paste0(sim_location,gtsavename,".RData"))
# 
# 
# # Narrow window analysis
# alpha <- alpha_narrow
# kernel_list <- define_kernels(alpha[1], alpha[2], alpha_star)
# for (s in names(kernel_list)) assign(s, kernel_list[[s]])
# 
# lr_test_list <- define_lr_tests(alpha[1], alpha[2], alpha_star)
# for (s in names(lr_test_list)) assign(s, lr_test_list[[s]])
# 
# savefilename <- paste(sim_location,"Z-LRcomparison_narrow_", i, ".rds",sep="")
# res <- doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
#                  doOne=doOne, block.size = blk_size, cores = num_cores,
#                  exports = ls(), extraPkgs = (.packages()))
# narrow <- getArray(res)
# 
# # Wide window analysis
# alpha <- alpha_wide
# kernel_list <- define_kernels(alpha[1], alpha[2], alpha_star)
# for (s in names(kernel_list)) assign(s, kernel_list[[s]])
# 
# lr_test_list <- define_lr_tests(alpha[1], alpha[2], alpha_star)
# for (s in names(lr_test_list)) assign(s, lr_test_list[[s]])
# 
# res <- doForeach(varList, seed="seq", cluster = NULL,
#                  doOne=doOne, block.size = blk_size, cores = num_cores,
#                  exports = ls(), extraPkgs = (.packages()))
# wide <- getArray(res)
# 
# 
# val <- abind(narrow=narrow, wide=wide, along=0, use.dnns=TRUE)
# names(dimnames(val)) <- c("window",names(dimnames(narrow)))
# dimnames(val)
# 
# rv <- c("window","F") # row variables
# cv <- c("test") # column variables
# reject <- reject.rate(val)
# dimnames(reject)$test[dimnames(reject)$test=="ZLp"] <- "\\ZLp"
# dimnames(reject)$test[dimnames(reject)$test=="ZLn"] <- "\\ZLn"
# reject.f <- formatC(reject, digits=1, format="f")
# fres <- ftable(reject.f, row.vars = rv, col.vars = cv)
# 
# tabL <- toLatex(fres,
#                 fontsize="normalsize",
#                 caption=paste("Estimated size and power of unconditional Z-tests and LR-tests. Replications = ",
#                               nsims,
#                               ". The number of days in the backtest sample is n = ",
#                               i,
#                               ". The narrow window is ",narrow_name," and the wide window is ",
#                               wide_name,".",sep=""),
#                 vlist = varList,label=paste0("table:unconditional-comparison_", i))
# 
# sink(file=paste(table_location,"unconditional-comparison_", i, ".tex",sep=""))
# print(tabL)
# sink()