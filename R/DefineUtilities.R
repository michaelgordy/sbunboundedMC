choose_dist <- function(F, n){
  switch(F, 
         "Normal" = rnorm(n), 
         "t10" = rt(n,df=10), 
         "t5" = rt(n,df=5),
         "t3" = rt(n,df=3), 
         "Scaled t10" = sqrt(8)*rt(n, df = 10)/sqrt(10),
         "Scaled t5" = sqrt(3)*rt(n,df=5)/sqrt(5),
         "Scaled t3" = rt(n,df=3)/sqrt(3))
}

change_dist <- function(F, U){
  switch(F, 
         "Normal" = qnorm(U), 
         "t10" = qt(U,df=10), 
         "t5" = qt(U,df=5),
         "t3" = qt(U,df=3), 
         "Scaled t10" = sqrt(8)*qt(U, df = 10)/sqrt(10),
         "Scaled t5" = sqrt(3)*qt(U,df=5)/sqrt(5),
         "Scaled t3" = qt(U,df=3)/sqrt(3))
}


reject.rate <- function(val,percentage =TRUE,level=0.05){
  rejectfunc <- function(v){mean(as.numeric(v <= level))}
  non.sim.margins <- setdiff(names(dimnames(val)), "n.sim")
  reject <- apply(val, non.sim.margins, rejectfunc) 
  if (percentage)
    reject <- 100*reject
  reject
}

