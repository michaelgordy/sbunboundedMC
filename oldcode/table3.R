#setwd(here::here())

source("./helperFiles/simSetup.R")

#setwd(here::here())


i <- 750

doOne <- function(n,F,test){
  data <- choose_dist(F,n)
  PITs <- pnorm(data)
  purrr::map_dbl(test,
                 function(krn) {
                   if (substr(krn,1,1) != "L")
                     spectral_Ztest(eval(parse(text=krn)),PITs)
                   else
                     spectral_LRtest(eval(parse(text=krn)),PITs)
                 })
}

varList <-
  varlist( # constructor for an object of class 'varlist'
    n.sim = list(type="N", expr = quote(m), value = nsims),
    n = list(type="frozen", value = c(i)),
    F=list(type="grid", expr = quote(F),
           value = c("Normal", "Scaled t5", "Scaled t3")),
    test = list(type="grid",  expr = quote(Kernel),
                value=c("BIN","LR1","PE2","LR2","PE3","LR3","ZPP","LRB")))

# Narrow window analysis
alpha <- alpha_narrow
kernel_list <- define_kernels(alpha[1], alpha[2], alpha_star)
for (s in names(kernel_list)) assign(s, kernel_list[[s]])

lr_test_list <- define_lr_tests(alpha[1], alpha[2], alpha_star)
for (s in names(lr_test_list)) assign(s, lr_test_list[[s]])

savefilename <- paste(sim_location,"Z-LRcomparison_narrow_", i, ".rds",sep="")
res <- doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
                 doOne=doOne, block.size = blk_size, cores = num_cores,
                 exports = ls(), extraPkgs = (.packages()))
narrow <- getArray(res)

# Wide window analysis
alpha <- alpha_wide
kernel_list <- define_kernels(alpha[1], alpha[2], alpha_star)
for (s in names(kernel_list)) assign(s, kernel_list[[s]])

lr_test_list <- define_lr_tests(alpha[1], alpha[2], alpha_star)
for (s in names(lr_test_list)) assign(s, lr_test_list[[s]])

res <- doForeach(varList, seed="seq", cluster = NULL,
                 doOne=doOne, block.size = blk_size, cores = num_cores,
                 exports = ls(), extraPkgs = (.packages()))
wide <- getArray(res)


val <- abind(narrow=narrow, wide=wide, along=0, use.dnns=TRUE)
names(dimnames(val)) <- c("window",names(dimnames(narrow)))
dimnames(val)

rv <- c("window","F") # row variables
cv <- c("test") # column variables
reject <- reject.rate(val)
dimnames(reject)$test[dimnames(reject)$test=="ZLp"] <- "\\ZLp"
dimnames(reject)$test[dimnames(reject)$test=="ZLn"] <- "\\ZLn"
reject.f <- formatC(reject, digits=1, format="f")
fres <- ftable(reject.f, row.vars = rv, col.vars = cv)

tabL <- toLatex(fres,
                fontsize="normalsize",
                caption=paste("Estimated size and power of unconditional Z-tests and LR-tests. Replications = ",
                              nsims,
                              ". The number of days in the backtest sample is n = ",
                              i,
                              ". The narrow window is ",narrow_name," and the wide window is ",
                              wide_name,".",sep=""),
                vlist = varList,label=paste0("table:unconditional-comparison_", i))

sink(file=paste(table_location,"unconditional-comparison_", i, ".tex",sep=""))
print(tabL)
sink()