
source("./helperFiles/simSetup.R")

i <- 750

doOne <- function(n,F,kernel){
  data <- choose_dist(F,n)
  PITs <- pnorm(data)
  purrr::map_dbl(as.list(kernel), ~spectral_Ztest(eval(parse(text=.x)), PITs))
}


# Table 2: unconditional-continuous_750.tex

#kern_vec2 <- c("BIN","ZU3","ZU","ZA","ZE","ZLp","ZLn","PE2","ZLL","ZPP","PE3","ZPUP")

kern_vec2 <- c("ZU","ZA","ZE","ZLL","ZPP","ZPUP","PNS","LLS","GS","GcS")


varList <-
  varlist( # constructor for an object of class 'varlist'
    n.sim = list(type="N", expr = quote(m), value = nsims),
    n = list(type="frozen", value = c(i)),
    F=list(type="grid", expr = quote(F),
           value = c("Normal", "Scaled t5", "Scaled t3")),
    kernel = list(type="inner",  expr = quote(Kernel),
                  value=kern_vec2))

# Narrow window analysis
alpha <- alpha_narrow
kernel_list <- define_kernels(alpha[1], alpha[2], alpha_star)
for (s in names(kernel_list)) assign(s, kernel_list[[s]])

savefilename <- paste(sim_location, "continuous_unconditional_narrow_", i, ".rds",sep="")
res <- doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
                 doOne=doOne, block.size = blk_size, cores = num_cores,
                 exports = ls(), extraPkgs = (.packages()))
narrow <- getArray(res)

# Wide window analysis
alpha <- alpha_wide
kernel_list <- define_kernels(alpha[1], alpha[2], alpha_star)
for (s in names(kernel_list)) assign(s, kernel_list[[s]])

savefilename <- paste(sim_location,"continuous_unconditional_wide_", i, ".rds",sep="")
res <- doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
                 doOne=doOne, block.size = blk_size, cores = num_cores,
                 exports = ls(), extraPkgs = (.packages()))
wide <- getArray(res)


val <- abind(narrow=narrow, wide=wide, along=0, use.dnns=TRUE)
names(dimnames(val)) <- c("window",names(dimnames(narrow)))
dimnames(val)

rv <- c("window","F") # row variables
cv <- c("kernel") # column variables
reject <- reject.rate(val)
dimnames(reject)$kernel <- kern_vec2
dimnames(reject)$kernel[dimnames(reject)$kernel=="ZLp"] <- "\\ZLp"
dimnames(reject)$kernel[dimnames(reject)$kernel=="ZLn"] <- "\\ZLn"
reject.f <- formatC(reject, digits=1, format="f")
fres <- ftable(reject.f, row.vars = rv, col.vars = cv)

tabL <- toLatex(fres,
                fontsize="normalsize",
                caption=paste("Estimated size and power of unconditional continuous tests. Replications = ",
                              nsims,
                              ". The number of days in the backtest sample is n = ",
                              i,
                              ". The narrow window is ",narrow_name," and the wide window is ",
                              wide_name,".",sep=""),
                vlist = varList,label=paste0("table:unconditional-continuous_", i))

sink(file=paste(table_location,"unconditional-continuous_", i, ".tex",sep=""))
print(tabL)
sink()
