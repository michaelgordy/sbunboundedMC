#setwd(here::here())

source("./helperFiles/simSetup.R")

#setwd(here::here())

i <- 750
# Part 1. Size

doOne <- function(n,CVT,kernel){
  PITs <- runif(n)
  kernfunc <- eval(parse(text=kernel))
  if (CVT == "None")
    return(spectral_Ztest(kernfunc,PITs))
  else{
    CVTfunc <- switch(kernfunc$type,
                      mono = eval(parse(text=paste(CVT,"_mono",sep=""))),
                      bi = eval(parse(text=paste(CVT,"_bi",sep=""))))
    return(spectral_MDtest(kernfunc,CVTfunc,PITs))
  }
}

varList <-
  varlist( # constructor for an object of class 'varlist'
    n.sim = list(type="N", expr = quote(m), value = nsims),
    n = list(type="frozen", value = c(i)),
    CVT = list(type="grid", expr = quote(ITT), value=c("None","EM","V.BIN","V.4","V.half")),
    kernel = list(type="frozen",  expr = quote(Test), value=c("ZU")))

# Narrow window analysis
alpha <- alpha_narrow

kernel_list <- define_kernels(alpha[1], alpha[2], alpha_star)
for (s in names(kernel_list)) assign(s, kernel_list[[s]])

cvt_list <- define_cvts(4)
for (s in names(cvt_list)) assign(s, cvt_list[[s]])

savefilename <- paste(sim_location,"conditional_narrow_size_EXTRACT.rds",sep="")
res <- doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
                 doOne=doOne, block.size = blk_size, cores = num_cores,
                 exports = ls(), extraPkgs = (.packages()))
narrow <- getArray(res)

val <- abind(narrow=narrow, along=0, use.dnns=TRUE)
names(dimnames(val)) <- c("window",names(dimnames(narrow)))

cv <- c("CVT") # column variables
reject <- reject.rate(val)
colnames_reject <- colnames(reject) # Extracting column names
values_reject <- as.numeric(reject)
values_reject_app <- append("None", values_reject)
colnames_reject_app <- append("Serial Dependence | CVT", colnames_reject)


# Creating the array, naming the second element in list to be CVT
reject_app <- array(c(values_reject_app), dim = c(1, 6), dimnames = list("F" = "Normal", "CVT" = colnames_reject_app))
format(as.numeric(reject_app[1, 2:6]), digits = 2)
reject.f <- formatC(reject_app, digits=1, format="f")
format(reject_app, digits = 2)
fres_1 <- ftable(reject.f, col.vars = cv)
fres_1[1, 2:6] <- format(as.numeric(fres_1[1, 2:6]), digits = 2)
fres_1

fres_1_df <- as.data.frame(fres_1)

fres_1_df_values <- fres_1_df %>%
  tidyr::spread(CVT, Freq)

# Part 2. Power

varList <-
  varlist( # constructor for an object of class 'varlist'
    n.sim = list(type="N", expr = quote(m), value = nsims),
    n = list(type="frozen", value = c(n_days)),
    cormodel = list(type="frozen",expr = quote(cormodel), value = list(ar=0.95,ma=-0.85)),
    F=list(type="grid", expr = quote(F),
           value = c("Normal", "Scaled t5")),
    CVT = list(type="grid", expr = quote(ITT), value=c("None","EM","V.BIN","V.4","V.half")),
    kernel = list(type="frozen",  expr = quote(Test), value=c("ZU")))

doOne <- function(n,cormodel,F,CVT,kernel){
  U <- rHS_arma(n,cormodel)
  X <- change_dist(F,U)
  PITs <- pnorm(X)
  kernfunc <- eval(parse(text=kernel))
  if (CVT == "None")
    return(spectral_Ztest(kernfunc,PITs))
  else{
    CVTfunc <- switch(kernfunc$type,
                      mono = eval(parse(text=paste(CVT,"_mono",sep=""))),
                      bi = eval(parse(text=paste(CVT,"_bi",sep=""))))
    return(spectral_MDtest(kernfunc,CVTfunc,PITs))
  }
}

# Narrow window analysis
alpha <- alpha_narrow

kernel_list <- define_kernels(alpha[1], alpha[2], alpha_star)
for (s in names(kernel_list)) assign(s, kernel_list[[s]])

cvt_list <- define_cvts(4)
for (s in names(cvt_list)) assign(s, cvt_list[[s]])

savefilename <- paste(sim_location,"conditional_narrow_power_EXTRACT.rds",sep="")
res <- doForeach(varList, seed="seq", sfile=savefilename, cluster = NULL,
                 doOne=doOne, block.size = blk_size, cores = num_cores,
                 exports = ls(), extraPkgs = (.packages()))

narrow <- getArray(res)

val <- abind(narrow=narrow, along=0, use.dnns=TRUE)
names(dimnames(val)) <- c("window",names(dimnames(narrow)))

cv <- c("CVT")
reject <- reject.rate(val)
reject.f <- formatC(reject, digits=1, format="f")


fres_2 <- ftable(reject.f, col.vars = cv)

### Create another column variable ("Serial Dependence")

augmented_fres2 <- as.numeric(c(1, fres_2[1 ,], 1, fres_2[2, ]))

fres_2_app <- array(matrix(as.numeric(c(1, fres_2[1 ,], 1, fres_2[2, ])), nrow = 2, byrow = TRUE),
                    dim = c(2, 6),
                    dimnames = list("F" = c("Normal", "Scaled t5"), "CVT" = colnames_reject_app))


fres2_df <- as.data.frame(fres_2)
fres2_df %<>% dplyr::mutate("Serial Dependence | CVT" = "ARMA(1, 1)")
fres2_df %<>% dplyr::select(-window)

fres2_df_values <- fres2_df %>%
  tidyr::spread(CVT, Freq) %>%
  dplyr::mutate_all(as.character)


combined_table <- rbind(fres_1_df_values, fres2_df_values)
print(xtable(combined_table,
             caption=paste("Extracted Conditional Table. Kernel = ZU. Window = narrow. Replications = ",
                                          nsims,
                                          ". The number of days in the backtest sample is n = ",
                                          n_days,
                                          ". The narrow window is ", narrow_name,".", sep=""),
             label="table:extracted_table"),
      include.rownames = FALSE,
      booktabs = T,
      file = paste0(getwd(), "/tables/extracted_test.tex"))
