setwd("C:/Users/Nicole/Documents/Graduate coursework/data analysis course/mtct")
mtct = read.table("childhiv.txt", header = T, na = ".")

setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("BART_cont.R")

Z = with(mtct, potato)
Y = with(mtct, childhiv)
X = with(mtct, cbind(caesar, everbreast, mnvp, bnvp, income, educ, age, bwt))

test.run <- BART.sens.cont(Y[!(is.na(Y)|is.na(Z))]~Z[!(is.na(Y)|is.na(Z))]+X[!(is.na(Y)|is.na(Z)),], grid.dim = c(3,3), standardize = T,
		U.model = "normal",
		verbose = T,
		nsim = 20)

save(test.run, file = "test.run.mtct.RData")