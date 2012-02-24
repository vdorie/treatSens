library(foreign)
SWAY<-read.dta("C:/Users/Nicole/Documents/causalSA/R_package/trunk/data/SWAY_I_27Jan2009.dta")

wage_mo <- with(SWAY, wage_mo)
comm_mobil <-with(SWAY, comm_mobil)

Y<-wage_mo
Z <- with(SWAY, abd)
X <- with(SWAY, cbind(mthr_ed,#fthr_ed,no_fthr96,no_mthr96,orphan96,hh_size96,
		#hh_land,hh_cattle,hh_stock,hh_plow,A14,A15,A16,A17,A18,A19,A20,
		#A21,A22,A23,A24,A25,A26,A27,A28,A29,#A30,A31,A32,A33,A34,A35,
		C_ach,C_akw))#,C_ata,C_kma,C_oro,C_pad,C_paj))#,C_pal))

setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("BART_sens.R")
test.run.swayI <- BART.sens(Y[!is.na(Y)]~Z[!is.na(Y)] + X[!is.na(Y),],
					grid.dim = c(10,10),
					est.type = "ATE",
					U.model = "normal",
					standardize = TRUE,
					nsim = 10)					


test.run.swayIII <- GLM.sens(Y[!is.na(Y)]~Z[!is.na(Y)] + X[!is.na(Y),], resp.rho.vals = 100, trt.rho.vals = 100, standardize = T, resp.rho.range = c(-0.5,0.5), trt.rho.range = c(-0.5,0.5),
		trt.family = binomial,
		resp.family = binomial,
		U.model = "normal",
		nsim = 20)
rvals = seq(-.5, .5, length.out = 100)
tauS <- apply(test.run.swayII[[1]], c(1,2), mean)
contour(rvals, rvals, tauS)

