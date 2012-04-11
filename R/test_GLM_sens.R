setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("GLM_sens.R")
source("BART_sens.R")

#load LaLonde data (or replace with your own test case!)
library(foreign)
lalonde<-read.dta("C:/Users/Nicole/Documents/causalSA/R_package/trunk/data/lalonde data.dta")

re74per <- with(lalonde,re74/1000)
re75per <- with(lalonde,re75/1000)
X <- with(lalonde,cbind(married,age,black,hisp,educ,re74per,reo74,re75per,reo75))
Z <- with(lalonde,t)
Y <- with(lalonde,re78/1000)


#test grid analysis
#function arguments:	formula, 			#formula: assume treatment is 1st term on rhs
#				resp.family = gaussian,	#family for GLM of model for response
#				trt.family = gaussian,	#family for GLM of model for treatment
#				U.model = "binomial",	#form of model for confounder: can be one of "binomial" and "normal"
#				grid.dim = c(20,20),	#final dimensions of output grid
#				standardize = TRUE,	#Logical: should variables be standardized?
#				nsim = 20,			#number of simulated Us to average over per cell in grid
#				data = NULL			#data object containing variables, if applicable

ZN = 1-Z
test.norm <- GLM.sens(Y~Z+X, grid.dim = c(50,50), standardize = T,
		trt.family = binomial,
		resp.family = gaussian,
		U.model = "normal",
		verbose = T,
		nsim = 50)
test.neg  <- GLM.sens(Y~ZN+X, grid.dim = c(20,20), standardize = T,
		trt.family = binomial,
		resp.family = gaussian,
		U.model = "normal",
		verbose = T,
		nsim = 20)

#function arguments:	formula, 			#formula: assume treatment is 1st term on rhs
#				grid.dim = c(20,20),	#final dimensions of output grid				
#				est.type = "ATE",		#type of estimand desired: one of ATT, ATC, ATE.  Assumes binary Z with 1 = treatment
#				U.model = "binomial",	#form of model for confounder: can be one of "binomial" and "normal"
#				standardize = TRUE,	#Logical: should values be standardized?  If FALSE, force specification of ranges?
#				nsim = 20,			#number of simulated Us to average over per cell in grid
#				data = NULL

test.run <- BART.sens(Y~Z+X, grid.dim = c(20,20), standardize = T,
		est.type = "ATE",
		U.model = "binomial",
		nsim = 50)
save(test.bin, test.norm, file = "test.GLM.Ucomparison.RData")

#Check out processing functions:
summary(test.run)
print(test.run)
plot(test.run)
slotNames(test.run)
test.run

pdf("U model comparison.pdf", width = 12, height = 6)
par(mfrow = c(1,2))
plot(test.norm, main = "Continuous U")
debug(plotSA)
plot(test.bin, main = "Binary U")

dev.off()

##########
#test binary U
###########

tau  = seq(-3, 3, by = 0.5)

p_y = 0.3
p_z = 0.2

cYZ = cYU = cZU = vector()

for(j in 1:length(tau)) {
Z = rnorm(500, 10, 5)  #rbinom(500, 1, .5) 
Y = tau[j]*Z+ rnorm(500, 10, 20)

corr.sim = matrix(NA, nrow = 1000, ncol = 2)
for(i in 1:1000){
		U <- contYZbinaryU(Y, Z, p_y, p_z)
		corr.sim[i,] = cor(cbind(U, Y, Z))[1,2:3]
}

cYZ[j] = cor(Y,Z)
cYU[j] = apply(corr.sim,2,mean)[1]
cZU[j] = apply(corr.sim,2,mean)[2]
apply(corr.sim,2,sd)
}

##########
#Plots
##########
pdf("BART_SA_contYZbinaryU.pdf")

plot(cYZ, cYU, ylim = c(0.1, 0.4), type = "l")
lines(cYZ, cZU, col = "red")

abline(h = p_y)
abline(h = p_z, col = "red")

boxplot(corr.sim, names = c(paste("rho_y =", p_y), paste("rho_z =", p_z)),
	main = "Continuous U &Y")

dev.off()