setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("BART_cont.R")

#load LaLonde data (or replace with your own test case!)
library(foreign)
lalonde<-read.dta("C:/Users/Nicole/Documents/causalSA/R_package/trunk/data/lalonde data.dta")

re74per <- with(lalonde,re74/1000)
re75per <- with(lalonde,re75/1000)
X <- with(lalonde,cbind(married,age,black,hisp,educ,re74per,reo74,re75per,reo75))
Z <- with(lalonde,t)
Y <- with(lalonde,re78/1000)

#Generate a continuous treatment from Z
#Z <- (Z+1)*log(Y+1)+rnorm(length(Z),0,5)
Z <-  (Y+apply(2*X,1,sum)+rexp(length(Z), 1/10))^2

test.run <- BART.sens.cont(Y[!is.na(Y)]~Z[!is.na(Y)]+X[!is.na(Y),], grid.dim = c(3,3), standardize = T,
		U.model = "normal",
		verbose = T,
		nsim = 10)

save(test.run, file = "test.run.BART.RData")
#Check out processing functions:
summary(test.run)
print(test.run)
plot(test.run)
slotNames(test.run)
test.run

