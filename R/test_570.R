hypo = dget("C:/Users/Nicole/Documents/Graduate coursework/Stat Methods 570s/570/midterm/hypogl.R")

names(hypo)
attach(hypo)
fit1 = glm(hyporate~age+gender+bmi+bg+thercode, family = poisson)
summary(fit1)
tc1 = thercode==1
tc2 = thercode==2

#note thercode is actual treatment, using age and bg here for purposes of testing continuous trt

setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("BART_cont.R")
load("test.run.570.RData")

test.age <- BART.sens.cont(hyporate~age+bmi+gender+bg+tc1+tc2, grid.dim = c(3,3), standardize = F,
		U.model = "normal",
		verbose = T,
		nsim = 20)

test.bg <- BART.sens.cont(hyporate~bg+age+gender+bmi+tc1+tc2, grid.dim = c(3,3), standardize = F,
		U.model = "normal",
		verbose = T,
		nsim = 20)

#test with actual treatment
source("BART_sens.R")
source("GLM_sens.R")
test.tc.glm <- GLM.sens(hyporate~thercode+bmi+age+gender+bg, grid.dim = c(20,20), standardize = F,
		trt.family = poisson,
		resp.family = poisson,
		U.model = "normal",
		verbose = T,
		nsim = 50)
test.tc <- BART.sens(hyporate~thercode+bmi+age+gender+bg, grid.dim = c(20,20), standardize = F,
		est.type = "ATE",
		U.model = "normal",
		nsim = 20)

save(test.age, test.bg, test.tc.glm, test.tc, 
		file = "test.run.570.RData")