##############
#Parameter values - will be filled in by user
##############

n11 = 200
n10 = 50
n01 = 100
n00 = 150

pu = .5

py = 0.3
pz = 0.2

#############
#Demonstration code
#############

n = n11+n10+n01+n00
Y = c(rep(0, n01+n00), rep(1, n11+n10))
Z = c(rep(0, n00), rep(1, n01), rep(0, n10), rep(1, n11))

corr.sim = matrix(NA, nrow = 1000, ncol = 2)
for(i in 1:1000){
		U <- binaryYZU(n11, n10, n01, n00, pu, py, pz)
		corr.sim[i,] = cor(cbind(U, Y, Z))[1,2:3]
	}
	alpha[1] = 1.2
}

apply(corr.sim,2,mean)
apply(corr.sim,2,sd)

###could also generate U with constant number of 1's that is closest to prop. given by alpha:

corr.sim2 = matrix(NA, nrow = 1000, ncol = 2)

for(i in 1:1000){
	ct = 0
	while(sum(alpha > 1)+sum(alpha<0) & ct < 1e4) {
	alpha = prob.calc(n11, n10, n01, n00, pu, py, pz, runif(1))
	ct = ct +1
	}

	if(ct == 1e4){
		print("Error: no viable solution after 1e4 iterations")
	}else{
		u11 = rep(0, n11)
		u10 = rep(0, n10)
		u01 = rep(0, n01)
		u00 = rep(0, n00)

		u1.cts = round(alpha*c(n11,n10,n01,n00))
		
		u11[sample(1:n11, u1.cts[1])] =1
		u10[sample(1:n10, u1.cts[2])] =1
		u01[sample(1:n01, u1.cts[3])] =1
		u00[sample(1:n00, u1.cts[4])] =1
	
		U = c(u00, u01, u10, u11)
		corr.sim2[i,] = cor(cbind(U, Y, Z))[1,2:3]
	}
	alpha[1] = 1.2
}

apply(corr.sim2,2,mean)
apply(corr.sim2,2,sd)


##########
#Plots
##########
setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
pdf("BART_SA_boxplots.pdf")

boxplot(corr.sim, names = c(paste("rho_y =", py), paste("rho_z =", pz)),
	main = "U random by probability")

boxplot(corr.sim2, names = c(paste("rho_y =", py), paste("rho_z =", pz)),
	main = "U random by permutation")


dev.off()