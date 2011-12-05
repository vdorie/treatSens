##############
#Parameter values - will be filled in by user
##############

n11 = 200
n10 = 50
n01 = 100
n00 = 150

pu = .75

py = 0.3
pz = 0.2

############
#Calculate the maximum possible values of correlations
#n_ij: count of observations with Y = i and Z = j
#pU: marginal Pr(U = 1)
############

max.phi <- function(n11, n10, n01, n00, pU){
	n = n11+n10+n01+n00
	nz1 = n11 + n01
	nz0 = n10 + n00
	ny1 = n11 + n10
	ny0 = n01 + n00

	py = ny1/n
	pz = nz1/n

	a_z = min(pz, pU)
	a_y = min(py, pU)
	b_z = max(pz, pU)
	b_y = max(py, pU)

	return(list(z_max = sqrt(a_z*(1-b_z)/(b_z*(1-a_z))), y_max = sqrt(a_y*(1-b_y)/(b_y*(1-a_y)))))
}

###############
#Calculate a_yz = pr(U = 1| Y = y, Z = z) by fixing a_10
#n_ij: count of observations with Y = i and Z = j
#pU: marginal Pr(U = 1)
#rho_y, rho_z: desired correlations between U and Y or Z
#p10: fixed Pr(U = 1|Y = 1, Z = 0)
###############

prob.calc <- function(n11, n10, n01, n00, pU, rho_y, rho_z, p10) {

	#Check that specified correlations are valid
	aaa = max.phi(n11, n10, n01, n00, pU)
	if(abs(pz) > aaa$z_max)
		stop("Z-U correlation outside maximum range: (", -aaa$z_max, ",", aaa$z_max,")\n")
	if(abs(py) > aaa$y_max)
		stop("Y-U correlation outside maximum range: (", -aaa$y_max, ",", aaa$y_max,")\n")

	#calculate marginal totals
	n = n11+n10+n01+n00
	nz1 = n11 + n01
	nz0 = n10 + n00
	ny1 = n11 + n10
	ny0 = n01 + n00

	#system of linear equations to solve
	A = matrix(c(0, 1, 0,0,
			n11, n10, n01, n00,
			n11*(n10+n00), -n10*(n11+n01), n01*(n10+n00), -n00*(n11+n01),
			n11*(n01+n00), n10*(n01+n00), -n01*(n11+n10), -n00*(n11+n10)),
		nrow = 4, byrow = T)

	b = matrix(c(p10, n*pU, rho_z*n*sqrt(nz1*nz0*pU*(1-pU)), rho_y*n*sqrt(ny1*ny0*pU*(1-pU))),
		nrow = 4)
	
	#return vector (a11, a10, a01, a00)
	return(solve(A,b))
}


###############
#Generate U 
#n_ij: count of observations with Y = i and Z = j
#pU: marginal Pr(U = 1)
#rho_y, rho_z: desired correlations between U and Y or Z
###############

binaryYZU <- function(n11, n10, n01, n00, pU, rho_y, rho_z)
	ct = 0
	while(sum(alpha > 1)+sum(alpha<0) & ct < 1e4) {
	alpha = prob.calc(n11, n10, n01, n00, pU, rho_y, rho_z, runif(1))
	ct = ct +1
	}

	if(ct == 1e4){
		stop("No viable solution after 1e4 iterations")
	}else{
		U = runif(n) < c(rep(alpha[4], n00), rep(alpha[3], n01), rep(alpha[2], n10), rep(alpha[1], n11))
	}
	return(U)
}

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