
n11 = c(50, 100, 150, 150, 200)
n10 = c(200, 150, 100, 100, 50)
n01 = c(150, 150, 150, 100, 100)
n00 = c(100, 100, 100, 150, 150)

pu = .5

py = c(-.5, -0.3, -0.1, 0.1, 0.3, .5)
pz = c(-.5, -0.3, -0.1, 0.1, 0.3, .5)

meanY <- sdY <-meanZ <- sdZ <- NAY <- NAZ <- array(NA, dim = c(5, 6, 6), 
	dimnames = list(paste("OR",round(n11*n00/(n01*n10),2)), paste("rhoY",py), paste("rhoZ",pz)))

for(i in 1:5) { for(j in 1:6) { for(k in 1:6) {

n = n11[i]+n10[i]+n01[i]+n00[i]
Y = c(rep(0, n01[i]+n00[i]), rep(1, n11[i]+n10[i]))
Z = c(rep(0, n00[i]), rep(1, n01[i]), rep(0, n10[i]), rep(1, n11[i]))

corr.sim = matrix(NA, nrow = 100, ncol = 2)
for(s in 1:100){
		U <- try(binaryYZU(n11[i], n10[i], n01[i], n00[i], pu, py[j], pz[k]))
		if(is.logical(U)) corr.sim[s,] = cor(cbind(U, Y, Z))[1,2:3]
	}

meanY[i,j,k] = apply(corr.sim,2,mean, na.rm = T)[1]
meanZ[i,j,k] = apply(corr.sim,2,mean, na.rm = T)[2]
sdY[i,j,k] = apply(corr.sim,2,sd, na.rm = T)[1]
sdZ[i,j,k] = apply(corr.sim,2,sd, na.rm = T)[2]
NAY[i,j,k] = sum(is.na(corr.sim[,1]))
NAZ[i,j,k] = sum(is.na(corr.sim[,2]))

cat("completing i =", i, ", j =", j, ", k =",k,"\n")
}}}
