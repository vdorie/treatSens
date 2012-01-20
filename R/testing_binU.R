setwd("C:/Users/Nicole/Documents/causalSA/R_package/trunk/R")
source("GLM_sens.R")

n = 500
rho = seq(-0.3, 0.3, by = 0.05)
Z = rnorm(n)

y.means <- z.means <- y.nas <- z.nas <- matrix(NA, nrow = length(rho), ncol = 20, dimnames = list(round(rho,2),round(seq(-0.5, 0.5, length.out = 20),2)))
ym <- zm <- yna <- zna <- array(NA, dim = c(20,20,length(rho)), dimnames = list(round(seq(-0.5, 0.5, length.out = 20),2),round(seq(-0.5, 0.5, length.out = 20),2),round(rho,2)))
rhoyz = vector()

for(i in 1:length(rho)) {
	Y = Z + rnorm(n, mean = 0, sd = sqrt((1-rho[i]^2)/rho[i]^2))
	rhoyz[i] = cor(Y,Z)

	test.run <- GLM.sens(Y~Z, resp.rho.vals = 20, trt.rho.vals = 20, standardize = F, resp.rho.range = c(-0.5,0.5), trt.rho.range = c(-0.5,0.5),
			trt.family = gaussian,
			resp.family = gaussian,
			conf.model = "binomial",
			nsim = 200)
	
	ym[,,i] <- apply(test.run[[5]], c(1,2), mean, na.rm = T)
	y.means[i,] <- apply(ym[,,i], 1, mean, na.rm = T)
	zm[,,i] <- apply(test.run[[6]], c(1,2), mean, na.rm = T)
	z.means[i,] <- apply(zm[,,i], 2, mean, na.rm = T)
	
	yna[,,] <- apply(is.na(test.run[[5]]), c(1,2), sum)
	y.nas[i,] <- apply(yna[,,i], 1, sum)

	zna[,,i] <- apply(is.na(test.run[[6]]), c(1,2), sum)
	z.nas[i,] <- apply(zna[,,i], 2, sum)
}

save(ym, zm, y.means, z.means, yna, zna, y.nas, z.nas, file = "testing_binU_adj.RData")

#load("testing_binU.RData")

yzvals = seq(-0.5, .5, length.out = 20)
image(yzvals, rho, ym[7,,])
image(yzvals, rho, zm[,2,])

pdf("testing_binU_adj.pdf")
plot(yzvals, apply(y.means,2,mean), xlab = "Target correlations", ylab = "Observed correlations", ylim = c(-0.5, 0.5))
fit = lm(apply(y.means,2,mean)~yzvals)
abline(fit$coef)
summary(fit)

points(yzvals, apply(z.means,2,mean), col = "red")
fit = lm(apply(z.means,2,mean)~yzvals)
abline(fit$coef, col = "red")
summary(fit)

legend("topleft", legend = c("Y", "Z"), col = c("black", "red"))
dev.off()