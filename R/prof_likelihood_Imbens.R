p = seq(.1, .9, by = .1)
linecols = palette(rainbow(9))

par(mfrow = c(2,2))

alpha = seq(-10,10, length.out = 1000)
prof.lik = matrix(NA, nrow = 1000, ncol = length(p))
for(i in 1:length(p)) {
for(j in 1:1000){
	prof.lik[j,i] = llik(startvals, Y, W, X, alpha[j], 0.5,p[i])
	
}}

plot(alpha, prof.lik[,1], type = "l", col = linecols[1],
		xlab = "alpha", ylab = "Profile log-likelihood", ylim = c(-2400, -1700))
for(i in 2:length(p))
	lines(alpha, prof.lik[,i], col = linecols[i])



delta = seq(-10,10, length.out = 1000)
prof.lik = matrix(NA, nrow = 1000, ncol = length(p))
for(i in 1:length(p)) {
for(j in 1:1000){
	prof.lik[j,i] = llik(startvals, Y, W, X, 0.5, delta[j],p[i])
	
}}

plot(delta, prof.lik[,1], type = "l", col = linecols[1],
		xlab = "delta", ylab = "Profile log-likelihood", ylim = c(-2200, -1700))
for(i in 2:length(p))
	lines(delta, prof.lik[,i], col = linecols[i])



tau = seq(0,5, length.out = 1000)
prof.lik = matrix(NA, nrow = 1000, ncol = length(p))
for(i in 1:length(p)) {
for(j in 1:1000){
	prof.lik[j,i] = llik(c(startvals[-length(startvals)], tau[j]), Y, W, X, 0.5, 0.5,p[i])
	
}}

plot(tau, prof.lik[,1], type = "l", col = linecols[1],
		xlab = "tau", ylab = "Profile log-likelihood", ylim = c(-2000, -1700))
for(i in 2:length(p))
	lines(tau, prof.lik[,i], col = linecols[i])




s2 = seq(20,60, length.out = 1000)
prof.lik = matrix(NA, nrow = 1000, ncol = length(p))
for(i in 1:length(p)) {
for(j in 1:1000){
	startvals[length(startvals)-1] = s2[j]
	prof.lik[j,i] = llik(startvals, Y, W, X, 0.5, 0.5,p[i])
	
}}

plot(s2, prof.lik[,1], type = "l", col = linecols[1],
		xlab = "sigma^2", ylab = "Profile log-likelihood", ylim = c(-2000, -1700))
for(i in 2:length(p))
	lines(s2, prof.lik[,i], col = linecols[i])