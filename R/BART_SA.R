n11 = 1200
n10 = 250
n01 = 200
n00 = 100

pu = .5

py = 0.3
pz = 0.2

alpha = matrix(rep(.10, 4), nrow = 4)

prob.calc <- function(alpha, n11, n10, n01, n00, pu, rho_y, rho_z) {

	n = n11+n10+n01+n00
	nz1 = n11 + n01
	nz0 = n10 + n00
	ny1 = n11 + n10
	ny0 = n01 + n00

	A = matrix(c(n11, n10, n01, n00,
			n11*(n10+n00), -n10*(n11+n01), n01*(n10+n00), -n00*(n11+n01),
			n11*(n01+n00), n10*(n01+n00), -n01*(n11+n10), -n00*(n11+n10)),
		nrow = 3, byrow = T)

	b = matrix(c(n*pu, pz*n*sqrt(nz1*nz0*pu*(1-pu)), py*n*sqrt(ny1*ny0*pu*(1-pu))),
		nrow = 3)

	return(sum(abs(A%*%alpha-b)))
}


a.vals = optim(alpha, prob.calc, gr = NULL, method = "L-BFGS-B", 
		lower = 0, upper = 1, control = list(maxit = 1000, factr = 1e5),
		n11, n10, n01, n00, pu, py, pz)

alpha = a.vals$par

n = n11+n10+n01+n00
Y = c(rep(0, n01+n00), rep(1, n11+n10))
Z = c(rep(0, n00), rep(1, n01), rep(0, n10), rep(1, n11))
U = runif(n) < c(rep(alpha[4], n00), rep(alpha[3], n01), rep(alpha[2], n10), rep(alpha[1], n11))
cor(cbind(U, Y, Z))