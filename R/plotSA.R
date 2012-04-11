#############
#plotSA - plot results of sensitivity analysis
#############

plotSA = function(x, coef.axes = F,
			contour.levels = NULL,
			zero.col = "red",
			lty.zero = 1,
			insig.col = "blue",
			lty.insig = 1,
			data.line = TRUE,
			X.pch = 19,
			signif.level = 0.05,
			labcex = 0.75,
			...) {
#note in help: if contours are too rough, up nsim in sens fn
	if(!coef.axes) {
		Zcors = apply(x@trt.cor, 2, mean, na.rm = T)
		Ycors = apply(x@resp.cor, 1, mean, na.rm = T)
		xlab = "Treatment partial correlation"
	 	ylab = "Response partial correlation"
		Xpart = x@Xpartials
	}else{
		Zcors = apply(x@alpha, 2, mean, na.rm = T)
		Ycors = apply(x@delta, 1, mean, na.rm = T)
		xlab = "Alpha"
		ylab = "Delta"
		Xpart = x@Xcoef
		#note that due to correlation among Xs, some may not appear on plot
		#because observed partial cors don't map directly to coefs in this case
		#forcing inclusion can lead to difficult to read plot.
	}

	taus = t(apply(x@tau, c(1,2), mean, na.rm = T))

	if(is.null(contour.levels)){
#need to robustify if selected cell is NA
		exTau = ifelse(sign(x@tau0)==1, taus[dim(taus)[1], dim(taus)[2]], taus[1,dim(taus)[2]])
		clevels = round(seq(x@tau0*0.8, exTau*.8, length.out = 8),2)
	}else{
		clevels = contour.levels
	} 
	
	par(mgp = c(2,.5,0))
	plot(Xpart, xlim = c(min(Zcors, na.rm = T),max(Zcors, na.rm = T)), ylim = c(min(Ycors, na.rm = T),max(Ycors, na.rm = T)),
		pch = X.pch, xlab = xlab, ylab = ylab,...)

	abline(h = 0)
	abline(v = 0)

	contour(Zcors, Ycors, taus, levels = clevels, 
		add = T, labcex = labcex,...)

	contour(Zcors, Ycors, taus, levels = 0, 
		add = T, col = zero.col,lty = lty.zero,labcex = labcex,...)

	contour(Zcors, Ycors, taus/apply(x@se.tau, c(1,2), mean), labels = "N.S.",
		levels = -sign(x@tau0)*qnorm(signif.level/2), add = T, col = insig.col,
		lty = lty.insig, labcex = labcex,...)
	
	legend(0.8*max(Zcors), 0, legend = round(x@tau0,2), cex = labcex,
		yjust = 0.5, x.intersp = 0,
		bg = ifelse(par("bg")== "transparent", "white", par("bg")), box.lty = 0)

	if(data.line) {
		proj.pts = apply(Xpart, 1, mean)
		max.pt = Xpart[proj.pts == max(proj.pts),]
		zcor = (1:length(Zcors))[abs(Zcors-max.pt[1]) ==  min(abs(Zcors-max.pt[1]))]
		if(zcor > max.pt[1]){
		 zpts = c(zcor-1, zcor)
		}else{
		 zpts = c(zcor, zcor+1)
		}
		ycor = (1:length(Ycors))[abs(Ycors-max.pt[2]) ==  min(abs(Ycors-max.pt[2]))]
		if(ycor > max.pt[2]){
		 ypts = c(ycor-1, ycor)
		}else{
		 ypts = c(ycor, ycor+1)
		}
		clevel = ((Zcors[zpts[2]] - Zcors[zpts[1]])*(Ycors[ypts[2]] - Ycors[ypts[1]]))^(-1)*
			sum(taus[zpts, ypts]*
			matrix(c(-(Zcors[zpts[2]] - max.pt[1])*(Ycors[ypts[1]] - max.pt[2]), 
				(Zcors[zpts[1]] - max.pt[1])*(Ycors[ypts[1]] - max.pt[2]),
				(Zcors[zpts[2]] - max.pt[1])*(Ycors[ypts[2]] - max.pt[2]),
				-(Zcors[zpts[1]] - max.pt[1])*(Ycors[ypts[2]] - max.pt[2])), 
				nrow = 2, byrow = T))
		contour(Zcors, Ycors, taus, levels = round(clevel,2),
			add = T, col = "grey",labcex = labcex,...)
	}
}

############
#plotSA.cont
############
plotSA.cont = function(x, coef.axes = F,
			signif.level = 0.05,
			labcex = 0.75,
			...) {
require(lattice)
#should do this with lattice to reduce whitespace
	par(mgp = c(2,.5,0))

	trt.ests <- sd.ests <- LCI <- UCI <- Z <- tau0 <- cell <- NULL
	ct = 0
	for(i in 1:dim(x@tau)[1]) {
	for(j in 1:dim(x@tau)[2]) {
		ct = ct+1
		trt.temp = apply(x@tau[i,j,,], 2, mean)
		sd.temp = apply(x@se.tau[i,j,,],2,mean)

		trt.ests = c(trt.ests, trt.temp)
		sd.ests = c(sd.ests, sd.temp)

		LCI = c(LCI,trt.temp + qnorm(signif.level/2)*sd.temp)
		UCI = c(UCI,trt.temp + qnorm(1-signif.level/2)*sd.temp)
		
		tau0 = c(tau0, x@tau0)
		Z = c(Z, x@Z)
		trt.cor = round(apply(x@trt.cor, 2, mean, na.rm = T),2)
		resp.cor = round(apply(x@resp.cor, 1, mean, na.rm = T),2)
		cell = c(cell, rep(paste("(", trt.cor[i], ",", resp.cor[j],")",sep = ""),length(x@Z)))
	}} 
		xyplot(trt.ests + LCI + UCI + tau0 ~ Z | cell, xlab = "Treatment", ylab = "Treatment effect", 
			layout = dim(x@tau)[c(1,2)], ylim = c(.95*min(c(LCI, x@tau0)), 1.05*max(c(UCI, x@tau0))),
			type = "smooth", col = c("black", "blue", "blue", "grey"),
			sub = "(partial r^2 treatment, partial r^2 response)", cex.sub = 0.5) 
}

#############
#Null plot for continuous treatment in BART
#############

plot.null = function(x, coef.axes = F,
			signif.level = 0.05,
			labcex = 0.75,
			...) {
	require(lattice)
	require(latticeExtra)
	par(mgp = c(2,.5,0))

	tau0 = x@tau0
	LCI = tau0 + qnorm(signif.level/2)*x@se.tau0
	UCI = tau0 + qnorm(1-signif.level/2)*x@se.tau0

	line.plots <- xyplot(tau0 + LCI + UCI ~ x@Z, xlab = "Treatment", ylab = "Treatment effect", 
		#ylim = c(.95*min(c(LCI, x@tau0)), 1.05*max(c(UCI, x@tau0))),
		type = "smooth", col = c("black", "blue", "blue"),
		scales = list(tck = c(1,0))
		) 
	hist.plot <- histogram(~ x@Z, freq = F, type = "density", col = "white",
			ylim = c(0, 4*max(hist(x@Z,plot=F)$density)))

	line.plots + as.layer(hist.plot, y.same = FALSE, axes = NULL, opposite = FALSE)
}




plotSA.smoothed = function(x,...) {
	Zcors = apply(x@trt.cor, 2, mean)
	Ycors = apply(x@resp.cor, 1, mean)

	taus = apply(x@tau, c(1,2), mean)

	clevels = round(seq(x@tau0*0.8, taus[dim(taus)[1], dim(taus)[2]]*.8, length.out = 8),2)
	contourCoords = contourLines(Zcors, Ycors, apply(x@tau, c(1,2), mean, na.rm =T), levels = clevels,...)

	zeroCoords = contourLines(Zcors, Ycors, taus, levels = 0)
	noSigCoords = contourLines(Zcors, Ycors, taus/apply(x@se.tau, c(1,2), mean), levels = sign(x@tau0)*2)
	
	fixVertical <- function(coords) {
		xRange = range(coords[[2]])
		if(xRange[2]-xRange[1] < 0.1){   #X range small means nearly vertical contour - smoothing won't work properly
			temp = coords[[2]]
			coords[[2]] = coords[[3]]
			coords[[3]] = temp
			return(list(coords = coords, switch = TRUE))
		}
		return(list(coords = coords, switch = FALSE))
	}

	lowessContour <- function(coords) {
		x = fixVertical(coords)
		lowCoords = lowess(x$coords[[3]] ~ x$coords[[2]], f = 0.5)
		if(x$switch){
			temp = lowCoords$x
			lowCoords$x = lowCoords$y
			lowCoords$y = temp
		}
		return(lowCoords)
	}

	plot(x@Xpartials, type = "p", xlim = c(min(Zcors),max(Zcors)), ylim = c(min(Ycors),max(Ycors)),
		pch = 19, xlab = "Treatment correlation", ylab = "Response correlation")

	invisible(lapply(contourCoords, function(x) lines(lowessContour(x))))

	invisible(lapply(noSigCoords, function(x) lines(lowessContour(x), col = "blue")))
	invisible(lapply(zeroCoords, function(x) lines(lowessContour(x), col = "red")))

	abline(h = 0)
	abline(v = 0)

	contour(Zcors, Ycors, apply(x@tau, c(1,2), mean), levels = clevels, 
		lty = 0, 
		add = T,...)

	contour(Zcors, Ycors, apply(x@tau, c(1,2), mean), levels = 0, 
		lty =0, 
		add = T, col = "red",...)

	contour(Zcors, Ycors, taus/apply(x@se.tau, c(1,2), mean), labels = "Insignificant",
		lty = 0, 
		levels = sign(x@tau0)*2, add = T, col = "blue")

}
