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
