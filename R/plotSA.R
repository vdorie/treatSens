#############
#plotSA - plot results of sensitivity analysis
#############

plotSA = function(x,...) {
	Zcors = apply(x@trt.cor, 2, mean)
	Ycors = apply(x@resp.cor, 1, mean)
	clevels = c(0,round(x@tau0, 3))
	contourCoords = contourLines(Zcors, Ycors, apply(x@tau, c(1,2), mean), levels = clevels,...)
	
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

	#invisible(lapply(contourCoords, function(x) lines(x[[3]] ~ x[[2]], col = "red")))

	contour(Zcors, Ycors, apply(x@tau, c(1,2), mean), levels = clevels, 
		lty = 0,add = T,method = "edge",...)
}

plot(test.run)