#############
#plotSA - plot results of sensitivity analysis
#############

plotSA = function(x,
			contour.levels = NULL,
			zero.col = "red",
			lty.zero = 1,
			insig.col = "blue",
			lty.insig = 1,
			data.line = TRUE,
			X.pch = NULL,
			signif.level = 0.05,
			labcex = 0.75,
      limit.Xplot = F, #MH: limit plotting covariates to enlarge contour
			...) {
  #note in help: if contours are too rough, up nsim in sens fn
  Zcors = as.numeric(dimnames(x$alpha)[[2]]) #horizontal grids of U
  Ycors = as.numeric(dimnames(x$delta)[[1]]) #vertical grids of U
  xlab = "Alpha"
  ylab = "Delta"
  Xpart = x$Xcoef[!is.na(x$Xcoef[,1]) & !is.na(x$Xcoef[,2]),] #coefficients of null model.
  Xpart.plot = x$Xcoef.plot[!is.na(x$Xcoef.plot[,1]) & !is.na(x$Xcoef.plot[,2]),]
  Xpart.plot2 = cbind(Xpart.plot[,1],Xpart.plot[,2], ifelse(Xpart[,2]>=0,1,0)) #MH: add sign of coef of X on Y to Xpart  
  #note that due to correlation among Xs, some may not appear on plot
  #because observed partial cors don't map directly to coefs in this case
  #forcing inclusion can lead to difficult to read plot.  
  if(is.null(X.pch)){
    X.pch = ifelse(Xpart[,2]>=0,3,6) # plus sign for non-transformed plots, reverse triangle for transformed plots
  }
  
  taus = t(apply(x$tau, c(1,2), mean, na.rm = T))
  
  if(is.null(contour.levels)){
    exTau = c(taus[dim(taus)[1], dim(taus)[2]], taus[1,dim(taus)[2]]) #MH: extreme values of tau at right end
    clevels = round(seq(exTau[2]*.8, exTau[1]*.8, length.out = 14), 2) #MH: vals at which contours are drawn. round(seq(x$tau0*0.8, exTau*.8, length.out = 8),2)
  }else{
    clevels = contour.levels
  } 
  
  if (any(Xpart[,2]<0)) {
    cat("Note: Predictors with negative coefficients for the response surface have been transformed through multiplication by -1 and are displayed as inverted triangles.", "\n")    
  }

  par(mgp = c(2,.5,0)) #dist of axis label, tick mark label, tick mark
  
  if (limit.Xplot) {
    #old codes
    plot(Xpart.plot2[,1], Xpart.plot2[,2], col=c("red","blue")[as.factor(Xpart.plot2[,3])], xlim = c(min(Zcors, na.rm = T),max(Zcors, na.rm = T)), 
         ylim = c(min(Ycors, na.rm = T),max(Ycors, na.rm = T)), pch = X.pch, xlab = xlab, ylab = ylab,...)    
  } else {
    #MH: define max, min of plots
    xplot.min = ifelse(min(Zcors, na.rm = T)<min(Xpart.plot[,1]),min(Zcors, na.rm = T),min(Xpart.plot[,1]))
    xplot.max = ifelse(max(Zcors, na.rm = T)>max(Xpart.plot[,1]),max(Zcors, na.rm = T),max(Xpart.plot[,1]))
    yplot.max = ifelse(max(Ycors, na.rm = T)>max(Xpart.plot[,2]),max(Ycors, na.rm = T),max(Xpart.plot[,2]))  
    plot(Xpart.plot2[,1], Xpart.plot2[,2], col=c("red","blue")[as.factor(Xpart.plot2[,3])], xlim = c(xplot.min,xplot.max),
         ylim = c(0,yplot.max), pch = X.pch, xlab = xlab, ylab = ylab,...)
  }

  abline(h = 0)
  abline(v = 0)
  
  contour(Zcors, Ycors, taus, levels = clevels, 
          add = T, labcex = labcex,...)
  
  contour(Zcors, Ycors, taus, levels = 0, 
          add = T, col = zero.col,lty = lty.zero,labcex = labcex,...)
  
  contour(Zcors, Ycors, taus/t(apply(x$se.tau, c(1,2), mean)), labels = "N.S.",
          levels = -sign(x$tau0)*qnorm(signif.level/2), add = T, col = insig.col,
          lty = lty.insig, labcex = labcex,...)
  
  legend(0.8*max(Zcors), 0, legend = round(x$tau0,2), cex = labcex,
         yjust = 0.5, x.intersp = 0,
         bg = ifelse(par("bg")== "transparent", "white", par("bg")), box.lty = 0)
  
  if (all(sign(Xpart[,1])!=sign(x$tau0))){
    warning("Cannot add data line because XXXXX.")
  }else{
    if(data.line & length(Xpart)>1){
      proj.pts = apply(Xpart, 1, mean)
      max.pt = Xpart[proj.pts == max(proj.pts[sign(Xpart[,1])==sign(x$tau0)]),]
      zcor = (1:length(Zcors))[abs(Zcors-max.pt[1]) ==  min(abs(Zcors-max.pt[1]))]
      if((Zcors[zcor] > max.pt[1] & zcor > 1)||(zcor==length(Zcors))){ #MH: ||(zcor==length(Zcors)) is added to avoid error
        zpts = c(zcor-1, zcor)
      }else{
        zpts = c(zcor, zcor+1)
      }
      ycor = (1:length(Ycors))[abs(Ycors-max.pt[2]) ==  min(abs(Ycors-max.pt[2]))]
      if((Ycors[ycor] > max.pt[2] & ycor > 1)||(ycor==length(Ycors))){ #MH: ||(ycor==length(Ycors)) is added to avoid error
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
    }else{
      if(data.line)
        warning("Cannot add data line because there are no non-treatment covariates.")
    }
  }
} #end of plotSA

############
#plotSA.cont
############
plotSA.cont = function(x, coef.axes = F,
			signif.level = 0.05,
			labcex = 0.75,
			nullCI = TRUE,
			...) {
require(lattice)

	trt.ests <- sd.ests <- LCI <- UCI <- LCI0 <- UCI0 <- Z <- tau0 <- cell <- NULL
	ct = 0
	for(i in 1:dim(x$tau)[1]) {
	for(j in 1:dim(x$tau)[2]) {
		ct = ct+1
		trt.temp = apply(x$tau[i,j,,], 2, mean)
		sd.temp = apply(x$se.tau[i,j,,],2,mean)

		trt.ests = c(trt.ests, trt.temp)
		sd.ests = c(sd.ests, sd.temp)

		LCI = c(LCI,trt.temp + qnorm(signif.level/2)*sd.temp)
		UCI = c(UCI,trt.temp + qnorm(1-signif.level/2)*sd.temp)

		tau0 = c(tau0, x$tau0)
		LCI0 = c(LCI0, loess(x$tau0 + qnorm(signif.level/2)*x$se.tau0~x$Z)$fitted)
		UCI0 = c(UCI0, loess(x$tau0 + qnorm(1-signif.level/2)*x$se.tau0~x$Z)$fitted)

			
		Z = c(Z, x$Z)

		trt.cor = round(apply(x$trt.cor, 2, mean, na.rm = T),2)
		resp.cor = round(apply(x$resp.cor, 1, mean, na.rm = T),2)
		cell = c(cell, rep(paste("(", trt.cor[i], ",", resp.cor[j],")",sep = ""),length(x$Z)))
	}} 

	panel.bands <- function(x,y,upper,lower,subscripts,...,font,fontface) {
		cat("length x", length(x), "\n")
		cat("length lower", length(lower), "\n")
		resort <- order(x)
		upper <- upper[subscripts]
		lower <- lower[subscripts]
		panel.polygon(c(x[resort],rev(x[resort])), c(upper[resort], rev(lower[resort])),border = NA,...)
	}

	xyplot(tau0 + trt.ests + LCI + UCI ~ Z | cell, xlab = "Treatment", ylab = "Treatment effect", 
		layout = dim(x$tau)[c(1,2)], ylim = c(.95*min(c(LCI, x$tau0)), 1.05*max(c(UCI, x$tau0))),
		type = "smooth", col.line = c("white", "black", "blue", "blue"),
		sub = "(partial r^2 treatment, partial r^2 response)", cex.sub = 0.5,
		upper = UCI0, 
		lower = LCI0,
		col = "gray80",
		panel = function(x,y,upper,lower,...){
			resort <- order(x)
			if(nullCI) {
				panel.polygon(c(x[resort],rev(x[resort])), c(upper[resort], rev(lower[resort])),border = NA,...)
			}
			panel.xyplot(x, y,...)
		}) 
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

	tau0 = x$tau0
	LCI = tau0 + qnorm(signif.level/2)*x$se.tau0
	UCI = tau0 + qnorm(1-signif.level/2)*x$se.tau0

	line.plots <- xyplot(tau0 + LCI + UCI ~ x$Z, xlab = "Treatment", ylab = "Treatment effect", 
		#ylim = c(.95*min(c(LCI, x$tau0)), 1.05*max(c(UCI, x$tau0))),
		type = "smooth", col = c("black", "blue", "blue"),
		span = .75,
		scales = list(tck = c(1,0))
		) 
	hist.plot <- histogram(~ x$Z, freq = F, type = "density", col = "white",
			ylim = c(0, 4*max(hist(x$Z,plot=F)$density)))

	line.plots + as.layer(hist.plot, y.same = FALSE, axes = NULL, opposite = FALSE)
}




