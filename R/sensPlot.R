#############
#sensPlot - plot results of sensitivity analysis
#############

sensPlot = function(x, 
                  contour.levels = NULL,
                  col.zero = "red",
                  lty.zero = 1,
                  col.insig = "blue",
                  lty.insig = 1,
                  data.line = TRUE,
                  X.pch = NULL,
                  signif.level = 0.05,
                  labcex = 0.75,
                  limit.Xplot = FALSE, #MH: limit plotting covariates to enlarge contour
                  txtlab = FALSE,  #add text label to the plots of covariates.
                  which.txtlab = NULL, #enter numeric vector to specify which label to show. e.g. c(1:3) shows first 3 covariates.
                  ...) {
  #note in help: if contours are too rough, up nsim in sens fn
  
  ##Add row/column for zeta = 0 if not included in grid
  null.tau=x$tau0
  null.se=x$se.tau0
  part.cors = x$sensParam == "cor"
 Zcors = as.numeric(dimnames(x$sp.z)[[2]]) #horizontal grids of U
 Ycors = as.numeric(dimnames(x$sp.y)[[1]]) #vertical grids of U
 if(part.cors){
	print.text <- "Partial correlations with U"
  }else {
	print.text <- "Coefficients on U" 
  }

  addYdim = sum(Ycors==0)==0
  Zcors = sort(c(Zcors,0))	#Note don't need if stmt for Z because main code excludes zeta.z = 0
  if(addYdim==TRUE){
    Ycors = c(0,Ycors)
  }
  
  Xpart = x$Xcoef[!is.na(x$Xcoef[,1]) & !is.na(x$Xcoef[,2]),] #coefficients of null model.  
  Xpart.plot = x$Xcoef.plot[!is.na(x$Xcoef.plot[,1]) & !is.na(x$Xcoef.plot[,2]),]
  Xpart.plot2 = cbind(Xpart.plot[,1],Xpart.plot[,2], ifelse(Xpart[,2]>=0,1,0)) #MH: add sign of coef of X on Y to Xpart  
  
  #note that due to correlation among Xs, some may not appear on plot
  #because observed partial cors don't map directly to coefs in this case
  #forcing inclusion can lead to difficult to read plot.  
  if(is.null(X.pch)){
    X.pch = ifelse(Xpart[,2]>=0,3,6) # plus sign for non-transformed plots, reverse triangle for transformed plots
  }else{
    X.pch = ifelse(Xpart[,2]>=0,X.pch[1],X.pch[2])	
  }
  
  nr = length(Zcors); nc = length(Ycors)
  taus.est = t(apply(x$tau, c(1,2), mean, na.rm = T))
  K = dim(x$se.tau)[3]
  W = apply(x$se.tau^2, c(1,2), mean, na.rm = T)
  B = apply(x$tau, c(1,2), sd, na.rm = T)^2
  se.est = t(sqrt(W+(1+1/K)*B))
  taus = matrix(null.tau,nr,nc)
  se.taus = matrix(null.se,nr,nc)
  row0 = c(1:length(Zcors))%*%((Zcors==0)*1) 
  if(addYdim==TRUE){
    taus[1:(row0-1),2:nc] = taus.est[1:(row0-1),]
    taus[(row0+1):length(Zcors),2:nc] = taus.est[(row0):nrow(se.est),]
    se.taus[1:(row0-1),2:nc] = se.est[1:(row0-1),]
    se.taus[(row0+1):length(Zcors),2:nc] = se.est[(row0):nrow(se.est),]
  }
  else{
    taus[1:(row0-1),] = taus.est[1:(row0-1),]
    taus[(row0+1):length(Zcors),] = taus.est[(row0):nrow(taus.est),]
    se.taus[1:(row0-1),] = se.est[1:(row0-1),]
    se.taus[(row0+1):length(Zcors),] = se.est[(row0):nrow(taus.est),]
  }
  dimnames(taus)[[1]] = Zcors
  dimnames(taus)[[2]] = Ycors
  taus <<- taus
  dimnames(se.taus)[[1]] = Zcors
  dimnames(se.taus)[[2]] = Ycors
  se.taus <<- se.taus
  
  if(is.null(contour.levels)){
    exTau = c(taus[dim(taus)[1], dim(taus)[2]], taus[1,dim(taus)[2]]) #extreme values of tau at right end
    clevels = round(seq(exTau[2]*.8, exTau[1]*.8, length.out = 14), 2) #vals at which contours are drawn
  }else{
    clevels = contour.levels
  } 
  
  if (any(Xpart[,2]<0)) {
    cat("Note: Predictors with negative coefficients for the response surface have been transformed through multiplication by -1 and are displayed as inverted triangles.", "\n")    
  }
  
  par(mgp = c(2,.5,0)) #dist of axis label, tick mark label, tick mark
  if(part.cors){
  	xlab = expression(paste("Partial cor. with U in model for treatment, ", rho^zu))
  	ylab = expression(paste("Partial cor. with U in model for response, ", rho^yu))
  }else{
  	xlab = expression(paste("Coef. on U in model for treatment, ", zeta^z))
  	ylab = expression(paste("Coef. on U in model for response, ", zeta^y))
  }

  if (limit.Xplot) {
    #old codes
    plot(Xpart.plot2[,1], Xpart.plot2[,2], col=c("red","blue")[as.factor(Xpart.plot2[,3])], xlim = c(min(Zcors, na.rm = T),max(Zcors, na.rm = T)), 
         ylim = c(min(Ycors, na.rm = T),max(Ycors, na.rm = T)), pch = X.pch, xlab = xlab, ylab = ylab)    
  } else {
    #MH: define max, min of plots
    xplot.min = ifelse(min(Zcors, na.rm = T)<min(Xpart.plot[,1]),min(Zcors, na.rm = T),min(Xpart.plot[,1]))
    xplot.max = ifelse(max(Zcors, na.rm = T)>max(Xpart.plot[,1]),max(Zcors, na.rm = T),max(Xpart.plot[,1]))
    yplot.max = ifelse(max(Ycors, na.rm = T)>max(Xpart.plot[,2]),max(Ycors, na.rm = T),max(Xpart.plot[,2]))  
    plot(Xpart.plot2[,1], Xpart.plot2[,2], col=c("red","blue")[as.factor(Xpart.plot2[,3])], xlim = c(xplot.min,xplot.max),
         ylim = c(0,yplot.max), pch = X.pch, xlab = xlab, ylab = ylab)
  }
  
  #codes for txtlab
  if (txtlab) {
    if (is.null(which.txtlab)) { #show all text label
      text(Xpart.plot2[,1], Xpart.plot2[,2], labels=x$varnames[-c(1,2)], cex=labcex, pos=1)
    } else { #show selected label
      which.txtlab2 = which.txtlab + 2
      varnames2 = x$varnames
      varnames2[as.numeric(paste(- which.txtlab2))] = ""
      text(Xpart.plot2[,1], Xpart.plot2[,2], labels=varnames2[-c(1,2)], cex=labcex, pos=1)
    }
  }
  
  abline(h = 0)
  abline(v = 0)
  
  legend(0.8*max(Zcors), 0, legend = round(x$tau0,2), cex = labcex,
         yjust = 0.5, x.intersp = 0, y.intersp = 0,
         bg = ifelse(par("bg")== "transparent", "white", par("bg")), 
         box.lty = 0)
  
  box()
  
  contour(Zcors, Ycors, taus, levels = clevels, 
          add = T, labcex = labcex, ...)
  
  contour(Zcors, Ycors, taus, levels = 0, lwd = 2,
          add = T, col = col.zero,lty = lty.zero,labcex = labcex,...)
  
  contour(Zcors, Ycors, taus/se.taus, labels = "N.S.",
          levels = c(-1, 1)*qnorm(signif.level/2), add = T, col = col.insig,
          lty = lty.insig, labcex = labcex, lwd = 2,...)
  
  
  if (all(sign(Xpart[,1])!=sign(x$tau0))){
    warning("Cannot add data line because XXXXX.")
  }else{
    if(data.line & length(Xpart)>1){
      proj.pts = apply(Xpart.plot, 1, mean)
      max.pt = Xpart.plot[proj.pts == max(proj.pts[sign(Xpart.plot[,1])==sign(x$tau0)]),]
      zcor = (1:length(Zcors))[abs(Zcors-max.pt[1]) ==  min(abs(Zcors-max.pt[1]))]
      if((Zcors[zcor] > max.pt[1] & zcor > 1)||(zcor==length(Zcors))){ 
        zpts = c(zcor-1, zcor)
      }else{
        zpts = c(zcor, zcor+1)
      }
      ycor = (1:length(Ycors))[abs(Ycors-max.pt[2]) ==  min(abs(Ycors-max.pt[2]))]
      if((Ycors[ycor] > max.pt[2] & ycor > 1)||(ycor==length(Ycors))){ 
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
              add = T, col = "grey",labcex = labcex, lwd = 2,...)
    }else{
      if(data.line)
        warning("Cannot add data line because there are no non-treatment covariates.")
    }
  }
} #end of sensPlot



