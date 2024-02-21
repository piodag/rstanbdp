###############################################################################
##
## Copyright: Giorgio Pioda, 2024
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the license, or
## any later version.
##
## This software is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY, without even the implied
## warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the general R package license for details.
##
###############################################################################

#' Plot Bayesian Deming regression and the confidence intervals from the full posterior distribution
#'
#' @export
#' @param bdpreg bdpreg object created with bdpreg
#' @param ci Probability for the HDI credibility interval. Default 0.95
#' @return no return
#'

bdpPlot <- function(bdpreg,ci=0.95){

  stanRegr <- bdpreg$out
  dat <- bdpreg$standata

  extr_r <- rstan::extract(stanRegr,pars=c("intercept","slope"))
  coef.ab<-rstan::summary(stanRegr)$summary[,1]

  if (dat$heteroscedastic == "linear") {
    heteroscedastic.text <- "Linear heteroscedastic model with n="
  }else{
    heteroscedastic.text <- "Homoscedastic model with n="
  }

  #data <- cbind(X,Y)
  #dat <- as.data.frame(na.omit(data))

  plot(dat$X,dat$Y,xlab="X",ylab="Y")

  abline(a=coef.ab["intercept"],b=coef.ab["slope"],col="blue",lwd=2)
  abline(a=0,b=1,col="red",lty=2)

  minX<-min(dat$X)-IQR(dat$X)/10
  maxX<-max(dat$X)+IQR(dat$X)/10

  xvall<-data.frame(X=seq(minX,maxX,length.out=100))

  #pred_h<-t(apply(xvall,1,function(x) quantile(extr_r$intercept + x * extr_r$slope , c(0.025,0.975))))
  #pred_h<-t(apply(xvall,1, function(x) HDInterval::hdi(extr_r$intercept + x * extr_r$slope)))
  #pred_m<-t(apply(xvall,1, function(x) median(extr_r$intercept + x * extr_r$slope)))
  pred_h<-t(apply(xvall,1, function(x) bayestestR::hdi(extr_r$intercept + x * extr_r$slope,ci=ci)))
  pred_h<-matrix(unlist(pred_h),ncol=3,byrow = T)

  polygon(c(xvall$X,rev(xvall$X)),c(pred_h[,3],rev(pred_h[,2])), col=rgb(0,0,1,0.1), border=NA)

  #lines(xvall,pred_m,col="blue",lty=2)


  #lines(xvall$X,pred_h[,1],col="blue",lty=2)
  #lines(xvall$X,pred_h[,2],col="blue",lty=2)

  legend("bottomright",legend=c("Regression","Identity"),
         lty=c(1,2),lwd=c(2,1),col=c("blue","red"))

  mtext(paste0(heteroscedastic.text,dat$N," data points. \n y = ",signif(coef.ab["slope"],5),"*x",ifelse(coef.ab["intercept"] > 0,"+","-"),
              abs(signif(coef.ab["intercept"],5))),
        side=3, line=-2,adj=0.1,font=1)

  mtext(paste0(signif(ci,3)*100,"% HDI-CI simulated from the full Bayesian pairs distribution"),
        side=1, line=2,
        #adj=0.1,
        font=1)

  title(paste("Plot of the Bayesian Deming regression"))

}
