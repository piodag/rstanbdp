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

#' Plot the calculated Y response with CI from the full Bayesian posterior distribution
#'
#' @export
#' @param bdpreg bdpreg object
#' @param Xval Reference method data
#' @param ci Probability for the HDI credibility interval. Default 0.95.
#' @param ... Arguments passed to `hist` (e.g. `breaks`, `xlim`, ...).
#' @return no return
#'

bdpCalcResponse<-function(bdpreg,Xval,ci=0.95,...){

  stanRegr <- bdpreg$out
  dat <- bdpreg$standata


  extr_r <- rstan::extract(stanRegr,pars=c("intercept","slope"))
  coef.ab<-rstan::summary(stanRegr)$summary[,1]

  if (dat$heteroscedastic == "linear") {
    het.text <- "Linear heteroscedastic model"
  } else if (dat$heteroscedastic == "exponential") {
    het.text <- "Exponential heteroscedastic model"
  } else {
    het.text <- "Homoscedastic model"
  }


  pred_r <- extr_r$intercept + Xval * extr_r$slope

  iqr <- IQR(pred_r)

  xmin <- ifelse(Xval < min(pred_r),Xval - iqr/3, min(pred_r) - iqr/3)
  xmax <- ifelse(Xval > max(pred_r),Xval + iqr/3, max(pred_r) + iqr/3)

  #ci.hdi<-HDInterval::hdi(pred_r)
  ci.hdi<-bayestestR::hdi(pred_r,ci=ci)
  ci.hdi<-matrix(unlist(ci.hdi),ncol=3,byrow = T)

  pred.h <- hist(pred_r,xlim=c(xmin,xmax),main="",xlab="predicted Y",...)

  abline(v=rstan::summary(pred_r)[4],col="red",lty=2)
  abline(v=ci.hdi[2],col="red",lty=2)
  abline(v=ci.hdi[3],col="red",lty=2)
  abline(v=Xval, col="black",lty=2)

  legend("topright",legend=c("Res, HDI-CI","X value"),
         lty=c(2,2),lwd=c(1,1),col=c("red","black"),
         title=paste0("n = ",dat$N,", d.f. = ",dat$df),
         cex=1)

  mtext(paste0(het.text,"\n Res Y = ",signif(rstan::summary(pred_r)[4],5),"\n HDI-CI ",signif(ci,3)*100,"% \n [ ",signif(ci.hdi[2],5)," ; ",signif(ci.hdi[3],5)," ]"),
        side=3, line=-3,adj=0.02,cex=1)

  mtext(paste("Y predictions simulated from the full Bayesian pairs distribution"),
        side=1, line=2,
        #adj=0.1,
        cex=1)

  title(paste("Histogram of the full response for X =",Xval),cex.main=1)

}
