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

#' Plot studentized residuals from the Bayesian Deming regression
#'
#' @export
#' @param bdpreg bdpreg object created with bdpreg
#' @return no return
#'

bdpPlotResiduals <- function(bdpreg){

  extr <- bdpExtract(bdpreg)

  dat<-bdpreg$standata

  if(dat$heteroscedastic == "linear") {
    d.text <- "Heteroscedastic linear model"
  }else{
      d.text <- "Homoscedastic linear model"
    }

  ymin <- qt(0.0001,df=dat$N-2)
  ymax <- qt(0.9999,df=dat$N-2)

  ymin <- if(ymin > min(extr$OptStandardRes) ) {
    ymin <-min(extr$OptStandardRes) } else{
      ymin <- ymin
    }

  ymax <- if(ymax < max(extr$OptStandardRes) ) {
    ymax <- max(extr$OptStandardRes) } else{
      ymax <- ymax
    }

  plot(extr$avgXY,extr$OptStandardRes, xlab="avgXY", ylab="Studentized residuals",
       main="Standardized residuals - TA", ylim=c(ymin,ymax))

  abline(h=0, lty=1)

  abline(h=qt(c(0.05,0.95),df=dat$N-2),lty=2)

  abline(h=qt(c(0.01,0.99),df=dat$N-2),col="blue",lty=3)

  abline(h=qt(c(0.001,0.999),df=dat$N-2),col="purple",lty=4)

  mtext(paste0(d.text),
        side=3, line=-1,adj=0.1,font=1)

  legend("topright",legend=c("95%","99%","99.9%"),col=c("black","blue","purple"),lty=2:4)

}
