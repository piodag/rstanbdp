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

#' Extract relevant data from the Bayesian Deming regression
#'
#' @export
#' @param bdpreg bdpreg object created with bdpreg
#' @return A data frame with extracted data
#'

bdpExtract <- function(bdpreg){

  stanRegr <- bdpreg$out
  dat <- bdpreg$standata

  extr_r <- rstan::extract(stanRegr,pars=c("intercept","slope"))
  coef.ab<-rstan::summary(stanRegr)$summary[,1]

  mX <- mean(dat$X)
  mY <- mean(dat$Y)

  D <- dat$Y-(coef.ab[1]+coef.ab[2]*dat$X)
  Xhat <- dat$X+(dat$ErrorRatio*coef.ab[2]*D/(1+dat$ErrorRatio*coef.ab[2]^2))
  Yhat <- dat$Y-(D/(1+dat$ErrorRatio*coef.ab[2]^2))

  Xres <- dat$X-Xhat
  Yres <- dat$Y-Yhat

  avgXY = dat$avgXY

  OptRes <- sqrt(Xres^2+Yres^2)*sign(Yres)

  if (dat$heteroscedastic == "linear") {
    Sigma <- coef.ab[3]+coef.ab[4]*avgXY

    return(data.frame(X = dat$X, Y = dat$Y, avgXY = (dat$X+dat$Y)/2, diffXY = dat$Y-dat$X,
                      Xhat=Xhat, Yhat=Yhat, Xres=Xres, Yres=Yres, OptRes = OptRes, linSigma = Sigma,
                      OptStandardRes = OptRes/Sigma))

  } else if (dat$heteroscedastic == "exponential"){

    Sigma <- coef.ab[3]*exp(coef.ab[4]*avgXY)

    return(data.frame(X = dat$X, Y = dat$Y, avgXY = (dat$X+dat$Y)/2, diffXY = dat$Y-dat$X,
                      Xhat=Xhat, Yhat=Yhat, Xres=Xres, Yres=Yres, OptRes = OptRes, linSigma = Sigma,
                      OptStandardRes = OptRes/Sigma))
  } else {

    return( data.frame(X = dat$X, Y = dat$Y, avgXY = avgXY, diffXY = dat$Y-dat$X,
                       Xhat=Xhat, Yhat=Yhat, Xres=Xres, Yres=Yres, OptRes = OptRes, OptStandardRes = OptRes/coef.ab[3]))
  }

}
