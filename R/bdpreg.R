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

#' Bayesian Deming Pioda Regression for two method comparison with Rstan
#'
#' \code{bdpreg} is used to compare two measurement methods by means of a Bayesian regression analysis.
#'
#' The Bayesian Deming regression can be run in a traditional fashion. In this case the error term is
#' sampled from a \eqn{T} distribution with \eqn{N-2} degree of freedom (\eqn{N} sample size).
#'
#' The Bayesian Deming regression can be run as a robust regression specifying a decreased \eqn{df} parameter.
#' It is possible to set \eqn{df = 1} and perform the sampling from an extremely robust Cauchy distribution
#' to suppress leveraged outliers. For moderate robustness a reasonably low value of \eqn{df} in the interval
#' \eqn{[6;10]} can be an appropriated choice.
#'
#' \code{ErrorRatio} can be set as usual for classical Deming regression. Default is 1. Strong
#' \code{ErrorRatio} can lead to instability in the chains that may not converge after the burn in. For
#' this purpose the \code{trunc} parameter can be used. In this way the normal distribution for the
#' slope gets truncated at a minimum of 0.3333 (default). The parameter \code{slopeTruncMin}
#' can override this value.
#'
#' With the parameter \code{heteroschedastic} it is possible to use an alternative regression which
#' models the heteroscedasticity with a linear growing variance. \code{Alpha} and \code{Beta} are the
#' intercept and the slope for the variance variation. \code{Alpha} must be > 0. \code{Beta} is usually
#' zero if no real heteroscedasticity is detected. Alternatively \code{Beta} shows low positive values,
#' typically below 0.5 if heteroscedasticity is successfully modeled. The CI of \code{Beta} could indeed
#' act as a test for heteroscedasticity. According to these empiric observations, \code{Beta} is also
#' truncated to avoid erratic behavior of the Hamiltonian sampler.
#'
#' The Bayesian Deming regression is recommended in many cases where traditional and non parametric
#' method fail. It is particularly convenient with very small data set and/or with data set with low
#' digit precision. In fact Bayesian Deming regression has no problem with ties.
#'
#' The method with linear heteroscedastic fitting can be a meaningful answer to heteroscedastic
#' data set. The CI are much narrower and the trade off between robustness and power can find
#' a natural solution. It must be considered as highly experimental but also highly promising
#' method. Users are advised to carefully check the sampled output for undesirable correlation
#' between Alpha and/or Beta vs the slope and/or intercept. A plot with \code{pairs()} highly
#' recommended.
#'
#' Stan is usually good enough that init values for the chains must not be specified. In extreme cases
#' it is anyway possible to set init values as a list of list.
#'
#'
#' @export
#' @param X Numeric vector of input values.
#' @param Y Numeric vector of output values.
#' @param trunc Boolean. Default TRUE. Use truncated slope prior for stability with extreme ErrorRatios.
#' See \code{slopeTruncMin}.
#' @param ErrorRatio Deming variance ratio between reference and test method. Default = 1.
#' @param df Degree of freedom. Must be df >= 1 (robust Cauchy regression). Default is \eqn{N-2}, For robust
#' regression set it to \eqn{df < N-2}
#' @param heteroscedastic Bayesian Deming model choice. Alternatives are:
#'          \code{"homo"} - Homoscedastic model. Default.\cr
#'          \code{"linear"} - Heteroscedastic with linear growth of the variance. Highly experimental model.\cr
#'          \code{"exponential"} - Heteroscedastic with exponential growth of the variance. Highly experimental model.\cr
#' @param slopeMu Slope normal Mu prior value. Default 1.
#' @param slopeSigma Slope normal Sigma prior value. Default 0.3.
#' @param slopeTruncMin slope normal lower truncation limit. Default 0.3333.
#' @param slopeTruncMax slope normal higher truncation limit. Default 10.
#' @param interceptMu Intercept normal Mu prior value. Default 0.
#' @param interceptSigma Intercept normal Sigma prior value. Default 30.
#' @param sigmaLambda sigma exponential prior lambda. Default 0.3.
#' @param AlphaMu Lin. heterosc. intercept normal mu prior. Must be > 0. Default 1.
#' @param AlphaSigma Lin. heterosc. intercept normal sigma prior. Default 10.
#' @param BetaMu Lin. heterosc. slope normal prior. Default 0.1.
#' @param BetaSigma Lin. heterosc. slope normal prior. Default 0.5.
#' @param BetaTruncMin Lin. heterosc. slope normal prior truncation min. Default -1.
#' @param BetaTruncMax  Lin. heterosc. slope normal prior truncation min. Default 1.
#' @param ... Arguments passed to `rstan::sampling` (e.g. `iter`, `chains`)
#' @return An object of class `bdpreg`  which contains `out` a `stanfit` object returned by `rstan::sampling` and `standata` as list of input parameters.
#' @references G. Pioda (2014) <https://piodag.github.io/bd1/>
#' @examples
#'library(rstanbdp)
#'data(glycHem)
#'
#'# Bayesian Deming Regression, for example with  df=10
#'fit.1 <-bdpreg(glycHem$Method1,glycHem$Method2,heteroscedastic="homo",
#'               df=10,chain=1,iter=1000)
#'
#'# Print results
#'bdpPrint(fit.1,digits_summary = 4)
#'
#'# Plot 2D intercepts /slopes pairs with CI and MD distance
#'bdpPlotBE(fit.1,cov.method="MCD",ci=0.95)
#'
#'# Plot regression with CI
#'bdpPlot(fit.1,ci=0.95)
#'
#'# Calculate response, plot histogram and CI
#'bdpCalcResponse(fit.1,Xval = 6)
#'
#'# Extract Xhat, Yhat and Residuals
#'bdpExtract(fit.1)
#'
#'# Plot a traceplot of the sampled chains
#'bdpTraceplot(fit.1)
#'
#'# Plot standardized residuals
#'bdpPlotResiduals(fit.1)
#'
#'# Plot posterior samples pairwise
#'bdpPairs(fit.1)

bdpreg <- function(X, Y, ErrorRatio = 1, df = NULL, trunc = TRUE,
                   heteroscedastic = c("homo","linear"),
                   slopeMu = 1, slopeSigma = 0.3,
                   slopeTruncMin = 0.3333, slopeTruncMax = 10,
                   interceptMu = 0, interceptSigma = 30,
                   sigmaLambda = 0.3,
                   AlphaMu = 1, AlphaSigma = 10,
                   BetaMu = 0.1, BetaSigma = 0.5,
                   BetaTruncMin = -1, BetaTruncMax = 1,...) {

  stopifnot(length(X)==length(Y))

  # Dirty workaround to have heteroscedastic="homo" as default

  if (missing(heteroscedastic) ) {
    heteroscedastic <- "homo"
  }

  if(heteroscedastic %in% c("homo","linear","exponential")){
       heteroscedastic <- heteroscedastic
  }else{ heteroscedastic <-"homo"}

  data <- cbind(X,Y)
  dat <- na.omit(data)
  if(nrow(data)>nrow(dat)){warning(paste("Some NA removed, data set contains",nrow(dat),"valid pairs"))}

  stopifnot(ErrorRatio > 0)
  if(!is.null(df)){df = df} else {df = nrow(dat) - 2 }
  stopifnot(df >= 1)

  avgXY <- (dat[,1] +  ErrorRatio * dat[,2]) / (1 + ErrorRatio)

  standata <- list(X = dat[,1], Y = dat[,2], avgXY = avgXY, N = nrow(dat), df = df, trunc = trunc,
                   ErrorRatio = ErrorRatio, heteroscedastic = heteroscedastic,
                   slopeMu = slopeMu, slopeSigma = slopeSigma,
                   slopeTruncMin = slopeTruncMin, slopeTruncMax = slopeTruncMax,
                   interceptMu = interceptMu, interceptSigma = interceptSigma,
                   sigmaLambda = sigmaLambda,
                   AlphaMu = AlphaMu, AlphaSigma = AlphaSigma,
                   BetaMu = BetaMu, BetaSigma = BetaSigma,
                   BetaTruncMin = BetaTruncMin, BetaTruncMax = BetaTruncMax)

  if (heteroscedastic == "linear") {
    if (trunc == TRUE){
      out <- rstan::sampling(stanmodels$bdpreg_linhettrunc, data = standata, ...)
    } else {
      out <- rstan::sampling(stanmodels$bdpreg_linhet, data = standata, ...)
      }
} else if (heteroscedastic == "exponential"){
    out <- rstan::sampling(stanmodels$bdpreg_exphettrunc, data = standata, ...)
} else {
  if (trunc == TRUE){
    out <- rstan::sampling(stanmodels$bdpreg_homotrunc, data = standata, ...)
  }else{
    out <- rstan::sampling(stanmodels$bdpreg_homo, data = standata, ...)
  }
}

  ret<-list(out=out,standata=standata)
  attr(ret,"class") <- "bdpreg"
  return(ret)
}
