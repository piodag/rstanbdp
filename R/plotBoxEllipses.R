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

#' Plot regression posterior pairs with CI Box and MD ellipses
#' @export
#' @param bdpreg bdpreg object created with bdpreg
#' @param cov.method rrcov covariance method ("SDe", "MCD", or "Classical"). Default MCD.
#' @param ci Probability for the HDI credibility interval. Default 0.95.
#' @param ... Arguments passed to `plot` (e.g. `xlim`, `xlab`, `main`)
#' @return no return

bdpPlotBE<-function(bdpreg,cov.method="MCD",ci=0.95,...){

  stanRegr <- bdpreg$out
  dat <- bdpreg$standata

  extr.pairs <- rstan::extract(stanRegr,pars=c("intercept","slope"))
  extr.pairs <- data.frame(B0=extr.pairs$intercept,B1=extr.pairs$slope)

  coef.ab<-rstan::summary(stanRegr)$summary[,1]

  if (dat$heteroscedastic == "linear") {
    het.text <- "Linear heteroscedastic model"
  } else if (dat$heteroscedastic == "exponential") {
    het.text <- "Exponential heteroscedastic model"
  } else {
    het.text <- "Homoscedastic model"
  }

    minAlfa<-ifelse(min(extr.pairs$B0) > -0.1, -0.1, min(extr.pairs$B0))
    maxAlfa<-ifelse(max(extr.pairs$B0) < 0.1, 0.1, max(extr.pairs$B0))

    minBeta<-ifelse(min(extr.pairs$B1) > 0.99, 0.99, min(extr.pairs$B1))
    maxBeta<-ifelse(max(extr.pairs$B1) < 1.01, 1.01, max(extr.pairs$B1))


    if(cov.method == "MCD") {
      t.mcd<-rrcov::CovMcd(cbind(extr.pairs$B0,extr.pairs$B1))
      text.label<-"robust MCD"
    } else if(cov.method == "Classical") {
      t.mcd<-rrcov::CovClassic(cbind(extr.pairs$B0,extr.pairs$B1))
      text.label<-"Classical"
    }else if ((cov.method == "SDe")){
      t.mcd<-rrcov::CovRobust(cbind(extr.pairs$B0,extr.pairs$B1))
      text.label<-"robust SDe"
    }

    t.md<-mahalanobis(c(0,1),center=t.mcd$center,cov=t.mcd$cov)
    res.p<-pchisq(t.md,df=2,lower.tail=F)
    names(res.p)<-"Chisq global p-value"
    plot(extr.pairs$B0,extr.pairs$B1,
         xlim=c(minAlfa,maxAlfa),
         ylim=c(minBeta,maxBeta),
         xlab="intercept",ylab="slope",
         pch=16,col=rgb(0,0,0,alpha=0.1),...)

    points(0,1,pch=4,cex=2,col="red",lwd=3)

    extr.pairs<-as.data.frame(extr.pairs)

    grid()

    res.cx<-bayestestR::hdi(extr.pairs[,1],ci=ci)
    res.cy<-bayestestR::hdi(extr.pairs[,2],ci=ci)

    rect(res.cx[2],res.cy[2],res.cx[3],res.cy[3],

    #res.cx<-extr.pairs@para[1,1:4]
    #res.cy<-extr.pairs@para[2,1:4]
    #rect(res.cx[3],res.cy[3],res.cx[4],res.cy[4],

         border="purple",lty=3,lwd=2)



    med.cx<-mean(extr.pairs[,1])
    med.cy<-mean(extr.pairs[,2])
    points(med.cx,med.cy,col="purple",pch=3,cex=2,lwd=3)

    mixtools::ellipse(mu=t.mcd$center,sigma=t.mcd$cov,alpha=0.05,
                      npoints=250,newplot=FALSE,
                      draw=TRUE, col="blue",lwd=1,lty=1)

    mixtools::ellipse(mu=t.mcd$center,sigma=t.mcd$cov,alpha=0.01,
                      npoints=250,newplot=FALSE,
                      draw=TRUE, col="blue",lwd=1,lty=2)

    points(t.mcd$center[1],t.mcd$center[2],col="blue",pch=8,cex=2,lwd=3)

    mtext(paste(cov.method,"center \n","intercept:",
                signif(t.mcd$center[1],5),"\n slope:",
                signif(t.mcd$center[2],5)),
          side=1, line=-3,
          adj=0.02,
          cex=1)

    mtext(paste("Chisq. p-value, 2 d.f.: ",signif(res.p,4)),side=1,
          line=-1,adj=0.02,cex=1)

    mtext(paste0(het.text,", HDi-CI ",signif(ci,3)*100,"%\n",
                "intercept: ",signif(as.numeric(med.cx),5)," [",signif(as.numeric(res.cx[2]),5) ,";",signif(as.numeric(res.cx[3]),5) ,"] \n",
                "slope: ",signif(as.numeric(med.cy),5)," [",signif(as.numeric(res.cy[2]),5),";",signif(as.numeric(res.cy[3]),5),"]"),
          side=3, line=-3,adj=0.02,cex=1)

    legend("topright",legend=c("e95%","e99%",paste0("HDI-CI ",signif(ci,3)*100,"%")),
           lty=c(1,2,3),lwd=c(1,1,2),col=c("blue","blue","purple"),
           title=paste0("n = ",dat$N, ", d.f = ",dat$df), cex=1)

    title(paste("Box & ellipses of the Bayesian samples with",text.label,"covariance"),
          cex.main=1)
}

#' Mahalanobis distance for the posterior pairs
#' @param stanRegr Rstan rstanbdp object
#' @param cov.method rrcov covariance method ("SDe", "MCD", or "Classical"). Default MCD.
#' @return Chi squared probability of the MD

bmpMD<-function(stanRegr,cov.method){

  extr.pairs <- rstan::extract(stanRegr,pars=c("intercept","slope"))
  extr.pairs <- data.frame(B0=extr.pairs$intercept,B1=extr.pairs$slope)

  if(cov.method == "MCD") {
    t.mcd<-rrcov::CovMcd(cbind(extr.pairs$B0,extr.pairs$B1))
    text.label<-"robust MCD"
  } else if(cov.method == "Classical") {
    t.mcd<-rrcov::CovClassic(cbind(extr.pairs$B0,extr.pairs$B1))
    text.label<-"Classical"
  }else if ((cov.method == "SDe")){
    t.mcd<-rrcov::CovRobust(cbind(extr.pairs$B0,extr.pairs$B1))
    text.label<-"robust SDe"
  }

  t.md<-mahalanobis(c(0,1),center=t.mcd$center,cov=t.mcd$cov)
  res.p<-pchisq(t.md,df=2,lower.tail=F)
  names(res.p)<-"Chisq global p-value"
  return(res.p)
}

