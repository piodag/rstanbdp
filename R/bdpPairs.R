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

#' Pairs plot of the full posterior predictors
#'
#' @export
#' @param bdpreg bdpreg object created with bdpreg
#' @param ... Arguments passed to `rstan::pairs`.
#' @return Pairs plot of the `stanfit` object
#'

bdpPairs <- function(bdpreg,...){

  stanRegr <- bdpreg$out

  pairs(stanRegr,...)

}
