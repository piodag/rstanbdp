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

#' Print summary of sampled data
#'
#' @export
#' @param bdpreg bdpreg object created with bdpreg
#' @param digits_summary number of digits for the results
#' @param ... Arguments passed to `rstan::print`
#' @return Print of the `stanfit` object
#'

bdpPrint <- function(bdpreg,digits_summary = 4,...){

  stanRegr <- bdpreg$out

  print(stanRegr,digits_summary = digits_summary,...)

}
