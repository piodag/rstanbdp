///////////////////////////////////////////////////////////////////////////////
//
// Copyright: Giorgio Pioda, 2024
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the license, or
// any later version.
//
// This software is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY, without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the general R package license for details.
//
//////////////////////////////////////////////////////////////////////7////////
//
// The input data are vectors 'X' and 'Y' of length 'N'.
// df must also be provided
// ErrorRatio accounts for Deming's fixed variance ratio
// Slope truncated prior at 0.3333. Less generalistic regression
// Stable if ErrorRatio is much different from 1

data {
  int <lower=0> N; // nr of observ
  real <lower=0> ErrorRatio; // Deming error ratio
  int<lower=0> df; // respose variable size
  vector[N] X; // predictor vector
  vector[N] Y;  // response vector
  vector[N] avgXY; // Inverse variance weighted average

   // Priors tweaking
  real slopeMu;
  real slopeSigma;
  real slopeTruncMin;
  real slopeTruncMax;

  real interceptMu;
  real interceptSigma;

  real AlphaMu;
  real AlphaSigma;

  real BetaMu;
  real BetaSigma;
  real BetaTruncMin;
  real BetaTruncMax;


}


transformed data {

}


parameters {
  real intercept; // intercept
  real slope; // coef for predictors
  real <lower=0> Alpha;
  real Beta;
}

transformed parameters {


}

model {
  vector[N] Pred;
  vector[N] Dis;
  vector[N] HatX;
  vector[N] HatY;
  vector[N] Opti;


     Pred = X * slope + intercept;
     Dis = Y-Pred;
     HatX = X+(ErrorRatio*Dis*slope/(1+ErrorRatio*slope^2));
     HatY = Y-(Dis/(1+ErrorRatio*slope^2));

     for (n in 1:N) {

     Opti[n] = sqrt((X[n]-HatX[n])^2+(Y[n]-HatY[n])^2);

}


  intercept ~ normal(interceptMu,interceptSigma); // flexible prior for intercept
  slope ~ normal(slopeMu,slopeSigma) T[slopeTruncMin,slopeTruncMax]; // truncated normal prior to avoid dummy convergence to zero
  Alpha ~ normal(AlphaMu,AlphaSigma) T[0,]; // half normal prior for the heteroscedastic variance (intercept)
  Beta ~ normal(BetaMu,BetaSigma) T[BetaTruncMin,BetaTruncMax]; // prior slope for the heteroscedastic (variance slope)

  Opti ~ student_t(df,0, Alpha * exp(Beta * avgXY));

}

