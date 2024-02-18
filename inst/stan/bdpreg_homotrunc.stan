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
// Bayesian Deming regression with homoscedastic variance
//
// The input data are vectors 'X' and 'Y' of length 'N'.
// df must also be provided
// ErrorRatio accounts for Deming's fixed variance ratio
// slope truncated prior at 0.3333. Less generalistic regression
// Stable if ErrorRatio is much different from 1

data {
  int<lower=0> N; // nr of observ
  int<lower=0> df; // respose variable size
  real<lower=0> ErrorRatio;
  vector[N] X; // X vector
  vector[N] Y;  // Y vector

  // Priors tweaking
  real slopeMu;
  real slopeSigma;
  real slopeTruncMin;
  real slopeTruncMax;

  real interceptMu;
  real interceptSigma;

  real sigmaLambda;

}


transformed data {

}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.

parameters {
  real intercept; // intercept
  real slope; //  slope
  real<lower=0> sigma; // sigma


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



  intercept ~ normal(interceptMu,interceptSigma);
  slope ~ normal(slopeMu,slopeSigma) T[slopeTruncMin,slopeTruncMax];
  sigma ~ exponential(sigmaLambda);

  Opti ~ student_t(df,0,sigma);

}

