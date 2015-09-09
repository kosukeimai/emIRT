// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>
//#include <RcppTN.h>
#include "vtn1.h"

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getVzstar_ordIRT(const arma::mat &Edd,
                    const arma::mat &beta,
                    const arma::mat &x,
                    const arma::mat &Etau,
                    const arma::mat &y,
                    const int N,
                    const int J
                    ) {

  double q1;
  signed int i, j;

  arma::mat Vzstar(N,J);

#pragma omp parallel for private(i,j,q1)
  for(i=0; i < N; i++){

	for(j=0; j < J; j++){

    q1 = x(i, 0) * beta(j, 0) + Etau(j,0);

//    if(y(i,j)==1)     Vzstar(i,j) = RcppTN::vtn1(q1, 1/sqrt(Edd(j,0)), R_NegInf, 0);
//    if(y(i,j)==2)     Vzstar(i,j) = RcppTN::vtn1(q1, 1/sqrt(Edd(j,0)), 0, 1);
//    if(y(i,j)==3)     Vzstar(i,j) = RcppTN::vtn1(q1, 1/sqrt(Edd(j,0)), 1, R_PosInf);
//    if(y(i,j)==0)     Vzstar(i,j) = RcppTN::vtn1(q1, 1/sqrt(Edd(j,0)), R_NegInf, R_PosInf);

    if(y(i,j)==1)     Vzstar(i,j) = vtn1(q1, 1/sqrt(Edd(j,0)), R_NegInf, 0);
    if(y(i,j)==2)     Vzstar(i,j) = vtn1(q1, 1/sqrt(Edd(j,0)), 0, 1);
    if(y(i,j)==3)     Vzstar(i,j) = vtn1(q1, 1/sqrt(Edd(j,0)), 1, R_PosInf);
    if(y(i,j)==0)     Vzstar(i,j) = vtn1(q1, 1/sqrt(Edd(j,0)), R_NegInf, R_PosInf);

	}
  }

  return(Vzstar) ;
} 



