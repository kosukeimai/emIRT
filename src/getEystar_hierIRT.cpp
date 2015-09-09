// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>
//#include <RcppTN.h>
#include "etn1.h"

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEystar_hierIRT(arma::mat &Eystar,
					const arma::mat &y,
					const arma::mat &z,
                    const arma::mat &g,
                    const arma::mat &i,
                    const arma::mat &j,
                    const arma::mat &Ea,
                    const arma::mat &Eb,
                    const arma::mat &Egamma,
                    const arma::mat &Eta,
                    const int ND,
                    const int NG,
                    const int NI,
                    const int NJ,
                    const int NL
                    ) {

	double ml;
	signed int l;

    // Main Calculation
#pragma omp parallel for
  	for(l=0; l < NL; l++){

  		ml = Ea(j(l,0),0) + Eb(j(l,0),0)*accu(Egamma.row(g(i(l,0),0)) % z.row(i(l,0))) + Eta(i(l,0),0)*Eb(j(l,0),0);

//		if(y(l,0)==1)     Eystar(l,0) = RcppTN::etn1(ml, 1.0, 0.0, R_PosInf);
//		if(y(l,0)==-1)    Eystar(l,0) = RcppTN::etn1(ml, 1.0, R_NegInf, 0.0);
//		if(y(l,0)==0)     Eystar(l,0) = RcppTN::etn1(ml, 1.0, R_NegInf, R_PosInf);

		if(y(l,0)==1)     Eystar(l,0) = etn1(ml, 1.0, 0.0, R_PosInf);
		if(y(l,0)==-1)    Eystar(l,0) = etn1(ml, 1.0, R_NegInf, 0.0);
		if(y(l,0)==0)     Eystar(l,0) = etn1(ml, 1.0, R_NegInf, R_PosInf);

		// Note: Taking etn() of extreme truncated normals is a problem			    
		// etn(-9.49378,1,0,1000) gives Inf for Eystar, which crashes everything
		// In these cases, we should ignore the vote
		if( !(arma::is_finite(Eystar(l,0))) ) Eystar(l,0) = ml;

			    
  	}  // end for(l=0; l < NL; l++)

    return; 
} 
