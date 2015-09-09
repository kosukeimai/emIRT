// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>
//#include <RcppTN.h>
#include "etn1.h"

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEystar_dynIRT(arma::mat &Eystar,
					const arma::mat &alpha,
                    const arma::mat &beta,
                    const arma::mat &x,
                    const arma::mat &y,
                    const arma::mat &bill_session,
                    const arma::mat &startlegis,
                    const arma::mat &endlegis,
                    const int N,
                    const int J
                    ) {

	double q1;
	signed int i, j;

    // Main Calculation
#pragma omp parallel for private(i,j,q1)
  	for(i=0; i < N; i++){

		for(j=0; j < J; j++){

			if( (bill_session(j,0) <= endlegis(i,0)) && (bill_session(j,0) >= startlegis(i,0)) ){

		    	q1 = alpha(j,0) + x(i,bill_session(j,0)) * beta(j,0);

//			    if(y(i,j)==1)     Eystar(i,j) = RcppTN::etn1(q1, 1.0, 0.0, R_PosInf);
//			    if(y(i,j)==-1)    Eystar(i,j) = RcppTN::etn1(q1, 1.0, R_NegInf, 0.0);
//			    if(y(i,j)==0)     Eystar(i,j) = RcppTN::etn1(q1, 1.0, R_NegInf, R_PosInf);

			    if(y(i,j)==1)     Eystar(i,j) = etn1(q1, 1.0, 0.0, R_PosInf);
			    if(y(i,j)==-1)    Eystar(i,j) = etn1(q1, 1.0, R_NegInf, 0.0);
			    if(y(i,j)==0)     Eystar(i,j) = etn1(q1, 1.0, R_NegInf, R_PosInf);

				// Note: Taking etn() of extreme truncated normals is a problem			    
			    // etn(-9.49378,1,0,1000) gives Inf for Eystar, which crashes everything
			    // In these cases, we should ignore the vote
			    if( !(arma::is_finite(Eystar(i,j))) ) Eystar(i,j) = q1;
		    
			} // end if( (bill_session(j,0) <= endlegis(i,0)) 

		}  // end for(j=0; j < J; j++)
  	}  // end for(i=0; i < N; i++)

    return; 
} 
