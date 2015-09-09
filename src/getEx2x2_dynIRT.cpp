// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEx2x2_dynIRT(arma::cube &Ex2x2,
				   const arma::mat &Ex,
                   const arma::mat &Vx,
                   const arma::mat &legis_by_session,
                   const arma::mat &Nlegis_session,
                   const int T
                   ) {

	int t;
	arma::mat Ex_sum;
	arma::mat Exx;

#pragma omp parallel for private(t,Ex_sum,Exx)
	for(t=0; t<T; t++){

		Ex2x2(0,0,t) = Nlegis_session(t,0);
		
		//Note that this accu() relies on assumption that Ex=0 for legislators not present that period
		Ex_sum = accu(Ex.col(t));
		Ex2x2(0,1,t) = Ex_sum(0,0);
		Ex2x2(1,0,t) = Ex_sum(0,0);

		//This also assumes that Vx=0 for all legislators not present that period
		Exx = trans(Ex.col(t)) * Ex.col(t) + sum(Vx.col(t),0);
		Ex2x2(1,1,t) = Exx(0,0);

	}

		//Rcout << "\n getEx2x2() complete..." << std::endl ;
		
    return;


/*
// Below is code for a version that does not assume missing values in Ex are zeroes

	int i,t;
	double Exx,Ex_sum;

    arma::cube Ex2x2(2, 2, T, arma::fill::zeros) ;

	for(t=0; t<T; t++){

		Ex2x2(0,0,t) = Nlegis_session(t,0);
		Ex_sum = 0.0;
		Exx = 0.0;

		//Non-matrix version, makes no assumption about missing values in Ex
		for(i=0; i < Nlegis_session(t,0); i++){
			Ex_sum += Ex(legis_by_session(t,i),t);
			Exx += Ex(legis_by_session(t,i),t)*Ex(legis_by_session(t,i),t) + Vx(legis_by_session(t,i),t);
		}

		Ex2x2(0,1,t) = Ex_sum;
		Ex2x2(1,0,t) = Ex_sum;
		Ex2x2(1,1,t) = Exx;

	}

    return(Ex2x2);
*/

}
