// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEta_hierIRT(arma::mat &Eta,
					arma::mat &Eta2,
					const arma::mat &Veta,
					const arma::mat &Eystar,
					const arma::mat &Eb,
					const arma::mat &Eba,
					const arma::mat &Ebb,
					const arma::mat &Egamma,
                    const arma::mat &z,
                    const arma::mat &g,
					const arma::mat &i,
                    const arma::mat &j,
                    const int ND,
                    const int NI,
                    const int NL
                    ) {

		
	signed int n, l;
	double tmp;
	
	// Note: This updates both eta and eta^2.  Standard notation is that eta^2 should be etaeta, but keep shorter inconsistent naming for now

    // Main Calculation
#pragma omp parallel for private(n,l,tmp)
  	for(n=0; n < NI; n++){

		// i = i(l,0);
		// g = g(n,0);

		tmp = 0;

		for(l=0; l < NL; l++){

			if( i(l,0)==n ) tmp += Eystar(l,0)*Eb(j(l,0),0) - Eba(j(l,0),0) - Ebb(j(l,0))*accu(Egamma.row(g(n,0)) % z.row(n)) ;

		}

		// Note this assumes Veta is A_inv, not A
		Eta(n,0) = tmp*Veta(n,0);
		Eta2(n,0) = Eta(n,0)*Eta(n,0) + Veta(n,0);

//		Rcout << "Eta: " << Eta(n,0) << " Veta: " << Veta(n,0) << " Eta2: " <<  Eta2(n,0) << std::endl;

  	}  // for(n=0; n < NI; n++){

    return; 
} 
