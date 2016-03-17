// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-
	
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getExi(arma::mat &exi,
					arma::mat &xi,
					arma::mat &exixi,
				   const arma::mat &Ealpha,
				   const arma::mat &Epsi,
                   const arma::mat &Ebeta,
                   const arma::mat &Ex,
                   const arma::mat &i,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) {

	int j, k;

#pragma omp parallel for private(k)
	for(j=0; j<NJ; j++){
		for(k=0; k<NK; k++){
			xi(j,k) = Ealpha(k,0) + Epsi(j,0) + Ebeta(j,0)*Ex(i(k,0),0);

			exi(j,k) = exp(xi(j,k));
			exixi(j,k) = exi(j,k)*xi(j,k);
		}
	}

    return;

}



