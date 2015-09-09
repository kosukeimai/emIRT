// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEtt_ordIRT(const arma::mat &Etau,
                const arma::mat &Vtau,
                const int nJ
                ){


	int j;
	arma::mat result(nJ, 1);

#pragma omp parallel for
    for(j=0; j < nJ; j++){

		// For variational update
		// result(j,0) = Etau(j,0)*Etau(j,0) + Vtau(j,0);

		// For ECM update
		result(j,0) = Etau(j,0)*Etau(j,0);

    }

    return(result) ;
}
