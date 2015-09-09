// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEb2_dynIRT(arma::mat &Eb2,
				 const arma::mat &Eystar,
                 const arma::mat &Ex,
                 const arma::cube &Vb2,
                 const arma::mat &bill_session,
                 const arma::mat &mubeta,
                 const arma::mat &sigmabeta,
                 const arma::mat &ones_col,
                 const int J
                 ) {

// In written implementation below, Ea=sum ystar_i and Eb = sum x_i*ystar_i

	int t,j;
	arma::mat Ex2;

#pragma omp parallel for private(j,t,Ex2)	
    for (j = 0; j < J; j++) {

		t = bill_session(j,0);
	    Ex2 = Ex.col(t);

		// We cannot just use ones here, as we have to zero out legislators not present
	    Ex2.insert_cols(0, ones_col.col(t));    
        Eb2.row(j) = trans(Vb2.slice(t) * (inv_sympd(sigmabeta) * mubeta + trans(Ex2) * Eystar.col(j))) ;

    }

    return;
    
}
