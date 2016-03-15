// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

int checkConv_poisIRT(const arma::mat &oldEalpha,
              const arma::mat &curEalpha,
			  const arma::mat &oldEpsi,
              const arma::mat &curEpsi,
			  const arma::mat &oldEbeta,
              const arma::mat &curEbeta,
			  const arma::mat &oldEx,
              const arma::mat &curEx,
              double thresh,
              int convtype
               ){

    double devEalpha, devEpsi, devEbeta, devEx;
    devEalpha = 100.0;
    devEpsi = 100.0;
    devEbeta = 100.0;
    devEx = 100.0;

  if (convtype == 1) {
        // correlation
        devEalpha = 1 - (cor(oldEalpha, curEalpha)).min() ;
        devEpsi = 1 - (cor(oldEpsi, curEpsi)).min() ;
        devEbeta = 1 - (cor(oldEbeta, curEbeta)).min() ;
        devEx = 1 - (cor(oldEx, curEx)).min() ;

    }

    if (convtype == 2) {
        // maximum absolute deviation
        devEalpha = (abs(curEalpha - oldEalpha)).max() ;
        devEpsi = (abs(curEpsi - oldEpsi)).max() ;
        devEbeta = (abs(curEbeta - oldEbeta)).max() ;
        devEx = (abs(curEx - oldEx)).max() ;
    }
    
   if( (devEalpha < thresh) & (devEpsi < thresh) & (devEbeta < thresh) & (devEx < thresh)) return(1);

   return(0) ;

}
