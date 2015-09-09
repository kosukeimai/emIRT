// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

int checkConv_ordIRT(const arma::mat &oldEx,
              const arma::mat &curEx,
			  const arma::mat &oldEb,
              const arma::mat &curEb,
			  const arma::mat &oldEtau,
              const arma::mat &curEtau,
			  const arma::mat &oldEdd,
              const arma::mat &curEdd,
              double thresh,
              int convtype
               ){

    double devEx, devEtau, devEb, devEdd;
    devEx=100.0;
    devEtau=100.0;
    devEb=100.0;
    devEdd=100.0;

    if (convtype == 1) {
        // correlation
        devEx = 1 - (cor(oldEx, curEx)).min() ;
        devEb = 1 - (cor(oldEb, curEb)).min() ;
        devEtau = 1 - (cor(oldEtau, curEtau)).min() ;
        devEdd = 1 - (cor(oldEdd, curEdd)).min() ;
    }

    if (convtype == 2) {
        // maximum absolute deviation
        devEx = (abs(curEx - oldEx)).max() ;
        devEb = (abs(curEb - oldEb)).max() ;
        devEtau = (abs(curEtau - oldEtau)).max() ;
        devEdd = (abs(curEdd - oldEdd)).max() ;
    }
    
    
   //Rcout << devEx << " " << devEtau << " " << devEb << " " << devEdd << std::endl ;
    
   if( (devEx < thresh) & (devEb < thresh) & (devEtau < thresh) & (devEdd < thresh)) return(1);
  
   return(0) ;

}
