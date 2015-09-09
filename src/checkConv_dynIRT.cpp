// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

int checkConv_dynIRT(const arma::mat &oldEx,
              const arma::mat &curEx,
			  const arma::mat &oldEb,
              const arma::mat &curEb,
			  const arma::mat &oldEa,
              const arma::mat &curEa,
              double thresh,
              int convtype
               ){

    double devEx, devEa, devEb;
    devEx=100.0;
    devEa=100.0;
    devEb=100.0;

	arma::vec oldEx_vec = vectorise(oldEx);
	arma::vec curEx_vec = vectorise(curEx);
	arma::vec oldEx_stripped = oldEx_vec(arma::find(oldEx_vec != 0));
	arma::vec curEx_stripped = curEx_vec(arma::find(curEx_vec != 0));

  if (convtype == 1) {
        // correlation
        devEx = 1 - (cor(oldEx_stripped, curEx_stripped)).min() ;
        devEb = 1 - (cor(oldEb, curEb)).min() ;
        devEa = 1 - (cor(oldEa, curEa)).min() ;
    }

    if (convtype == 2) {
        // maximum absolute deviation
        devEx = (abs(curEx_stripped - oldEx_stripped)).max() ;
        devEb = (abs(curEb - oldEb)).max() ;
        devEa = (abs(curEa - oldEa)).max() ;
    }
    
   if( (devEx < thresh) & (devEb < thresh) & (devEa < thresh)) return(1);

   return(0) ;

}
