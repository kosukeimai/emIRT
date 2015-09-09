// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

int checkConv_hierIRT(const arma::mat &oldEa,
              const arma::mat &curEa,
			  const arma::mat &oldEb,
              const arma::mat &curEb,
			  const arma::mat &oldEgamma,
              const arma::mat &curEgamma,
			  const arma::mat &oldEx,
              const arma::mat &curEx,
			  int ND,
              double thresh,
              int convtype
               ){


    double devEa, devEb, devEgamma, tmp;
    devEa=100.0;
    devEb=100.0;
    devEgamma=100.0;

	int i;


    if (convtype == 1) {
        // correlation
        devEa = 1 - (cor(oldEa, curEa)).min() ;
        devEb = 1 - (cor(oldEb, curEb)).min() ;
        devEgamma = 1 - (cor(oldEgamma.col(0), curEgamma.col(0))).min() ;
		for(i=1; i<ND; i++){
        	tmp = 1 - (cor(oldEgamma.col(i), curEgamma.col(i))).min() ;
			if(tmp > devEgamma) devEgamma = tmp;
		}
    }

    if (convtype == 2) {
        // maximum absolute deviation
        devEa = (abs(curEa - oldEa)).max() ;
        devEb = (abs(curEb - oldEb)).max() ;
        devEgamma = (abs(curEgamma - oldEgamma)).max() ;
    }
    

   if( (devEa < thresh) & (devEb < thresh) & (devEgamma < thresh)) return(1);

  
   return(0) ;

}
