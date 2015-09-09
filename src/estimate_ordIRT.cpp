// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#define DEBUG false

#ifdef _OPENMP
#include <omp.h>
#endif

#include <RcppArmadillo.h>
#include "getEzstar_ordIRT.h"
#include "getEbb_ordIRT.h"
#include "getVx_ordIRT.h"
#include "getEx_ordIRT.h"
#include "getVzstar_ordIRT.h"
#include "getEzzstar_ordIRT.h"
#include "getEtt_ordIRT.h"
#include "getExx_ordIRT.h"
#include "getEdd_ordIRT.h"
#include "getEd_ordIRT.h"
#include "getEx2x2_ordIRT.h"
#include "getVb2_ordIRT.h"
#include "getEb2_ordIRT.h"
#include "checkConv_ordIRT.h"

using namespace Rcpp ;

List estimate_ordIRT(arma::mat tau_start,
               arma::mat DD_start,
               arma::mat beta_start,
               arma::mat x_start,
               arma::mat y,
               arma::mat xmu,
               arma::mat xsigma,
               arma::mat betamu,
               arma::mat betasigma,
               unsigned int D = 1,
               unsigned int threads = 0,
               bool verbose = true,
               unsigned int maxit = 300,
               double thresh = 1e-6,
               unsigned int checkfreq = 50
               ) {

    //// Data Parameters
    unsigned int j;
    unsigned int nJ = y.n_cols ;
    unsigned int nN = y.n_rows ;

    //// Admin
    unsigned int threadsused = 0 ;
    unsigned int counter = 0 ;
    int isconv = 0 ;
	int convtype = 1;  // For correlation-based convergence criterion
    arma::mat convtrace(maxit - 1,
                        (2 * D) + 1,
                        arma::fill::zeros
                        ) ;
    arma::mat statstrace(maxit - 1,
                         2,
                         arma::fill::zeros
                         ) ;
    arma::mat probs1 ;
    arma::mat probsobs ;

//    checkInputs();

    //// Init "Current" Containers
    arma::mat curEzstar;
    arma::mat curVzstar;
    arma::mat curEzzstar;

    arma::mat curEb2(nJ, 2);
	arma::cube curVb2(2, 2, nJ, arma::fill::zeros);
    arma::mat curEx2x2(2, 2);
        
    arma::mat curVx(nN, 1);
    arma::mat curVb(nJ, 1, arma::fill::ones);
    arma::mat curVtau(nJ, 1);
    arma::mat curExx(nN, 1);
    arma::mat curEtt(nJ, 1);
	curVx.fill(xsigma(0,0));
	curVb.fill(betasigma(0,0));
	
    arma::mat curEbb(nJ, 1);
    arma::mat curEdd = DD_start;
    arma::mat finalEd;
    
    arma::mat curEx = x_start;
    arma::mat curEb = beta_start;
    arma::mat curEtau = tau_start ;    

    //// Init "Old" Containers
    arma::mat oldEb = beta_start ;
    arma::mat oldEx = x_start ;
    arma::mat oldEtau = tau_start;
    arma::mat oldEdd = DD_start;

    // OpenMP Support
#ifdef _OPENMP
    omp_set_num_threads(1) ;
    if (threads > 0) {
        omp_set_num_threads(threads) ;
        threadsused = omp_get_max_threads() ;
    }
#endif



    // Main Loop Until Convergence, commented out
		while (counter < maxit) {
		        counter++ ;

          //Rcout << counter << std::endl ;

	      curEzstar = getEzstar_ordIRT(curEdd, curEb, curEx, curEtau, y, nN, nJ);
	      
		  curEbb = getEbb_ordIRT(curEb, curVb, nJ);
		  curVx = getVx_ordIRT(curEbb, curEdd, xsigma);
		  curEx = getEx_ordIRT(curEzstar, curEb, curEtau, curVx, curEdd, xmu, xsigma, curVb2, nN, nJ);
		  curExx = getExx_ordIRT(curEx, curVx, nN); 

		  curEx2x2 = getEx2x2_ordIRT(curEx, curVx, nN);
		  curVb2 = getVb2_ordIRT(curEx, curEx2x2, betasigma, curEdd, nJ);
		  curEb2 = getEb2_ordIRT(curEzstar, curEx, curVb2, betamu, betasigma, curEdd, nJ);
		  curEtau = curEb2.col(0);
		  curEb = curEb2.col(1);
    	  for(j=0; j<nJ; j++){
			  curVtau(j,0) = curVb2.slice(j)(0,0);
			  curVb(j,0) = curVb2.slice(j)(1,1);
		  }

	      curVzstar = getVzstar_ordIRT(curEdd, curEb, curEx, curEtau, y, nN, nJ);
	      curEzzstar = getEzzstar_ordIRT(curEzstar, curVzstar, nN, nJ);

		  curEtt = getEtt_ordIRT(curEtau, curVtau, nJ);
		  curEdd = getEdd_ordIRT(curEx,curExx,curEb,curEbb,curEtau,curEtt,curEzstar,curEzzstar,nN,nJ);

        // Check for Interupt & Update Progress
        if (counter % checkfreq == 0) {
            R_CheckUserInterrupt() ;
            if (verbose) {
                Rcout << "Iteration: " << counter << std::endl ;
            }
        }

        // Check Convergence on Its >= 2
		if(counter > 2) isconv = checkConv_ordIRT(oldEx, curEx, oldEb, curEb, oldEtau, curEtau, oldEdd, curEdd, thresh, convtype);
        if (isconv) break;

        // Update Old Values If Not Converged
        oldEx = curEx ;
        oldEb = curEb ;
		oldEtau = curEtau;
		oldEdd = curEdd;

		}
// LOOP ENDING HERE

		// E[Delta] not needed for updates, only for auxiliary testing
		// Hence we only calculate it once at the end once converged
		finalEd = getEd_ordIRT(curEx,curExx,curEb,curEbb,curEtau,curEtt,curEzstar,curEzzstar,nN,nJ);

    List ret ;
    List means ;
    List vars ;
    List runtime ;

    means["x"] = curEx ;
    means["beta"] = curEb ;
    means["tau"] = curEtau ;
    means["Delta_sq"] = curEdd ;
    means["Delta"] = finalEd ;
	
    vars["x"] = curVx ;
    vars["beta"] = curVb;
    vars["tau"] = curVtau;
    
    runtime["iters"] = counter ;
    runtime["conv"] = isconv ;
    runtime["threads"] = threadsused ;
    runtime["tolerance"] = thresh ;

    ret["means"] = means ;
    ret["vars"] = vars ;
    ret["runtime"] = runtime ;

    ret["N"] = nN ;
    ret["J"] = nJ ;

    return(ret) ;
}
