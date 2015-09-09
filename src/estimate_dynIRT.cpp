// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#define DEBUG false

#ifdef _OPENMP
#include <omp.h>
#endif

#include <RcppArmadillo.h>
#include "getEystar_dynIRT.h"
#include "getLBS_dynIRT.h"
#include "getNlegis_dynIRT.h"
#include "getEx2x2_dynIRT.h"
#include "getVb2_dynIRT.h"
#include "getEb2_dynIRT.h"
#include "getVb_dynIRT.h"
#include "getVa_dynIRT.h"
#include "getEba_dynIRT.h"
#include "getEbb_dynIRT.h"
#include "getLast_dynIRT.h"
#include "getX_dynIRT.h"
#include "getOnecol_dynIRT.h"
#include "checkConv_dynIRT.h"

using namespace Rcpp ;

List estimate_dynIRT(arma::mat alpha_start,
               arma::mat beta_start,
               arma::mat x_start,
               arma::mat y,
               arma::mat startlegis,
               arma::mat endlegis,
               arma::mat bill_session,
               unsigned int T,
               arma::mat xmu0,
               arma::mat xsigma0,
               arma::mat betamu,
               arma::mat betasigma,
               arma::mat omega2,
               unsigned int threads = 0,
               bool verbose = true,
               unsigned int maxit = 2500,
               double thresh = 1e-6,
               unsigned int checkfreq = 50
               ) {

    //// Data Parameters
    unsigned int nJ = y.n_cols ;
    unsigned int nN = y.n_rows ;
   
    //// Admin
    unsigned int threadsused = 0 ;
	int convtype=1;
    unsigned int counter = 0 ;
    int isconv = 0;

    //// Init "Current" Containers
    arma::mat curEystar(nN, nJ, arma::fill::zeros);
    arma::mat curEa = alpha_start;
    arma::mat curEb = beta_start;
	arma::mat Nlegis_session;	// T x 1, each element has N number of legislators for session t
	arma::mat legis_by_session;	// T rows, each row has vector of legislators in session
	arma::cube curEx2x2(2, 2, T, arma::fill::zeros);
	arma::cube curVb2(2, 2, T, arma::fill::zeros);
    arma::mat curEb2(nJ, 2) ;
    arma::mat curVb;
    arma::mat curVa;
	arma::mat end_session;		// Filled by getLast()
	arma::mat ones_col;

    arma::mat curEx = x_start;					// (nNxT) matrix
	// Filling zeroes important for Vx, as getEx2x2() assumes Vx=0 for missing legislators for that period
	// Since Vx never gets updated for missing legislators, only need to do this once
    arma::mat curVx(nN, T, arma::fill::zeros);
	unsigned int i, j;

	// Clean curEx, setting all values of x_{it} that are not estimated (i.e. before startlegis and after endlegis) to 0
	// This guarantees that for final output, ideal points not estimated are output as 0 regardless of starting value
	for(i=0; i < nN; i++){
		for(j=0; j<T; j++){
			if(j < startlegis(i,0)) curEx(i,j) = 0;
			if(j > endlegis(i,0)) curEx(i,j) = 0;
		}
	}

	arma::mat curEbb(nJ,1);
	arma::mat curEba(nJ,1);
	for(j=0; j<nJ; j++){
		curEbb(j,0) = beta_start(j,0)*beta_start(j,0);
		curEba(j,0) = alpha_start(j,0)*beta_start(j,0);
	}
	
    
    //// Init "Old" Containers
    arma::mat oldEa = alpha_start;
    arma::mat oldEb = beta_start;
    arma::mat oldEx = curEx;


    // OpenMP Support
#ifdef _OPENMP
    omp_set_num_threads(1) ;
    if (threads > 0) {
        omp_set_num_threads(threads) ;
        threadsused = omp_get_max_threads() ;
    }
#endif

        // It turns out legis_by_session isn't necessary unless missing value in Ex is not 0
        // But only computed once, so no point changing it now
		legis_by_session = getLBS_dynIRT(startlegis, endlegis, T, nN);
    	Nlegis_session = getNlegis_dynIRT(legis_by_session, T, nN);
		end_session = getLast_dynIRT(bill_session, T, nJ);
		ones_col = getOnecol_dynIRT(startlegis, endlegis, T, nN);

// Main Loop Until Convergence
		while (counter < maxit) {
		        counter++ ;


        getEystar_dynIRT(curEystar,curEa,curEb,curEx,y,bill_session,startlegis, endlegis, nN,nJ);

		getX_dynIRT(curEx, curVx,
				curEbb, omega2,
				curEb, curEystar, curEba,
				startlegis, endlegis,
				xmu0, xsigma0,
				T, nN, end_session);

		getEx2x2_dynIRT(curEx2x2,curEx,curVx,legis_by_session,Nlegis_session,T);
    	getVb2_dynIRT(curVb2, curEx2x2, betasigma, T);
		getEb2_dynIRT(curEb2, curEystar,curEx,curVb2,bill_session,betamu,betasigma,ones_col,nJ);

		curEa = curEb2.col(0);
		curEb = curEb2.col(1);

		curEba = getEba_dynIRT(curEa,curEb,curVb2,bill_session,nJ);
		curEbb = getEbb_dynIRT(curEb,curVb2,bill_session,nJ);

        // Check for Interupt & Update Progress
        if (counter % checkfreq == 0) {
            R_CheckUserInterrupt() ;
            if (verbose) {
                Rcout << "Iteration: " << counter << std::endl ;
            }
        }

		// Counter>2 allows starts of curEx at 0
		if(counter > 2) isconv = checkConv_dynIRT(oldEx, curEx, oldEb, curEb, oldEa, curEa, thresh, convtype);
        if (isconv==1) break;

        // Update Old Values If Not Converged
        oldEx = curEx;
        oldEb = curEb;
        oldEa = curEa;

		}
// LOOP ENDING HERE

		// Only needed after convergence
    	curVb = getVb_dynIRT(curVb2, bill_session, nJ);
    	curVa = getVa_dynIRT(curVb2, bill_session, nJ);

// 	Rcout << "\n Completed after " << counter << " iterations..." << std::endl ;

    List ret ;
    List means ;
    List vars ;
    List runtime ;

    means["x"] = curEx;
    means["alpha"] = curEa;
    means["beta"] = curEb;
	
    vars["x"] = curVx ;
    vars["alpha"] = curVa ;
    vars["beta"] = curVb;
    
    runtime["iters"] = counter ;
    runtime["conv"] = isconv ;
    runtime["threads"] = threadsused ;
    runtime["tolerance"] = thresh ;

    ret["means"] = means ;
    ret["vars"] = vars ;
    ret["runtime"] = runtime ;

    ret["N"] = nN ;
    ret["J"] = nJ ;
    ret["T"] = T ;

    return(ret) ;
}
