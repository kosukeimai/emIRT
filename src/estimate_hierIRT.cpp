// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#ifdef _OPENMP
#include <omp.h>
#endif

#include <RcppArmadillo.h>

#include "getEystar_hierIRT.h"
#include "getEx2_hierIRT.h"
#include "getEx2x2_hierIRT.h"
#include "getVeta_hierIRT.h"
#include "getEta_hierIRT.h"
#include "getVb2_hierIRT.h"
#include "getEb2_hierIRT.h"
#include "getEba_hierIRT.h"
#include "getEbb_hierIRT.h"
#include "getVgamma_hierIRT.h"
#include "getEgamma_hierIRT.h"
#include "getEgg_hierIRT.h"
#include "getEsigma_hierIRT.h"
#include "checkConv_hierIRT.h"

using namespace Rcpp;


List estimate_hierIRT(arma::mat alpha_start,
               arma::mat beta_start,
               arma::mat gamma_start,
               arma::mat sigma_start,
               arma::mat eta_start,
               arma::mat y,
               arma::mat z,
               arma::mat g,
               arma::mat i,
               arma::mat j,
               arma::mat gammamu,
               arma::mat gammasigma,
               arma::mat betamu,
               arma::mat betasigma,
               arma::mat sigmav,
               arma::mat sigmas,
               unsigned int ND,
               unsigned int NG,
               unsigned int NI,
               unsigned int NJ,
               unsigned int NL,
               unsigned int threads = 0,
               bool verbose = true,
               unsigned int maxit = 2500,
               double thresh = 1e-9,
               unsigned int checkfreq = 50
               ) {
    // Init Qtys

	unsigned int jloop;
    //// Admin
    unsigned int threadsused = 0 ;
	int convtype=1;
    unsigned int counter = 0 ;
    int isconv = 0;

    arma::mat curEa = alpha_start;
    arma::mat curEb = beta_start;
    arma::mat curEgamma = gamma_start;
    arma::mat curEsigma = sigma_start;
    arma::mat curEta = eta_start;

    arma::mat oldEa = alpha_start;
    arma::mat oldEb = beta_start;
    arma::mat oldEgamma = gamma_start;

	arma::mat curEystar(NL,1);
	arma::mat curEta2(NI,1);
	for(jloop=0; jloop<NI; jloop++)	curEta2(jloop,0) = curEta(jloop,0)*curEta(jloop,0);

	arma::mat curEbb(NJ,1);
	arma::mat curEba(NJ,1);
	for(jloop=0; jloop<NJ; jloop++){
		curEbb(jloop,0) = beta_start(jloop,0)*beta_start(jloop,0);
		curEba(jloop,0) = alpha_start(jloop,0)*beta_start(jloop,0);
	}

	arma::cube curEgg(ND, ND, NG);
	for(jloop=0; jloop < NG; jloop++)  curEgg.slice(jloop) = trans(curEgamma.row(jloop)) * curEgamma.row(jloop);

	arma::mat curEx2(NI,2);
	for(jloop=0; jloop < NI; jloop++)  curEx2(jloop,0) = 1;
	getEx2_hierIRT(curEx2, curEgamma, curEta, g, z, NI);
	arma::cube curEx2x2(2,2,NI);
	curEx2x2.fill(1.0);		// For top-left diagonal element
	getEx2x2_hierIRT(curEx2x2, curEx2, curEgg, curEgamma, curEta, curEta2, g, z, NI);

	arma::mat curEx;
	arma::mat oldEx;
		
	arma::mat curVeta(NI,1);
	arma::cube curVb2(2, 2, NJ);
	arma::mat curEb2(NJ, 2);
	arma::cube curVgamma(ND, ND, NG);
			
    // OpenMP Support
#ifdef _OPENMP
    omp_set_num_threads(1) ;
    if (threads > 0) {
        omp_set_num_threads(threads) ;
        threadsused = omp_get_max_threads() ;
    }
#endif

	// Note: Veta, Vb, and Vg are actually A_inv, B_inv, and C_inv, not A, B, and C

    // Main Loop Until Convergence
		while (counter < maxit) {
		        counter++;

//      Rcout << "\n ######### ITERATION: " << counter << std::endl ;
                
		getEystar_hierIRT(curEystar, y, z, g, i, j, curEa, curEb, curEgamma, curEta, ND, NG, NI, NJ, NL);

		curVeta = getVeta_hierIRT(i, j, g, curEsigma, curEbb, NI, NL);
		getEta_hierIRT(curEta, curEta2, curVeta, curEystar, curEb, curEba, curEbb, curEgamma, z, g, i, j, ND, NI, NL);

		getEx2_hierIRT(curEx2, curEgamma, curEta, g, z, NI);
		getEx2x2_hierIRT(curEx2x2, curEx2, curEgg, curEgamma, curEta, curEta2, g, z, NI);
		curEx = curEx2.col(1);

		getVb2_hierIRT(curVb2, betasigma, curEx2x2, i, j, NL, NJ);
		getEb2_hierIRT(curEb2, curVb2, betasigma, betamu, curEystar, curEx2, i, j, NL, NJ);
		curEa = curEb2.col(0);
		curEb = curEb2.col(1);
		getEbb_hierIRT(curEbb, curEb, curVb2, NJ);
		getEba_hierIRT(curEba, curEb, curEa, curVb2, NJ);

		getVgamma_hierIRT(curVgamma, gammasigma, curEbb, g, i, j, z, NL, NG);
		getEgamma_hierIRT(curEgamma, curVgamma, gammasigma, gammamu, g, i, j, z, curEb, curEbb, curEystar, curEa, curEta, NL, NG);
		getEgg_hierIRT(curEgg, curEgamma, curVgamma, NG);

		getEsigma_hierIRT(curEsigma, curEta2, sigmav, sigmas, g, NG, NI);

        // Check for Interupt & Update Progress
        if (counter % checkfreq == 0) {
            R_CheckUserInterrupt() ;
            if (verbose) {
                Rcout << "Iteration: " << counter << std::endl ;
            }
        }

		// Counter>2 allows starts of curEx at 0
		if(counter > 2) isconv = checkConv_hierIRT(oldEa, curEa, oldEb, curEb, oldEgamma, curEgamma, oldEx, curEx, ND, thresh, convtype);
        if (isconv==1) break;

        // Update Old Values If Not Converged
		oldEx = curEx;
		oldEa = curEa;
		oldEb = curEb;
		oldEgamma = curEgamma;

		}
// LOOP ENDING HERE

		// One more update of ideal points based on last set of converged values
		getEx2_hierIRT(curEx2, curEgamma, curEta, g, z, NI);
		curEx = curEx2.col(1);

    List result ;
    List N;
    List means;
    List vars;
    List runtime;

	N["D"] = ND;
	N["G"] = NG;
	N["I"] = NI;
	N["J"] = NJ;
	N["L"] = NL;

    runtime["iters"] = counter ;
    runtime["conv"] = isconv ;
    runtime["threads"] = threadsused ;
    runtime["tolerance"] = thresh ;

	means["sigma"] = curEsigma;
	means["eta"] = curEta;
	means["x_implied"] = curEx;
	means["gamma"] = curEgamma;
	means["alpha"] = curEa;
	means["beta"] = curEb;

	vars["eta"] = curVeta;
	vars["gamma"] = curVgamma;
	vars["beta2"] = curVb2;

	result["means"] = means;
	result["vars"] = vars;
	result["runtime"] = runtime;
	result["N"] = N;

    return(result) ;
}
