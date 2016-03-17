
#ifdef _OPENMP
#include <omp.h>
#endif

#include "getEalpha_poisIRT.h"
#include "getEx_poisIRT.h"
#include "getVbeta_poisIRT.h"
#include "getVx_poisIRT.h"
#include "getExi_poisIRT.h"
#include "getEbeta_poisIRT.h"
#include "getEpsi_poisIRT.h"
#include "checkConv_poisIRT.h"

#include <RcppArmadillo.h>

using namespace Rcpp ;

List estimate_poisIRT (arma::mat alpha_start,
               arma::mat psi_start,
               arma::mat beta_start,
               arma::mat x_start,
               arma::mat Y,
               arma::mat i,
               int ni,
	           double psi_mu,
               double psi_sigma,
               double alpha_mu,
               double alpha_sigma,
               double beta_mu,
               double beta_sigma,
               double x_mu,
               double x_sigma,
               unsigned int threads = 0,
               bool verbose = true,
               unsigned int maxit = 2500,
               double thresh = 1e-9,
               unsigned int checkfreq = 50
               ) {

    //// Admin
    unsigned int threadsused = 0 ;
	int convtype=1;
    unsigned int counter = 0 ;
    int isconv = 0;
    
    int k;
    
    arma::mat curEalpha = alpha_start;
    arma::mat curEpsi = psi_start;
    arma::mat curEbeta = beta_start;
    arma::mat curEx = x_start;

    arma::mat oldEalpha = alpha_start;
    arma::mat oldEpsi = psi_start;
    arma::mat oldEbeta = beta_start;
    arma::mat oldEx = x_start;
    
    int NK = alpha_start.n_rows;
    int NJ = psi_start.n_rows;
	int NI = ni;

	arma::mat curExfull(NK, 1);

	arma::mat curVbeta(NJ, 1, arma::fill::ones);
 	arma::mat curVx(NI, 1, arma::fill::ones);

	arma::mat onesNK(NK, 1, arma::fill::ones);
	arma::mat onesNJ(NJ, 1, arma::fill::ones);

	// These are the variational parameters xi
    arma::mat exi(NJ, NK, arma::fill::ones);
    arma::mat xi(NJ, NK, arma::fill::zeros);
    arma::mat exixi(NJ, NK, arma::fill::zeros);

    // OpenMP Support
#ifdef _OPENMP
    omp_set_num_threads(1) ;
    if (threads > 0) {
        omp_set_num_threads(threads) ;
        threadsused = omp_get_max_threads() ;
    }
#endif


	while(counter < maxit){
		counter++;

		getExi(exi, xi, exixi, curEalpha, curEpsi, curEbeta, curEx, i, NI, NK, NJ);
		getEalpha(curEalpha, curEx, onesNJ, exi, xi, exixi, Y, i, curEpsi, curEbeta, alpha_mu, alpha_sigma, NI, NK, NJ);

		getVbeta(curVbeta, curEx, curVx, exi, i, beta_sigma, NI, NK, NJ);
		getVx(curVx, curEbeta, curVbeta, exi, i, x_sigma, NI, NK, NJ);
		getEx(curEx, curVx, curEalpha, onesNJ, exi, xi, exixi, Y, curEpsi, curEbeta, curVbeta, i, x_mu, x_sigma, NI, NK, NJ);

		for(k=0; k<NK; k++)	curExfull(k,0) = curEx(i(k,0),0);
		getEbeta(curEbeta, curVbeta, curEalpha, onesNK, exi, xi, exixi, Y, curEpsi, curExfull, curVx, i, beta_mu, beta_sigma, NI, NK, NJ);

		getEpsi(curEpsi, curExfull, onesNK, exi, xi, exixi, Y, i, curEalpha, curEbeta, psi_mu, psi_sigma, NI, NK, NJ);

        // Check for Interupt & Update Progress
        if (counter % checkfreq == 0) {
            R_CheckUserInterrupt() ;
            if (verbose) {
                Rcout << "Iteration: " << counter << std::endl ;
            }
        }

		// Counter>2 allows starts of curEx at 0
		if(counter > 2) isconv = checkConv_poisIRT(oldEalpha, curEalpha, oldEpsi, curEpsi, oldEbeta, curEbeta, oldEx, curEx, thresh, convtype);
        if (isconv==1) break;

		oldEalpha = curEalpha;
		oldEpsi = curEpsi;
		oldEbeta = curEbeta;
		oldEx = curEx;

	}

// 	Rcout << "\n Completed after " << counter << " iterations..." << std::endl ;

    List result;
    List means;
    List vars;
    List runtime;
    List N;
	List i_of_k;

    means["alpha"] = curEalpha;
    means["beta"] = curEbeta;
    means["psi"] = curEpsi;
    means["x"] = curEx;
	means["xi"] = xi;

    vars["x"] = curVx ;
    vars["beta"] = curVbeta;
	
    runtime["iters"] = counter ;
    runtime["conv"] = isconv ;
    runtime["threads"] = threadsused ;
    runtime["tolerance"] = thresh ;

    N["K"] = NK ;
    N["J"] = NJ ;
    N["I"] = ni ;

    result["means"] = means;
    result["vars"] = vars;
    result["runtime"] = runtime;
    result["N"] = N;
    result["i_of_k"] = i;

    return(result) ;
}
