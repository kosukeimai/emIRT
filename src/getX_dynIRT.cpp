// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getX_dynIRT(arma::mat &Ex,
                arma::mat &Vx,
                const arma::mat &Ebb,
                const arma::mat &omega2,
                const arma::mat &Eb,
                const arma::mat &Eystar,
                const arma::mat &Eba,
                const arma::mat &startlegis,
                const arma::mat &endlegis,
                const arma::mat &xmu0,
                const arma::mat &xsigma0,
                const int T,
                const int N,
                const arma::mat &end_session
                ) {


	int i, t;

	arma::mat betaDD(T,1,arma::fill::zeros);
	arma::mat Eba_sum(T,1,arma::fill::zeros);
	arma::mat Eby_sum, Eb_t, Eystar_t;

	arma::mat Ot(N,T,arma::fill::zeros);
	arma::mat Kt(N,T,arma::fill::zeros);
	arma::mat St(N,T,arma::fill::zeros);
	arma::mat Jt(N,T,arma::fill::zeros);
	arma::mat C_var(N,T,arma::fill::zeros);
	arma::mat c_mean(N,T,arma::fill::zeros);
	double yDD;

	// These quantities are called repeatedly, calculate and store for reuse
	//betaDD and yDD corresponde to beta.dot.dot and y.dot.dot respectively
	betaDD(0,0) = sqrt(accu(Ebb.submat(0,0,end_session(0,0)-1,0)));
	Eba_sum(0,0) = accu(Eba.submat(0,0,end_session(0,0)-1,0));
	if(T > 1){

	#pragma omp parallel for
		for(t = 1; t < T; t++){
			betaDD(t,0) = sqrt(accu(Ebb.submat(end_session(t-1,0),0,end_session(t,0)-1,0)));
			Eba_sum(t,0) = accu( Eba.submat(end_session(t-1,0),0,end_session(t,0)-1,0) );
		}
	}

#pragma omp parallel for private(t,Eystar_t,Eb_t,Eby_sum,yDD)
	for(i=0; i < N; i++){
		
		// Initialize first period forward filter using priors
		t = startlegis(i,0);

		if(t==0){
			Eystar_t = Eystar.submat(i,0,i,end_session(t,0)-1);
			Eb_t = Eb.submat(0,0,end_session(t,0)-1,0);
		}
		if(t != 0){
			Eystar_t = Eystar.submat(i,end_session(t-1,0),i,end_session(t,0)-1);
			Eb_t = Eb.submat(end_session(t-1,0),0,end_session(t,0)-1,0);
		}

		Eby_sum = Eystar_t * Eb_t;
		yDD = (Eby_sum(0,0) - Eba_sum(t,0))/betaDD(t,0);
		Ot(i,t) = omega2(i,0) + xsigma0(i,0);
		St(i,t) = betaDD(t,0)*betaDD(t,0)*Ot(i,t) + 1;
		Kt(i,t) = betaDD(t,0)*Ot(i,t)/St(i,t);
		C_var(i,t) = (1 - Kt(i,t)*betaDD(t,0))*Ot(i,t);
		c_mean(i,t) = xmu0(i,0) + Kt(i,t)*(yDD - betaDD(t,0)*xmu0(i,0));

		//Forward-filter test only
		//Vx(i,t) = C_var(i,t);
		//Ex(i,t) = c_mean(i,t);
			
		// If justice in only one period (never tested)
		if(startlegis(i,0) == endlegis(i,0)){
			Vx(i,t) = C_var(i,t);
			Ex(i,t) = c_mean(i,t);
		}

		// If justice in multiple periods (should be most instances)
		if(startlegis(i,0) != endlegis(i,0)){

			// Finish forward filtering
			for(t = startlegis(i,0) + 1; t <= endlegis(i,0); t++){

				Eystar_t = Eystar.submat(i,end_session(t-1,0),i,end_session(t,0)-1);
				Eb_t = Eb.submat(end_session(t-1,0), 0, end_session(t,0)-1, 0);

				Eby_sum = Eystar_t * Eb_t;
				yDD = (Eby_sum(0,0) - Eba_sum(t,0))/betaDD(t,0);
				Ot(i,t) = omega2(i,0) + C_var(i,t-1);
				St(i,t) = betaDD(t,0)*betaDD(t,0)*Ot(i,t) + 1;
				Kt(i,t) = betaDD(t,0)*Ot(i,t)/St(i,t);
				C_var(i,t) = (1 - Kt(i,t)*betaDD(t,0))*Ot(i,t);
				c_mean(i,t) = c_mean(i,t-1) + Kt(i,t)*(yDD - betaDD(t,0)*c_mean(i,t-1));

				//Forward-filter test only
				//Vx(i,t) = C_var(i,t);
				//Ex(i,t) = c_mean(i,t);

			}

			// Initialize backward sampling here
			Vx(i, endlegis(i,0)) = C_var(i,endlegis(i,0));
			Ex(i, endlegis(i,0)) = c_mean(i,endlegis(i,0));

			for(t = endlegis(i,0) - 1; t >= startlegis(i,0); t--){
				Jt(i,t) = C_var(i,t)/Ot(i,t+1);
				Vx(i,t) = C_var(i,t) + Jt(i,t)*Jt(i,t)*(Vx(i,t+1) - Ot(i,t+1));
				Ex(i,t) = c_mean(i,t) + Jt(i,t)*(Ex(i,t+1) - c_mean(i,t));
			}

		} //end if(startlegis(i,0) != endlegis(i,0));

		
	} 	// for(i=0; i < N; i++)

	return;

}
