\name{dynIRT}
\alias{dynIRT}
\title{ Dynamic IRT estimation via Variational Inference }
\description{ 
  \code{ordIRT} estimates an dynamic IRT model with two response categories per item, over several
  sessions. Ideal points over time follow a random walk prior, and the model originates from the
  work of Martin and Quinn (2002). Estimation is conducted using the variational EM algorithm described
  in the reference paper below. The algorithm will produce point estimates that are comparable to those
  of \code{\link[MCMCpack]{MCMCdynamicIRT1d}}, but will do so much more rapidly and also scale better
  with larger data sets.
}

\usage{
  dynIRT(.data, .starts = NULL, .priors = NULL, .control = NULL) 
}

\arguments{
  \item{.data}{ a list with the following items.
      \itemize{
      \item{\code{rc}}{  A (N x J) matrix of observed votes. `1' and `-1' are the yea and nay codes, while `0' is
      	a missing data code.}
      \item{\code{startlegis}}{  An (N x 1) matrix indicating the first session that each legislator serves. Justices
        are assumed to serve in all terms between (and including) `startlegis' to `endlegis'.  Terms start at term 0,
        and end at term T - 1.}
      \item{\code{endlegis}}{  An (N x 1) matrix indicating the last session that each legislator serves. Justices
        are assumed to serve in all terms between (and including) `startlegis' to `endlegis'.  Terms start at term 0,
        and end at term T - 1.}
      \item{\code{bill.session}}{  A (J x 1) matrix of integers indicating the session each bill occurred in. Session count
        begins at 0, so the maximum value of bill.session is T - 1.}
      \item{\code{T}}{  integer, indicating total number of consecutive terms in the data. Count starts from 1, so the
        maximum values of startlegis/endlegis/bill.session is T - 1, since they start from term 0.}
      }
    }
  \item{.starts}{ a list containing several matrices of starting values for the parameters. The
     list should contain the following matrices:
    \itemize{
    \item{\code{alpha}}{  A (J x 1) matrix of starting values for the item difficulty parameter \eqn{\alpha_j}.}
    \item{\code{beta}}{  A (J x 1) matrix of starting values for the item discrimination parameter \eqn{\beta_j}.}
    \item{\code{x}}{  An (N x T) matrix of starting values for the respondent ideal points \eqn{x_{it}}, with rows indicating the
      legislator and columns indicating the session.  Although not strictly necessary, it is generally good practice here to
      set the start values for legislators who are not serving in a particular session to 0, as that is what the point estimate for
      them will return.}
    }
  }
  \item{.priors}{  list, containing several matrices of starting values for the parameters. The
     list should contain the following matrices:
    \itemize{
    \item{\code{x.mu0}}{  A (N x 1) prior means matrix for respondent ideal points \eqn{c_{i0}}. These are generally
      set to be somewhat informative to resolve the standard rotational invariance problem in ideal point models.}
    \item{\code{x.sigma0}}{  A (N x 1) prior variance matrix for respondent ideal points \eqn{C_{i0}}.}
    \item{\code{beta.mu}}{  A (2 x 1) prior means matrix for all bill parameters \eqn{\alpha_j} and \eqn{\beta_j}.}
    \item{\code{beta.sigma}}{  A (2 x 2) prior covariance matrix for all bill parameters \eqn{\alpha_j} and \eqn{\beta_j}.}
    \item{\code{omega2}}{  A (N x 1) matrix with the evolutionary variance for each legislator \eqn{\omega^2_{ix}}.}
     }
  }
  \item{.control}{  list, specifying some control functions for estimation.  Options include the following:
    \itemize{
    \item{\code{threads}}{  integer, indicating number of cores to use. Default is to use a single core, but more can be
    	supported if more speed is desired.}
    \item{\code{verbose}}{  boolean, indicating whether output during estimation should be verbose or not.  Set FALSE by default.}
    \item{\code{thresh}}{  numeric. Algorithm will run until all parameters correlate at 1 - thresh across
    	consecutive iterations. Set at 1e-6 by default.}
    \item{\code{maxit}}{  integer. Sets the maximum number of iterations the algorithm can run. Set at 500 by default.}
    \item{\code{checkfreq}}{  integer. Sets frequency of verbose output by number of iterations. Set at 50 by default.}
    }
  }
}

\value{
    An object of class \code{dynIRT}.

  \item{means}{ list, containing several matrices of point estimates for the parameters corresponding
    to the inputs for the priors.  The list should contain the following matrices.
    \itemize{
    \item{\code{x}}{  A (N x T) matrix of point estimates for the respondent ideal points \eqn{x_{it}}, with rows
      indicating respondent and columns indicating session of the estimated ideal point. Estimates will equal
      exactly 0 for legislator/period combinations in which the legislator did not serve.}
    \item{\code{alpha}}{  A (J x 1) matrix of point estimates for the item difficulty parameter \eqn{\alpha_j}.}
    \item{\code{beta}}{  A (J x 1) matrix of point estimates for the item discrimination parameter \eqn{\beta_j}.}
    }                 
   }

  \item{vars}{ list, containing several matrices of variance estimates for parameters corresponding
    to the inputs for the priors.  Note that these variances are those recovered via variational approximation,
    and in most cases they are known to be far too small and generally unusable.  Better estimates of variances
    can be obtained manually via the parametric bootstrap.  The list should contain the following matrices:
    \itemize{
    \item{\code{x}}{  A (N x T) matrix of variance estimates for the respondent ideal points \eqn{x_{it}}, with rows
      indicating respondent and columns indicating session of the estimated ideal point. Estimates will equal
      exactly 0 for legislator/period combinations in which the legislator did not serve.}
    \item{\code{alpha}}{  A (J x 1) matrix of variance estimates for the item difficulty parameter \eqn{\alpha_j}.}
    \item{\code{beta}}{  A (J x 1) matrix of variance estimates for the item discrimination parameter \eqn{\beta_j}.}
    }                 
   }

  \item{runtime}{ A list of fit results, with elements listed as follows:}
    \itemize{
    \item{\code{iters}}{  integer, number of iterations run.}
    \item{\code{conv}}{  integer, convergence flag. Will return 1 if threshold reached, and 0 if maximum number of iterations reached.}
    \item{\code{threads}}{  integer, number of threads used to estimated model.}
    \item{\code{tolerance}}{  numeric, tolerance threshold for convergence.  Identical to thresh argument in input to .control list.}
    }                 
  \item{N}{ Number of respondents in estimation, should correspond to number of rows in roll call matrix.}
  \item{J}{ Number of items in estimation, should correspond to number of columns in roll call matrix.}
  \item{T}{ Number of time periods fed into the estimation, identical to argument input from .data list.}

  \item{call}{ Function call used to generate output.}
}


\author{ 
Kosuke Imai \email{kimai@princeton.edu}

James Lo \email{lojames@usc.edu}

Jonathan Olmsted \email{jpolmsted@gmail.com}
}

\references{

Original model and the example is based off of Andrew Martin and Kevin Quinn, ``Dynamic Ideal
Point Estimation via Markov Chain Monte Carlo for the U.S. Supreme Court, 1953-1999.'' Political
Analysis 10(2) 134-153.  Also see \url{http://mqscores.berkeley.edu}.

Variational model is described in Kosuke Imai, James Lo, and Jonathan Olmsted ``Fast Estimation
of Ideal Points with Massive Data.'' Working Paper.
Available at \url{http://imai.princeton.edu/research/fastideal.html}.

}

\seealso{
  '\link{mq_data}'.
}

\examples{

### Replication of Martin-Quinn Judicial Ideology Scores
### Based on July 23, 2014 (2014 Release 01) release of the Supreme Court Database
### Start values and priors based on replication code provided by Kevin Quinn

data(mq_data)

## Estimate dynamic variational model using dynIRT()
lout <- dynIRT(.data = mq_data$data.mq,
                    .starts = mq_data$cur.mq,
                    .priors = mq_data$priors.mq,
                    .control = {list(
                    threads = 1,
                    verbose = TRUE,
                    thresh = 1e-6,
		    maxit=500
                    )})

## Extract estimate from variational model
## Delete point estimates of 0, which are justices missing from that session
vi.out <- c(t(lout$means$x))
vi.out[vi.out==0] <- NA
vi.out <- na.omit(vi.out)

## Compare correlation against MCMC-estimated result
## Correlates at r=0.93 overall, and 0.96 when excluding Douglas
cor(vi.out, mq_data$mq_mcmc)	
cor(vi.out[mq_data$justiceName != "Douglas"],
 mq_data$mq_mcmc[mq_data$justiceName != "Douglas"])	


}
\keyword{ multivariate }
