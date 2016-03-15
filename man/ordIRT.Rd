\name{ordIRT}
\alias{ordIRT}
\title{ Two-parameter Ordinal IRT estimation via EM }
\description{ 
  \code{ordIRT} estimates an ordinal IRT model with three ordered response categories. Estimation
  is conducted using the EM algorithm described in the reference paper below. The algorithm will
  produce point estimates that are comparable to those of \code{\link[MCMCpack]{MCMCordfactanal}},
  but will do so much more rapidly and also scale better with larger data sets.
}

\usage{
  ordIRT(.rc, .starts = NULL, .priors = NULL, .D = 1L, .control = NULL) 
}

\arguments{
  \item{.rc}{  matrix of numeric values containing the data to be scaled.  Respondents are
    assumed to be on rows, and items assumed to be on columns, so the matrix is assumed to be
    of dimension (N x J).  For each item, only 3 ordered category responses are accepted, and the
    only allowable responses are `1', `2', and `3', with `0' as a missing data record. If data
    of more than 3 categories are to be rescaled, they should be collapsed into 3 categories and
    recoded accordingly before proceeding.}
  \item{.starts}{ a list containing several matrices of starting values for the parameters. Note that
  the parameters here correspond to the re-parameterized version of the model (i.e. alpha is \eqn{alpha^*},
  not the original \eqn{alpha_{1j}}  The list should contain the following matrices:
    \itemize{
    \item{\code{beta}}{  A (J x 1) matrix of starting values for the reparameterized item discrimination parameter \eqn{\beta^*}.}
    \item{\code{x}}{  An (N x 1) matrix of starting values for the respondent ideal points \eqn{x_t}.}
    \item{\code{tau}}{  A (J x 1) matrix of starting values for the bill cutpoint \eqn{\tau_j}.}
    \item{\code{DD}}{  A (J x 1) matrix of starting values for the squared bill cutpoint difference \eqn{\tau_j^2}.}
    }
  }
  \item{.priors}{  list, containing several matrices of starting values for the parameters. Note that
  the parameters here correspond to the re-parameterized version of the model (i.e. alpha is \eqn{\alpha^*},
  not the original \eqn{\alpha_{1j}}  The list should contain the following matrices:
    \itemize{
    \item{\code{x$mu}}{  A (1 x 1) prior means matrix for respondent ideal points \eqn{x_t}.}
    \item{\code{x$sigma}}{  A (1 x 1) prior covariance matrix for respondent ideal points \eqn{x_t}.}
    \item{\code{beta$mu}}{  A (2 x 1) prior means matrix for \eqn{\tau_j} and \eqn{\beta^*}.}
    \item{\code{beta$sigma}}{  A (2 x 2) prior covariance matrix for \eqn{\tau_j} and \eqn{\beta^*}.}
     }
  }
  \item{.D}{ integer, indicates number of dimensions. Only one dimension is implemented and this argument is ignored.}
  \item{.control}{list, specifying some control functions for estimation.  Options include the following:
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
    An object of class \code{ordIRT}.

  \item{means}{ list, containing several matrices of point estimates for the parameters corresponding
    to the inputs for the priors.  The list should contain the following matrices.
    \itemize{
    \item{\code{x}}{  A (N x 1) matrix of point estimates for the respondent ideal points \eqn{x_t}.}
    \item{\code{beta}}{  A (J x 1) matrix of point estimates for the reparameterized item discrimination parameter \eqn{\beta^*}.}
    \item{\code{tau}}{  A (J x 1) matrix of point estimates for the bill cutpoint \eqn{\alpha^*}.}
    \item{\code{Delta_sq}}{  A (J x 1) matrix of point estimates for the squared bill cutpoint difference \eqn{\tau_j^2}.}
    \item{\code{Delta}}{  A (J x 1) matrix of point estimates for the bill cutpoint difference \eqn{\tau_j}.}
    }                 
   }

  \item{vars}{ list, containing several matrices of variance estimates for parameters corresponding
    to the inputs for the priors.  Note that these variances are those recovered via variational approximation,
    and in most cases they are known to be far too small and generally unusable.  Better estimates of variances
    can be obtained manually via the parametric bootstrap.  The list should contain the following matrices:
    \itemize{
    \item{\code{x}}{  A (N x 1) matrix of variances for the respondent ideal points \eqn{x_t}.}
    \item{\code{beta}}{  A (J x 1) matrix of variances for the reparameterized item discrimination parameter \eqn{\beta^*}.}
    \item{\code{tau}}{  A (J x 1) matrix of variances for the bill cutpoint \eqn{\tau_j}.}
    }                 
   }

  \item{runtime}{ A list of fit results, with elements listed as follows:}
    \itemize{
    \item{\code{iters}}{  integer, number of iterations run.}
    \item{\code{conv}}{  integer, convergence flag. Will return 1 if threshold reached, and 0 if maximum number of iterations reached.}
    \item{\code{threads}}{  integer, number of threads used to estimated model.}
    \item{\code{tolerance}}{  numeric, tolerance threshold for convergence.  Identical to thresh argument in input to .control list.}
    }                 
  \item{n}{ Number of respondents in estimation, should correspond to number of rows in roll call matrix.}
  \item{j}{ Number of items in estimation, should correspond to number of columns in roll call matrix.}
  \item{call}{ Function call used to generate output.}
}


\author{ 
Kosuke Imai \email{kimai@princeton.edu}

James Lo \email{jameslo@princeton.edu}

Jonathan Olmsted \email{jpolmsted@gmail.com}
}

\references{

Kosuke Imai, James Lo, and Jonathan Olmsted ``Fast Estimation of Ideal Points with Massive Data.''
Working Paper. Available at \url{http://imai.princeton.edu/research/fastideal.html}.

}

\seealso{
  '\link{AsahiTodai}'.
}

\examples{

\dontrun{
### Real data example: Asahi-Todai survey (not run)
## Collapses 5-category ordinal survey items into 3 categories for estimation
data(AsahiTodai)
out.varinf <- ordIRT(.rc = AsahiTodai$dat.all, .starts = AsahiTodai$start.values,
					.priors = AsahiTodai$priors, .D = 1,
					.control = {list(verbose = TRUE,
                     thresh = 1e-6, maxit = 500)})

## Compare against MCMC estimates using 3 and 5 categories
cor(ideal3, out.varinf$means$x)
cor(ideal5, out.varinf$means$x)  
}


### Monte Carlo simulation of ordIRT() model vs. known parameters
## Set number of legislators and items
set.seed(2)
NN <- 500
JJ <- 100

## Simulate true parameters from original model
x.true <- runif(NN, -2, 2)
beta.true <- runif(JJ, -1, 1)
tau1 <- runif(JJ, -1.5, -0.5)
tau2 <- runif(JJ, 0.5, 1.5)
ystar <- x.true \%o\% beta.true + rnorm(NN *JJ)

## These parameters are not needed, but correspond to reparameterized model
#d.true <- tau2 - tau1
#dd.true <- d.true^2
#tau_star <- -tau1/d.true
#beta_star <- beta.true/d.true

## Generate roll call matrix using simulated parameters
newrc <- matrix(0, NN, JJ)
for(j in 1:JJ) newrc[,j] <- cut(ystar[,j], c(-100, tau1[j], tau2[j],100), labels=FALSE)

## Generate starts and priors
cur <- vector(mode = "list")
cur$DD <- matrix(rep(0.5,JJ), ncol=1)
cur$tau <- matrix(rep(-0.5,JJ), ncol=1)
cur$beta <- matrix(runif(JJ,-1,1), ncol=1) 
cur$x <- matrix(runif(NN,-1,1), ncol=1) 
priors <- vector(mode = "list")
priors$x <- list(mu = matrix(0,1,1), sigma = matrix(1,1,1) )
priors$beta <- list(mu = matrix(0,2,1), sigma = matrix(diag(25,2),2,2))

## Call ordIRT() with inputs
time <- system.time({
    lout <- ordIRT(.rc = newrc,
                    .starts = cur,
                    .priors = priors,
                    .control = {list(
                        threads = 1,
                        verbose = TRUE,
                        thresh = 1e-6,
			maxit=300,
			checkfreq=50
                        )})
})

## Examine runtime and correlation of recovered ideal points vs. truth
time
cor(x.true,lout$means$x)

}
\keyword{ multivariate }
