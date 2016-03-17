\name{poisIRT}
\alias{poisIRT}
\title{ Poisson IRT estimation via EM }
\description{ 
  \code{poisIRT} estimates an IRT model with count (usually word counts) in cells. Estimation
  is conducted using the EM algorithm described in the reference paper below. The algorithm
  generalizes a model by Slapin and Proksch (2009) that is commonly applied to manifesto
  data.
}

\usage{
    poisIRT(.rc, i = 0:(nrow(.rc)-1), NI = nrow(.rc), .starts = NULL, .priors = NULL,
    .control = NULL)
}

\arguments{
  \item{.rc}{  matrix, usually with unique words along the J rows and different documents across K columns.
     Each cell will contain a count of words.  There should be no NA values, so documents missing a particular
     word should list 0 in the cell.}
  \item{i}{  vector of length K, indicating for each of the K documents which actor it belongs to. Assignment
     of actors begins at actor 0.  If set to 0:(K-1), and NI=K below, then each document
     is assigned its own ideal point, and we get the Slapin and Proksch Wordfish model.}
  \item{NI}{  integer, number of unique actors.  Must be less than or equal to K.  If NI=K, then each document
     is assigned its own ideal point, and we get the Slapin and Proksch Wordfish model.}
  \item{.starts}{ a list containing several matrices of starting values for the parameters. 
      The list should contain the following matrices:
    \itemize{
    \item{\code{alpha}}{  A (J x 1) matrix of starting values for the word frequency parameter \eqn{alpha}.}
    \item{\code{psi}}{  A (K x 1) matrix of starting values for the document verbosity parameter \eqn{psi}.}
    \item{\code{beta}}{  A (J x 1) matrix of starting values for the word discrimination parameter \eqn{\beta}.}
    \item{\code{x}}{  An (NI x 1) matrix of starting values for the actor ideal points \eqn{x_i}.}
    }
  }
  \item{.priors}{  list, containing several matrices of starting values for the parameters.
     The list should contain the following matrices:
    \itemize{
    \item{\code{x$mu}}{  numeric, prior mean for actor ideal points \eqn{x_i}.}
    \item{\code{x$sigma2}}{   numeric, prior variance for actor ideal points \eqn{x_i}.}
    \item{\code{beta$mu}}{  numeric, prior mean for \eqn{\beta_j}.}
    \item{\code{beta$sigma2}}{  numeric, prior variance for \eqn{\beta_j}.}
    \item{\code{alpha$mu}}{  numeric, prior mean for \eqn{\alpha_j}.}
    \item{\code{alpha$sigma2}}{  numeric, prior variance for \eqn{\alpha_j}}
    \item{\code{psi$mu}}{  numeric, prior mean for \eqn{\psi_k}.}
    \item{\code{psi$sigma2}}{  numeric, prior variance for \eqn{\psi_k}.}

     }
  }
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
    An object of class \code{poisIRT}.

  \item{means}{ list, containing several matrices of point estimates for the parameters corresponding
    to the inputs for the priors.  The list should contain the following matrices.
    \itemize{
    \item{\code{alpha}}{  A (J x 1) matrix of point estimates for the word frequency parameter \eqn{alpha}.}
    \item{\code{psi}}{  A (K x 1) matrix of point estimates for the document verbosity parameter \eqn{psi}.}
    \item{\code{beta}}{  A (J x 1) matrix of point estimates for the word discrimination parameter \eqn{\beta}.}
    \item{\code{x}}{  An (NI x 1) matrix of point estimates for the actor ideal points \eqn{x_i}.}
    }                 
   }

  \item{vars}{ list, containing several matrices of variance estimates for parameters corresponding
    to the inputs for the priors.  Note that these variances are those recovered via variational approximation,
    and in most cases they are known to be far too small and generally unusable.  Better estimates of variances
    can be obtained manually via the parametric bootstrap.  The list should contain the following matrices:
    \itemize{
    \item{\code{beta}}{  A (J x 1) matrix of variational variance estimates for the word discrimination parameter \eqn{\beta}.}
    \item{\code{x}}{  An (NI x 1) matrix of variational variance estimates for the actor ideal points \eqn{x_i}.}
    }                 
   }

  \item{runtime}{ A list of fit results, with elements listed as follows:
    \itemize{
    \item{\code{iters}}{  integer, number of iterations run.}
    \item{\code{conv}}{  integer, convergence flag. Will return 1 if threshold reached, and 0 if maximum number of iterations reached.}
    \item{\code{threads}}{  integer, number of threads used to estimated model.}
    \item{\code{tolerance}}{  numeric, tolerance threshold for convergence.  Identical to thresh argument in input to .control list.}
    }                 
  }

  \item{N}{ A list of sizes, with elements listed as follow:
    \itemize{
     \item{K}{ Number of unique words in term-document matrix.}
     \item{J}{ Number of documents in term-document matrix.}
     \item{I}{ Number of actors in model, always less than or equal to J.}
     \item{call}{ Function call used to generate output.}
	}
  }

   \item{i_of_k}{ A copy of input for argument `i', which allows the J documents to be linked to I actors.}
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
  '\link{manifesto}'
}

\examples{

\dontrun{
## Load German Manifesto data
data(manifesto)

## Estimate variational Wordfish model
lout <- poisIRT(.rc = manifesto$data.manif,
		    i = 0:(ncol(manifesto$data.manif)-1),
		    NI=ncol(manifesto$data.manif),
                    .starts = manifesto$starts.manif,
                    .priors = manifesto$priors.manif,
                    .control = {list(
                        threads = 1,
                        verbose = TRUE,
                        thresh = 1e-6,
						maxit=1000
                        )})

## Positional Estimates for Parties
lout$means$x  
}

}
\keyword{ multivariate }
