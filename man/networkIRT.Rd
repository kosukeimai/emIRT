\name{networkIRT}
\alias{networkIRT}
\title{ Network IRT estimation via EM }
\description{ 
  \code{networkIRT} estimates an IRT model with network in cells. Estimation
  is conducted using the EM algorithm described in the reference paper below. The algorithm
  generalizes a model by Slapin and Proksch (2009) that is commonly applied to manifesto
  data.
}

\usage{
   networkIRT(.y, .starts = NULL, .priors = NULL, .control = NULL,
    .anchor_subject = NULL, .anchor_item = NULL)
}

\arguments{
  \item{.y}{  matrix, with 1 indicating a valid link and 0 otherwise.  Followers (usually voters) are on
     rows, elites are on columns.  No NA values are permitted.}
  \item{.starts}{ a list containing several matrices of starting values for the parameters. 
      The list should contain the following matrices:
    \itemize{
    \item{\code{alpha}}{  A (J x 1) matrix of starting values for politician propensity to be followed \eqn{alpha}.}
    \item{\code{beta}}{  A (N x 1) matrix of starting values for follower propensity to follow others \eqn{\beta}.}
    \item{\code{w}}{  An (J x 1) matrix of starting values for politician ideal points \eqn{z}.}
    \item{\code{theta}}{  An (N x 1) matrix of starting values for the follower ideal points \eqn{x}.}
    \item{\code{gamma}}{  An (1 x 1) matrix, should generally be fixed to be 1.}
    }
  }
  \item{.priors}{  list, containing several matrices of starting values for the parameters.
     The list should contain the following matrices (1x1) matrices:
    \itemize{
    \item{\code{alpha$mu}}{  prior mean for \eqn{\alpha}.}
    \item{\code{alpha$sigma}}{  prior variance for \eqn{\alpha}}
    \item{\code{beta$mu}}{ prior mean for \eqn{\beta}.}
    \item{\code{beta$sigma}}{ prior variance for \eqn{\beta}.}
    \item{\code{w$mu}}{ prior mean for z.}
    \item{\code{w$sigma}}{ prior variance for z}
    \item{\code{theta$mu}}{ prior mean for x.}
    \item{\code{theta$sigma}}{ prior variance for x.}
    \item{\code{gamma$mu}}{  Should be fixed to equal 1.}
    \item{\code{gamma$sigma}}{  Should be fixed to equal 1.}

     }
  }
  \item{.control}{list, specifying some control functions for estimation.  Options include the following:
    \itemize{
    \item{\code{threads}}{  integer, indicating number of cores to use. Default is to use a single core,
        but more can be supported if more speed is desired.}
    \item{\code{verbose}}{  boolean, indicating whether output during estimation should be verbose or not.
        Set FALSE by default.}
    \item{\code{thresh}}{  numeric. Algorithm will run until all parameters correlate at 1 - thresh across
    	consecutive iterations. Set at 1e-6 by default.}
    \item{\code{maxit}}{  integer. Sets the maximum number of iterations the algorithm can run.
        Set at 500 by default.}
    \item{\code{checkfreq}}{  integer. Sets frequency of verbose output by number of iterations.
        Set at 50 by default.}
    }
  }
  \item{.anchor_subject}{integer, specifying subject to use as identification anchor.}
    \item{.anchor_item}{integer, specifying item to use as identification anchor.}
}

\value{
    An object of class \code{networkIRT}.

  \item{means}{ list, containing several matrices of point estimates for the parameters corresponding
    to the inputs for the priors.  The list should contain the following matrices.
    \itemize{
    \item{\code{alpha}}{  A (J x 1) matrix of point estimates for politician propensity to be followed \eqn{alpha}.}
    \item{\code{beta}}{  A (N x 1) matrix of point estimates for follower propensity to follow others \eqn{\beta}.}
    \item{\code{w}}{  An (J x 1) matrix of point estimates for politician ideal points \eqn{z}.}
    \item{\code{theta}}{  An (N x 1) matrix of point estimates for the follower ideal points \eqn{x}.}                
   }
  }
  \item{vars}{ list, containing several matrices of variance estimates for parameters corresponding
    to the inputs for the priors.  Note that these variances are those recovered via variational approximation,
    and in most cases they are known to be far too small and generally unusable.  Better estimates of variances
    can be obtained manually via the parametric bootstrap.  The list should contain the following matrices:
    \itemize{
    \item{\code{alpha}}{  A (J x 1) matrix of variance estimates for politician propensity to be followed \eqn{alpha}.}
    \item{\code{beta}}{  A (N x 1) matrix of variance estimates for follower propensity to follow others \eqn{\beta}.}
    \item{\code{w}}{  An (J x 1) matrix of variance estimates for politician ideal points \eqn{z}.}
    \item{\code{theta}}{  An (N x 1) matrix of variance estimates for the follower ideal points \eqn{x}.}                
    }                 
   }

  \item{runtime}{ A list of fit results, with elements listed as follows:}
    \itemize{
    \item{\code{iters}}{  integer, number of iterations run.}
    \item{\code{conv}}{  integer, convergence flag. Will return 1 if threshold reached, and 0
       if maximum number of iterations reached.}
    \item{\code{threads}}{  integer, number of threads used to estimated model.}
    \item{\code{tolerance}}{  numeric, tolerance threshold for convergence.  Identical to thresh
       argument in input to .control list.}
    }                 
  \item{N}{ Number of followers in estimation, should correspond to number of rows in data matrix .y}
  \item{J}{ Number of politicians in estimation, should correspond to number of columns in data matrix .y}
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
  '\link{ustweet}'
}

\examples{

\dontrun{
data(ustweet)

## A ridiculously short run to pass CRAN
## For a real test, set maxit to a more reasonable number to reach convergence
lout <- networkIRT(.y = ustweet$data,
                   .starts = ustweet$starts,
                   .priors = ustweet$priors,
                   .control = {list(verbose = TRUE,
                                    maxit = 3,
                                    convtype = 2,
                                    thresh = 1e-6,
                                    threads = 1
                                    )
                           },
                   .anchor_item = 43
                   )

}

}
\keyword{ multivariate }
