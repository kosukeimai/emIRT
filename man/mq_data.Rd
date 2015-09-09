\name{mq_data}
\alias{mq_data}

\title{
  Martin-Quinn Judicial Ideology Scores
}

\description{

  This data set contains materials related to the Martin-Quinn measures of judicial ideology,
  estimated for every justice serving from October 1937 to October 2013.  The materials are
  based on the July 23, 2014 (2014 Release 01) release of the Supreme Court Database, which
  contain the votes of each Supreme Court justice on each case heard in the court.  The
  data is set up for input to \code{dynIRT()}, and also includes point estimates of the same
  model obtained using standard MCMC techniques.  Start values and priors input to this model
  are identical to those used in the MCMC estimates, and were provided by Kevin Quinn. 

}

\usage{
data(mq_data)
}

\value{

  \item{mq_data}{ list, containing the following elements:
    \itemize{
    \item{\code{data.mq}}{ Justice voting data, formatted for input to \code{dynIRT()}.}
    \item{\code{cur.mq}}{ Start values, formatted for input to \code{dynIRT()}. }
    \item{\code{priors.mq}}{ Priors, formatted for input to \code{dynIRT()}.}
    \item{\code{mq_mcmc}}{ Ideal point estimates with data via MCMC.}
    \item{\code{justiceName}}{ A vector of names identifying the justice that goes with each estimated ideal point.}
    }
  }
}

\references{

Original model and the example is based off of Andrew Martin and Kevin Quinn, ``Dynamic Ideal
Point Estimation via Markov Chain Monte Carlo for the U.S. Supreme Court, 1953-1999.'' Political
Analysis 10(2) 134-153.  Also see \url{http://mqscores.berkeley.edu}.

Variational model is described in Kosuke Imai, James Lo, and Jonathan Olmsted ``Fast Estimation
of Ideal Points with Massive Data.'' Working Paper. Available at
\url{http://imai.princeton.edu/research/fastideal.html}.

}

\seealso{
  '\link{dynIRT}'.
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

\keyword{datasets}
