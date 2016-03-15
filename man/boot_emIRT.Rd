\name{boot_emIRT}
\alias{boot_emIRT}
\title{ Parametric bootstrap of EM Standard Errirs }
\description{ 
  \code{boot_emIRT} take an emIRT() object (from binary, ordinal, dynamic, or hierarchical models) and implements
  a parametric bootstrap of the standard errors for the ideal points.  It assumes you have already run the model
  successfully on one of those functions, and takes the output from that estimate, along with the original arguments,
  as the arguments for the bootstrap function.
}

\usage{
  boot_emIRT(emIRT.out, .data, .starts, .priors, .control, Ntrials=50, verbose=10)
}

\arguments{
  \item{emIRT.out}{  an emIRT() object, which is output from a call to binIRT(), dynIRT(), ordIRT(), or hierIRT()}
  \item{.data}{ the data used to produce the emIRT object. }
  \item{.starts}{ the starts used to produce the emIRT object.  }
  \item{.priors}{ the priors used to produce the emIRT object. }
  \item{.control}{ the control arguments used to produce the emIRT object. }
  \item{Ntrials}{ Number of bootstrap trials to run. }
  \item{verbose}{ Number of trials before progress notification triggers. }
}

\value{
    An object of class \code{emIRT}.  The output takes the original emIRT.out object and appends the following:
  \item{bse}{ list, containing only the matrix:
    \itemize{
    \item{\code{x}}{  A (N x 1) matrix of bootstrapped standard errors}
    }                 
   }
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
  '\link{binIRT}', '\link{ordIRT}', '\link{hierIRT}', '\link{dynIRT}'.
}

\examples{

\dontrun{

### Binary IRT example
example(binIRT)
boot.bin <- boot_emIRT(lout, .data = rc, .starts = s, .priors = p,
  .control = list(threads = 1, verbose = FALSE, thresh = 1e-06), Ntrials=10, verbose=2)
boot.bin$bse$x

### Dynamic IRT example
example(dynIRT)
boot.dyn <- boot_emIRT(lout, .data = mq_data$data.mq, .starts = mq_data$cur.mq,
  .priors = mq_data$priors.mq, .control = list(threads = 1, verbose = FALSE,
  thresh = 1e-06), Ntrials=10, verbose=2)
boot.dyn$bse$x

### Ordinal IRT example
example(ordIRT)
boot.ord <- boot_emIRT(lout, .data=newrc, .starts=cur, .priors=priors,
  .control = list(threads = 1, verbose = TRUE, thresh = 1e-6, maxit=300,
  checkfreq=50), Ntrials=5, verbose=1)
boot.ord$bse$x

### Hierarhical IRT example
example(hierIRT, run.dontrun=TRUE)
boot.hier <- boot_emIRT((lout, .data=dwnom$data.in, .starts=dwnom$cur, .priors=dwnom$priors,
	.control=list(threads = 8, verbose = TRUE, thresh = 1e-4, maxit=200, checkfreq=1),
	Ntrials=5, verbose=1)
boot.hier$bse$x_implied

}

}
\keyword{ multivariate }
