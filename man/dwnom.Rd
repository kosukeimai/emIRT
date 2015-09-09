\name{dwnom}
\alias{dwnom}

\title{
  Poole-Rosenthal DW-NOMINATE data and scores, 80-110 U.S. Senate
}

\description{

  This data set contains materials related to the Poole-Rosenthal DW-NOMINATE measure of senator
  ideology.  The software (and other materials) is available at \url{ftp://voteview.com/dw-nominate.htm},
  which includes a simpler example of an application to the 80-110 U.S. Senate.  The data set here is
  derived from the data and estimates from that example, but are formatted to be run in \code{hierIRT()}.
  In particular, start values for estimation are identical to those provided by the example.

}

\usage{
data(dwnom)
}

\value{

  \item{dwnom}{ list, containing the following elements:
    \itemize{
    \item{\code{data.in}}{ Legislator voting data, formatted for input to \code{hierIRT()}.}
    \item{\code{cur}}{ Start values, formatted for input to \code{hierIRT()}. }
    \item{\code{priors}}{ Priors, formatted for input to \code{hierIRT()}.}
    \item{\code{legis}}{ data frame, containing contextual information about the legislators estimated.}
    \item{\code{nomres}}{ data frame, containing estimates from DW-NOMINATE on the same data. These
      are read from the file SL80110C21.DAT. }
    }
  }
}

\references{

DW-NOMINATE is described in Keith T. Poole and Howard Rosenthal. 1997. Congress: A Political Economic History
of Roll Call Voting. Oxford University Press.  See also \url{ftp://voteview.com/dw-nominate.htm}.

Variational model is described in Kosuke Imai, James Lo, and Jonathan Olmsted ``Fast Estimation
of Ideal Points with Massive Data.'' Working Paper. Available at
\url{http://imai.princeton.edu/research/fastideal.html}.

}

\seealso{
  '\link{hierIRT}'.
}

\examples{

### Real data example of US Senate 80-110 (not run)
### Based on voteview.com example of DW-NOMINATE (ftp://voteview.com/dw-nominate.htm)
### We estimate a hierarchical model without noise and a linear time covariate
### This model corresponds very closely to the DW-NOMINATE model

\dontrun{
data(dwnom)

## This takes about 10 minutes to run on 8 threads
## You may need to reduce threads depending on what your machine can support
lout <- hierIRT(.data = dwnom$data.in,
                    .starts = dwnom$cur,
                    .priors = dwnom$priors,
                    .control = {list(
                    threads = 8,
                    verbose = TRUE,
                    thresh = 1e-4,
				    maxit=200,
				    checkfreq=1
                        )})

## Bind ideal point estimates back to legislator data
final <- cbind(dwnom$legis, idealpt.hier=lout$means$x_implied)

## These are estimates from DW-NOMINATE as given on the Voteview example
## From file "SL80110C21.DAT"
nomres <- dwnom$nomres

## Merge the DW-NOMINATE estimates to model results by legislator ID
## Check correlation between hierIRT() and DW-NOMINATE scores
res <- merge(final, nomres, by=c("senate","id"),all.x=TRUE,all.y=FALSE)
cor(res$idealpt.hier, res$dwnom1d)

}
}

\keyword{datasets}
