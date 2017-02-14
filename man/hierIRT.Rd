\name{hierIRT}
\alias{hierIRT}
\title{ Hierarchichal IRT estimation via Variational Inference }
\description{ 
  \code{hierIRT} estimates an hierarchical IRT model with two response categories, allowing the use of covariates
  to help determine ideal point estimates. Estimation is conducted using the variational EM algorithm described
  in the reference paper below. A special case of this model occurs when time/session is used as the covariate ---
  this allows legislator ideal points to vary over time with a parametric time trend.  Notably, the popular
  DW-NOMINATE model (Poole and Rosenthal, 1997) is one such example, in which legislator ideal points shift
  by a constant amount each period, and the error term in the hierarchical model is set to 0.  In contrast to
  other functions in this package, this model does not assume a `rectangular' roll call matrix, and all data
  are stored in vector form.
}

\usage{
  hierIRT(.data, .starts = NULL, .priors = NULL, .control = NULL) 
}

\arguments{
  \item{.data}{ a list with the following items:
      \itemize{
      \item{\code{y}}{  A (L x 1) matrix of observed votes. `1' and `-1' are the yea and nay codes.}
      \item{\code{i}}{  A (L x 1) integer matrix of indexes of the ideal point i[l] linked to each observed
        vote l = 0 \eqn{\ldots} L.  Indexes begin at 0 and reach a maximum value of I - 1.}
      \item{\code{j}}{  A (L x 1) integer matrix of indexes of the bill/item j[l] linked to each observed
        vote l = 0 \eqn{\ldots} L.  Indexes begin at 0 and reach a maximum value of J - 1.}
      \item{\code{g}}{  A (I x 1) integer matrix of indexes of the group membership g[i[l]] linked to each
        ideal point i = 0 \eqn{\ldots} I.  Indexes begin at 0 and reach a maximum value of G - 1.}
      \item{\code{z}}{  A (I x D) numeric matrix of observed covariates.  Rows correspond to ideal points
        i = 0 \eqn{\ldots} I.  The columns correspond to the D different covariates.  Typically, the first
        column will be an intercept and fixed to 1, while other columns represent ideal point-specific covariates
        such as session.}
      }
    }
  \item{.starts}{ a list containing several matrices of starting values for the parameters. The list should
  	contain the following matrices:
    \itemize{
    \item{\code{alpha}}{  A (J x 1) matrix of starting values for the item difficulty parameter \eqn{\alpha_j}. }
    \item{\code{beta}}{  A (J x 1) matrix of starting values for the item discrimination parameter \eqn{\beta_j}. }
    \item{\code{gamma}}{  An (I x D) matrix of starting values for the group level coefficients \eqn{\gamma_{m}}. }
    \item{\code{eta}}{  An (I x 1) matrix of starting values for the ideal point error term \eqn{\eta_n}.}
    \item{\code{sigma}}{  An (G x 1) matrix of starting values for the group level variance parameter \eqn{\sigma^2_m}. }    
    }
  }
  \item{.priors}{  list, containing several matrices of starting values for the parameters. The list should
    contain the following matrices:
    \itemize{
    \item{\code{gamma.mu}}{  A (D x 1) prior means matrix for all group level coefficients \eqn{\gamma_{m}}.}
    \item{\code{gamma.sigma}}{  A (D x D) prior covariance matrix for all group level coefficients \eqn{\gamma_{m}}.}
    \item{\code{beta.mu}}{  A (2 x 1) prior means matrix for all bill parameters \eqn{\alpha_j} and \eqn{\beta_j}.}
    \item{\code{beta.sigma}}{  A (2 x 2) prior covariance matrix for all bill parameters \eqn{\alpha_j} and \eqn{\beta_j}.}
    \item{\code{sigma.v}}{  A (1 x 1) matrix containing the group level variance prior parameter \eqn{\nu_{\sigma}}.}
    \item{\code{sigma.s}}{  A (1 x 1) matrix containing the group level variance prior parameter \eqn{s^2_{\sigma}}.}
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
    An object of class \code{hierIRT}.

  \item{means}{ list, containing several matrices of point estimates for the parameters corresponding
    to the inputs for the priors.  The list should contain the following matrices.
    \itemize{
    \item{\code{alpha}}{  A (J x 1) matrix of point estimates for the item difficulty parameter \eqn{\alpha_j}. }
    \item{\code{beta}}{  A (J x 1) matrix of point estimates for the item discrimination parameter \eqn{\beta_j}. }
    \item{\code{gamma}}{  An (I x D) matrix of point estimates for the group level coefficients \eqn{\gamma_{m}}. }
    \item{\code{eta}}{  An (I x 1) matrix of point estimates for the ideal point error term \eqn{\eta_n}.}
    \item{\code{sigma}}{  An (G x 1) matrix of point estimates for the group level variance parameter \eqn{\sigma^2_m}. }    
    \item{\code{x_implied}}{  An (I x 1) matrix of the implied ideal point \eqn{x_i}, calculated as a function of
      gamma, z, and eta using the point estimates for those parameters. }    
    }                 
   }
  \item{vars}{ list, containing several matrices of variance estimates for several parameters of interest for diagnostic
    purposes. Note that these variances are those recovered via variational approximation, and in most cases they are
    known to be far too small and generally unusable. The list should contain the following matrices:
    \itemize{
    \item{\code{eta}}{  A (I x 1) matrix of variance estimates for the ideal point noise parameter \eqn{\eta_n}.}
    \item{\code{gamma}}{  A (G x D x D) cube of covariance estimates for the gamma parameters for each group.
      Each of the G items is a matrix with a single covariance matrix for the m-th group's D gamma parameters.}
    \item{\code{beta2}}{  A (J x 2 x 2) cube of covariance estimates for the item parameters \eqn{\alpha_j} and \eqn{\beta_j}.
      Each of the J items is a matrix with a single covariance matrix for the j-th item.}
    }                 
   }

  \item{runtime}{ A list of fit results, with elements listed as follows:}
    \itemize{
    \item{\code{iters}}{  integer, number of iterations run.}
    \item{\code{conv}}{  integer, convergence flag. Will return 1 if threshold reached, and 0 if maximum number of iterations reached.}
    \item{\code{threads}}{  integer, number of threads used to estimated model.}
    \item{\code{tolerance}}{  numeric, tolerance threshold for convergence.  Identical to thresh argument in input to .control list.}
    }                 

  \item{N}{ A list of counts of various items:}
    \itemize{
    \item{\code{D}}{  integer, number of dimensions (i.e. number of covariates, including intercept).}
    \item{\code{G}}{  integer, number of groups.}
    \item{\code{I}}{  integer, number of ideal points.}
    \item{\code{J}}{  integer, number of items/bill parameters.}
    \item{\code{L}}{  integer, number of observed votes.}    
    }                 

  \item{call}{ Function call used to generate output.}
}


\author{ 
Kosuke Imai \email{kimai@princeton.edu}

James Lo \email{lojames@usc.edu}

Jonathan Olmsted \email{jpolmsted@gmail.com}
}

\references{

Variational model is described in Kosuke Imai, James Lo, and Jonathan Olmsted ``Fast Estimation
of Ideal Points with Massive Data.'' Working Paper.
Available at \url{http://imai.princeton.edu/research/fastideal.html}.

}

\seealso{
  '\link{dwnom}'.
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
\keyword{ multivariate }
