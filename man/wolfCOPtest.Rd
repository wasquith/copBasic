\encoding{utf8}
\name{wolfCOPtest}
\alias{wolfCOPtest}
\title{A Schweizer and Wolff Sigma Test Against an Independence Copula}
\description{
\emph{EXPERIMENTAL}---Perform a \emph{Schweizer--Wolff Sigma} test for \eqn{\sigma_\mathbf{C}} stemming from \eqn{\mathrm{NULL} : \mathbf{C} = \mathbf{\Pi}} (\emph{Independence copula}, \eqn{\mathbf{\Pi}}, \code{\link{P}}). The Schweizer--Wolff Sigma is defined under \code{\link{wolfCOP}}. The algorithm is based on simulation of \eqn{\mathbf{\Pi}} and computation of the L-moments of a logit transform of the computed \eqn{\sigma_\mathbf{\Pi}} for which nonlinear regression of the L-moments to the sample size have been created as used through a generalized normal distribution model (GNO) (\code{lmomco::pargno}, \code{lmomco::cdfgno}, \code{lmomco::quagno}) of the \emph{logit transform} of \eqn{\sigma_\mathbf{\Pi}} to compute the p-value and supporting information.

Development scripts based on Monte Carlo simulation are located at \code{inst/make_wolfCOPtest/} of \pkg{copBasic} sources. A logit transformation is defined as
\deqn{\eta = \mathrm{log}\biggl(\frac{\sigma_\mathbf{C}}{1-\sigma_\mathbf{C}}\biggr)\mbox{\  and}
}
\deqn{
\eta^{(-1)} = \frac{1}{1 + \mathrm{exp}(-\sigma_\mathbf{C})}\mbox{.}
}
}
\usage{
wolfCOPtest(x, y, asuv=FALSE, aslist=TRUE, bylogit=TRUE, dtype="gno",
                  na.rm=TRUE, digits=6,
                  probs=c(0.90, 0.95, 0.98, 0.99, 0.995), ...)
}
\arguments{
  \item{x}{Nonexceedance probabilities \eqn{u} in the \eqn{X} direction if \code{asuv=TRUE}, else data in the \eqn{X} direction for which nonexceedance probabilities will be computed by plotting-position formula. If \code{x} is a single value, it is assumed to be a previously computed \eqn{\sigma_\mathbf{C}} and \code{y} will be considered as the sample size. Note, the statistical estimation is based on samples of \eqn{n \ge 9}, if the sample size is less than this, then internally \eqn{n = 9} is set for the computations;}
  \item{y}{Nonexceedance probabilities \eqn{v} in the \eqn{Y} direction if \code{asuv=TRUE}, else data in the \eqn{Y} direction for which nonexceedance probabilities will be computed by plotting-position formula (\code{lmomco::pp()} function);}
  \item{asuv}{Logical identifying whether the \code{x} and \code{y} are already in nonexceedance probabilities or rather such need to be computed internally by plotting-position formula (\code{lmomco::pp()} function);}
  \item{aslist}{Logical to return a list of the test results rather than a vector;}
  \item{bylogit}{Compute results using a logit transform of the \eqn{\sigma_\mathbf{C}} distribution through the L-moments of the logit in lieu of remaining in \dQuote{sigma} units. This argument is subject to deprecation as \code{inst/make_wolfCOPtest/} development to implementation provides evidence for preference of working through a \bold{logit} that honors the \eqn{\sigma_\mathbf{C} \in (0, 1)};}
  \item{dtype}{Distribution in the abbreviation recognized by the \pkg{lmomco} package that will be used to estimate quantiles and p-value based on fitted L-moments whether based on logit transformation (\code{bylogit=TRUE}) or just no transformation (\code{bylogit=FALSE});}
  \item{na.rm}{A logical to trigger a \code{complete.cases} operation on internal matrix of the data as probabilities;}
  \item{digits}{Digits for rounding before return;}
  \item{probs}{Quantiles by nonexceedance probability to extract from the fitted distribution---these are the critical values in the upper tail. If \code{probs} is empty, then 0.95 nonexceedance probability is silently used; and}
  \item{...}{Additional arguments to pass to \code{lmomco::pp}.}
}
\value{
  The value for \eqn{\sigma_\mathbf{C}} is returned.
}
\author{W.H. Asquith}
\seealso{\code{\link{wolfCOP}}, \code{\link{wolfCOPsamc}}
}
\examples{
uv <- simCOP(n=62, cop=PSP, graphics=FALSE, seed=6)
wolfCOPtest(uv[,1], uv[,2], asuv=TRUE, bylogit=TRUE,  aslist=FALSE)
wolfCOPtest(uv[,1], uv[,2], asuv=TRUE, bylogit=FALSE, aslist=FALSE)

uv <- simCOP(n=62, cop=P,   graphics=FALSE, seed=6)
wolfCOPtest(uv[,1], uv[,2], asuv=TRUE, bylogit=TRUE,  aslist=FALSE)
wolfCOPtest(uv[,1], uv[,2], asuv=TRUE, bylogit=FALSE, aslist=FALSE)
}
\keyword{copula (hypothesis test)}
\keyword{generalized normal distribution}
\keyword{GNO (see generalized normal distribution)}
\keyword{Schweizer and Wolff Sigma}
