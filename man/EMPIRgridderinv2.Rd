\encoding{utf8}
\name{EMPIRgridderinv2}
\alias{EMPIRgridderinv2}
\title{Derivative Inverses of the Grid of the Bivariate Empirical Copula for U with respect to V}
\description{
Generate a gridded representation of the inverse of the derivatives of the \emph{bivariate empirical copula} of \eqn{U} with respect to \eqn{V}.  This function is the empirical analog to \code{\link{derCOPinv2}}.
}
\usage{
EMPIRgridderinv2(empgrid=NULL, kumaraswamy=FALSE, dergrid=NULL, ...)
}
\arguments{
  \item{empgrid}{The grid from \code{\link{EMPIRgrid}};}
  \item{kumaraswamy}{A logical to trigger Kumaraswamy smoothing of the conditional quantile function;}
  \item{dergrid}{The results of \code{\link{EMPIRgridder2}} and if left \code{NULL} then that function is called internally. There is some fragility at times in the quality of the numerical derivative and the author has provided this argument so that the derivative can be computed externally and then fed to this inversion function; and}
  \item{...}{Additional arguments to pass.}
}
\value{
  The gridded values of the inverse of the derivative of \eqn{U} with respect to \eqn{V}.
}
\author{W.H. Asquith}
\seealso{\code{\link{EMPIRcop}}, \code{\link{EMPIRcopdf}}, \code{\link{EMPIRgrid}}, \code{\link{EMPIRgridder2}}, \code{\link{EMPIRgridderinv}}, \code{\link{EMPIRgridderinv2}}}
\examples{
# See examples under EMPIRgridderinv
}
\keyword{empirical copula}
\keyword{empirical copula (derivative inverses)}
\keyword{derivative}
