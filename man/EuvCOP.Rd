\encoding{utf8}
\name{EuvCOP}
\alias{EuvCOP}
\title{Expected value of U given V}
\description{
Compute the \emph{expected value} of \eqn{U} given a \eqn{V} through the emph{conditional distribution function} using the appropriate \emph{partial derivative} of a copula (\eqn{\mathbf{C}(u,v)}) with respect to \eqn{V}. The inversion of the partial derivative is the \emph{conditional quantile function}. Basic principles provide the expectation for a \eqn{x \ge 0} is
\deqn{E[X] = \int_0^\infty xf(x)\mathrm{d}x = \int_0^\infty \bigl(1-F_x(X)\bigr)\mathrm{d}x\mbox{,}}
which for the setting here becomes
\deqn{E[U|V = v] = \int_0^1 \bigl(1 - \frac{\delta}{\delta v} \mathbf{C}(u,v)\bigr)\mathrm{d}u\mbox{.}}
This function solves the integral using the \code{\link{derCOP2}} function.
}
\usage{
EuvCOP(v=seq(0.01, 0.99, by=0.01), cop=NULL, para=NULL, asuv=FALSE, ...)
}
\arguments{
  \item{v}{Nonexceedance probability \eqn{v} in the \eqn{Y} direction;}
  \item{cop}{A copula function with vectorization as in \code{asCOP};}
  \item{para}{Vector of parameters or other data structures, if needed, to pass to the copula;}
  \item{asuv}{Return a data frame of the \eqn{U} and \eqn{V}; and}
  \item{...}{Additional arguments to pass to \code{\link{derCOP2}}.}
}
\value{
  Value(s) for the expectation are returned.
}
\author{W.H. Asquith}
\seealso{\code{\link{EvuCOP}}, \code{\link{derCOP2}}}
\examples{
# See Examples under EvuCOP().
}
\keyword{copula (utility)}
\keyword{copula operator}
\keyword{copula (conditional distribution)}
\keyword{conditional expectation}
