\encoding{utf8}
\name{CLcop}
\alias{CLcop}
\title{ The Clayton Copula}
\description{
The Clayton copula (Joe, 2014, p. 176) is
\deqn{\mathbf{C}_{\Theta}(u,v) = \mathbf{CL}(u,v) = \mathrm{max}[(u^{-\Theta}+v^{-\Theta}-1; 0)]^{-1/\Theta}\mbox{,}}
where \eqn{\Theta \in [-1,\infty), \Theta \ne 0}
As \eqn{\Theta \rightarrow \infty^{+}}, the copula limits to , and the parameter \eqn{\Theta} is readily computed from Kendall's Tau:
\eqn{\tau = \frac{\Theta}{\Theta+2}\mbox{.}}
}
\usage{
CLcop(u, v, para=NULL, tau=NULL, ...)
}
\arguments{
  \item{u}{Nonexceedance probability \eqn{u} in the \eqn{X} direction;}
  \item{v}{Nonexceedance probability \eqn{v} in the \eqn{Y} direction;}
  \item{para}{A vector (single element) of parameters---the \eqn{\Theta} parameter of the copula;}
  \item{tau}{Optional Kendall's Tau; and}
  \item{...}{Additional arguments to pass.}
}
\value{
  Value(s) for the copula are returned otherwise if \code{tau} is given, then the \eqn{\Theta} is computed and a \code{list}
}
\references{
Joe, H., 2014, Dependence modeling with copulas: Boca Raton, CRC Press, 462 p.
}
\author{ W.H. Asquith}
\seealso{\code{\link{W}}
}
\examples{
\dontrun{
UV <- simCOP(10000, cop=CLcop, para=CLcop(tau=0.5)$para) #}
}
\keyword{Clayton copula}
\keyword{copula (formulas)}
\keyword{copula}
\concept{Clayton copula}
