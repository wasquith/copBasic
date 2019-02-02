\encoding{utf8}
\name{RFcop}
\alias{RFcop}
\title{ The Raftery Copula }
\description{
The \emph{Raftery copula} (Nelsen, 2006, p. 172) is
\deqn{\mathbf{C}_{\Theta}(u,v) = \mathbf{RF}(u,v) = \mathbf{M}(u,v) + \frac{1-\Theta}{1+\Theta}(uv)^{1/(1-\Theta)}\bigl[1-(\mathrm{max}\{u,v\})^{-(1+\Theta)/(1-\Theta)}\bigr]\mbox{,}}
where \eqn{\Theta \in (0,1)}. The copula, as \eqn{\Theta \rightarrow 0^{+}} limits, to the \emph{independence coupla} (\eqn{\mathbf{P}(u,v)}; \code{\link{P}}), and as \eqn{\Theta \rightarrow 1^{-}}, limits to the \emph{comonotonicity copula} (\eqn{\mathbf{M}(u,v)};  \code{\link{M}}). The parameter \eqn{\Theta} is readily computed from \emph{Spearman's Rho} (\code{\link{rhoCOP}}) by \eqn{\rho_\mathbf{C} = \Theta(4-3\Theta)/(2-\Theta)^2} or from \emph{Kendall's Tau} (\code{\link{tauCOP}}) by \eqn{\tau_\mathbf{C} = 2\Theta/(3-\Theta)}.
}
\usage{
RFcop(u, v, para=NULL, rho=NULL, tau=NULL, ...)
}
\arguments{
  \item{u}{Nonexceedance probability \eqn{u} in the \eqn{X} direction;}
  \item{v}{Nonexceedance probability \eqn{v} in the \eqn{Y} direction;}
  \item{para}{A vector (single element) of parameters---the \eqn{\Theta} parameter of the copula;}
  \item{rho}{Optional Spearmans's Rho;}
  \item{tau}{Optional Kendall's Tau; and}
  \item{...}{Additional arguments to pass.}
}
\value{
  Value(s) for the copula are returned. Otherwise if either \code{rho} or \code{tau} is given, then the \eqn{\Theta} is computed and a \code{list} having
  \item{para}{The parameter \eqn{\Theta}, and}
  \item{rho}{Spearmans's Rho if the \code{rho} is given.}
  \item{tau}{Kendall's Tau if the \code{tau} is given but also if both \code{rho} and \code{tau} are \code{NULL} as mentioned next.}
and if \code{para=NULL} and \code{rho} and \code{tau=NULL}, then the values within \code{u} and \code{v} are used to compute Kendall's Tau and then compute the parameter, and these are returned in the aforementioned list.
}
\references{
Nelsen, R.B., 2006, An introduction to copulas: New York, Springer, 269 p.
}
\author{ W.H. Asquith}
\seealso{\code{\link{M}}, \code{\link{P}}
}
\examples{
# Lower tail dependency of Theta = 0.5 --> 2*(0.5)/(1+0.5) = 2/3 (Nelsen, 2006, p. 214)
taildepCOP(cop=RFcop, para=0.5)$lambdaL # 0.66667

\dontrun{
# Simulate for a Spearman Rho of 0.7, then extract estimated Theta that
# internally is based on Kendall Tau of U and V, then convert estimate
# to equivalent Rho.
UV <- simCOP(1000, cop=RFcop, RFcop(rho=0.7)$para)
Theta <- RFcop(UV$U, UV$V)$para      # 0.605093
Rho <- Theta*(4-3*Theta)/(2-Theta)^2 # 0.679403 (nearly 0.7)#}
}
\keyword{Raftery copula}
\keyword{copula (formulas)}
\keyword{copula}
\concept{Raftery copula}
\keyword{Nelsen (2006) Examples and Exercises}