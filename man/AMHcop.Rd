\encoding{utf8}
\name{AMHcop}
\alias{AMHcop}
\title{ The Ali--Mikhail--Haq Copula}
\description{
The \emph{Ali--Mikhail--Haq copula} (Joe, 2014, p. 168) is
\deqn{\mathbf{C}_{\Theta}(u,v) = \mathbf{AMH}(u,v) = \frac{uv}{1 - \Theta(1-u)(1-v)}\mbox{,}}
where \eqn{\Theta \in [-1,+1)}. The copula \eqn{\Theta \rightarrow 0} becomes the \emph{independence copula} (\eqn{\mathbf{\Pi}(u,v)}; \code{\link{P}}). Finally, the parameter \eqn{\Theta} is readily computed from \emph{Kendall's Tau} (\code{\link{tauCOP}}) by
\deqn{\tau_\mathbf{C} = \frac{3\Theta - 2}{3\Theta} -   
                        \frac{2(1-\Theta)^2\log(1-\Theta)}{3\Theta^2}\mbox{,}}
and \emph{Spearman Rho} (\code{\link{rhoCOP}}) by
\deqn{\rho_\mathbf{C} = 12\biggl(\frac{\mathrm{dilog}(1-\Theta)\times(1+\Theta) - 2\log(1-\Theta) + 2\Theta\log(1-\Theta)-3\Theta}{\Theta^2}\biggl)\; -\; 3\mbox{,}}
where \eqn{\mathrm{dilog}(x)} is the \emph{dilogarithm} function (\emph{Spence's function}):
\deqn{\mathrm{dilog}(x) = \int_1^x \frac{\log(t)}{1-t}\, \mathrm{d}t\mbox{.}}
The support of \eqn{\tau_\mathbf{C}} is \eqn{\approx [-0.1817258, 0.3333333]} and the \eqn{\rho_\mathbf{C}} is \eqn{\approx [-0.2710647, 0.4784176]}, which shows that this copula has a limited range of dependency.
}
\usage{
AMHcop(u, v, para=NULL, tau=NULL, ...)
}
\arguments{
  \item{u}{Nonexceedance probability \eqn{u} in the \eqn{X} direction;}
  \item{v}{Nonexceedance probability \eqn{v} in the \eqn{Y} direction;}
  \item{para}{A vector (single element) of parameters---the \eqn{\Theta} parameter of the copula. However, if a second parameter is present, it is treated as a logical to reverse the copula (\eqn{u + v - 1 + \mathbf{AMH}(1-u,1-v; \Theta)});}
  \item{tau}{Optional Kendall's Tau; and}
  \item{...}{Additional arguments to pass.}
}
\value{
  Value(s) for the copula are returned. Otherwise if \code{tau} is given, then the \eqn{\Theta} is computed and a \code{list} having
  \item{para}{The parameter \eqn{\Theta}, and}
  \item{tau}{Kendall's tau.}
and if \code{para=NULL} and \code{tau=NULL}, then the values within \code{u} and \code{v} are used to compute Kendall Tau and then compute the parameter, and these are returned in the aforementioned list.
}
\references{
Pranesh, Kumar, 2010, Probability distributions and estimation of Ali--Mikhail--Haq copula: Applied Mathematical Sciences, v. 4, no. 14, p. 657--666.
}
\author{ W.H. Asquith}
\seealso{\code{\link{P}}
}
\examples{
\dontrun{
t <- 0.9 # The Theta of the copula and we will compute Spearman Rho.
di <- integrate(function(t) log(t)/(1-t), lower=1, upper=(1-t))$value
A <- di*(1+t) - 2*log(1-t) + 2*t*log(1-t) - 3*t
rho <- 12*A/t^2 - 3    # 0.4070369
rhoCOP(AMHcop, para=t) # 0.4070369}

\dontrun{
layout(matrix(1:2,byrow=TRUE)) # Note Kendall's Tau same on reversal.
UVn <- simCOP(3560, cop=AMHcop, para=c(-0.9, "FALSE" ), col=4)
mtext("Normal definition [default]") # '2nd' parameter could be skipped
UVr <- simCOP(3560, cop=AMHcop, para=c(-0.9, "TRUE"),   col=2)
mtext("Reversed definition")
AMHcop(UVn[,1], UVn[,2])$tau # Should be about -0.1663313 because
AMHcop(UVr[,1], UVr[,2])$tau # tauCOP(AMHcop, para=-0.9) is theoretically. }
}
\keyword{Ali-Mikhail-Haq copula}
\keyword{copula (formulas)}
\keyword{copula}
\concept{Ali-Mikhail-Haq copula}
