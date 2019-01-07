\encoding{utf8}
\name{GLcop}
\alias{GLcop}
\alias{JOcopBB4}
\title{ The Galambos Extreme Value Copula (with Gamma Power Mixture [Joe/BB4] and Lower Extreme Value Limit) }
\description{
The \emph{Galambos copula} (Joe, 2014, p. 174) is
\deqn{\mathbf{C}_{\Theta}(u,v) = \mathbf{GL}(u,v) = uv\,\mathrm{exp}\bigl[(x^{-\Theta} + y^{-\Theta})^{-1/\Theta}\bigr]\mbox{,}}
where \eqn{\Theta \in [0, \infty)}, \eqn{x = -\log u}, and \eqn{y = -\log v}. As \eqn{\Theta \rightarrow 0^{+}}, the copula limits to \emph{independence} (\eqn{\mathbf{\Pi}}; \code{\link{P}}) and as  \eqn{\Theta \rightarrow \infty}, the copula limits to perfect association (\eqn{\mathbf{M}}; \code{\link{M}}). The copula here is a \emph{bivariate extreme value copula} (\eqn{BEV}), and parameter estimation requires for \eqn{\Theta} numerical methods.

There are two other genetically related forms. Joe (2014, p. 197) describes an extension of the Galambos copula as a \emph{Galambos gamma power mixture} (GPM), which is Joe's \emph{BB4 copula}, with the following form
\deqn{\mathbf{C}_{\Theta,\delta}(u,v) = \mathbf{GPM}(u,v) =
\biggl(x + y - 1 - \bigl[(x - 1)^{-\delta} + (y - 1)^{-\delta} \bigr]^{-1/\delta} \biggr)^{-1/\Theta}\mbox{,}}
where \eqn{x = u^{-\Theta}}, \eqn{y = v^{-\Theta}}, and \eqn{\Theta \ge 0, \delta \ge 0}. (Joe shows \eqn{\delta > 0}, but zero itself seems to work without numerical problems in practical application.) As \eqn{\delta \rightarrow 0^{+}}, the \dQuote{MTCJ family} results (implemented internally with \eqn{\Theta} as the incoming parameter). As \eqn{\Theta \rightarrow 0^{+}} the Galambos above results with \eqn{\delta} as the incoming parameter.

This second copula in turn has a \emph{lower extreme value limit form} that leads to a \emph{min-stable bivariate exponential} having \emph{Pickand's dependence function} of
\deqn{A(x,y; \Theta, \delta) = x + y - \bigl[x^{-\Theta} + y^{-\Theta} - (x^{\Theta\delta} + y^{\Theta\delta})^{-1/\delta} \bigr]^{-1/\Theta}\mbox{,}}
where this third copula is
\deqn{\mathbf{C}^{LEV}_{\Theta,\delta}(u,v) = \mathbf{GLEV}(u,v) =
\mathrm{exp}[-A(-\log u, -\log v; \Theta, \delta)]\mbox{,}}
for \eqn{\Theta \ge 0, \delta \ge 0} and is known as the \emph{two-parameter Galambos}. (Joe shows \eqn{\delta > 0}, but zero itself seems to work without numerical problems in practical application.)
}
\usage{
GLcop(u, v, para=NULL, ...)
JOcopBB4(u, v, para=NULL, ...)
}
\arguments{
  \item{u}{Nonexceedance probability \eqn{u} in the \eqn{X} direction;}
  \item{v}{Nonexceedance probability \eqn{v} in the \eqn{Y} direction;}
  \item{para}{To trigger \eqn{\mathbf{GL}(u,v)}, a vector (single element) of \eqn{\Theta}, to trigger \eqn{\mathbf{GLEV}(u,v)}, a two element vector of \eqn{\Theta} and \eqn{\delta}, and to trigger \eqn{\mathbf{GPM}(u,v)}, a three element vector of \eqn{\Theta}, \eqn{\delta}, and any number (the presence of the third entry alone is the triggering mechanism); and}
  \item{...}{Additional arguments to pass.}
}
\value{
  Value(s) for the copula are returned.
}
\references{
Joe, H., 2014, Dependence modeling with copulas: Boca Raton, CRC Press, 462 p.
}
\author{ W.H. Asquith}
\seealso{\code{\link{M}}, \code{\link{P}}, \code{\link{GHcop}}, \code{\link{HRcop}}, \code{\link{tEVcop}}
}
\examples{
# Theta = pi for GLcop and recovery through Blomqvist Beta    (Joe, 2014, p. 175)
log(2)/(log(log(2)/log(1+blomCOP(cop=GLcop, para=pi))))

# Theta = 2 and delta = 3 for the GPM form and Blomqvist Beta (Joe, 2014, p. 197)
t <- 2; Btheo <- blomCOP(GLcop, para=c(t,3,1))
Bform <- (2^(t+1) - 1 - taildepCOP(GLcop, para=c(t,3,1))$lambdaU*(2^t -1))^(-1/t)
print(c(Btheo, 4*Bform-1)) # [1] 0.8611903 0.8611900
}
\keyword{Galambos copula}
\keyword{Galambos copula (lower extreme value form)}
\keyword{copula (formulas)}
\keyword{copula}
\concept{Galambos extreme value copula}
\keyword{Joe (2014) Examples and Exercises}
\keyword{copula (extreme value)}
\keyword{extreme value copula}
\keyword{Joe/BB4 copula}
\concept{Joe/BB4 copula}
