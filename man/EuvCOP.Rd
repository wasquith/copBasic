\encoding{utf8}
\name{EuvCOP}
\alias{EuvCOP}
\title{Expected value of U given V}
\description{
Compute the \emph{expected value} of \eqn{U} given a \eqn{V} through the emph{conditional distribution function} using the appropriate \emph{partial derivative} of a copula (\eqn{\mathbf{C}(u,v)}) with respect to \eqn{V}. The inversion of the partial derivative is the \emph{conditional quantile function}. Basic principles provide the expectation for a \eqn{x \ge 0} is

\deqn{E[X] = \int_0^\infty xf(x)\mathrm{d}x = \int_0^\infty \bigl(1-F_x(X)\bigr)\mathrm{d}x\mbox{,}}

which for the setting here becomes

\deqn{E[U|V = v] = \int_0^1 \bigl(1 - \frac{\delta}{\delta v} \mathbf{C}(u,v)\bigr)\mathrm{d}u\mbox{.}}

This function solves the integral using the \code{\link{derCOP2}} function. This avoids a call of the \code{\link{derCOPinv2}} through its \code{uniroot()} inversion of the derivative. The example shown for \code{EuvCOP()} below does a validation check using conditional simulation, which is dependence (of course) of the design of the \pkg{copBasic} package, as part of simple isolation of a horizontal slice of the simulation and computing the mean of the \eqn{V} within the slice.
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
\note{
The author is well aware that the name of this function does not contain the number \code{2} as the family of functions also sharing this \emph{with respect to \eqn{v}} nature. It was a design mistake in 2008 to have used the \code{2}. The \code{uv} in the function name is the moniker for this \emph{with respect to \eqn{v}}.
}
\author{W.H. Asquith}
\seealso{\code{\link{EvuCOP}}, \code{\link{derCOP2}}}
\examples{
# We can show algorithmic validation using a highly asymmetric case of a
# copula having its parameter also nearly generating a singular component.
v <- c(0.2, 0.8); n <- 5E2; set.seed(1)
para <- list(cop=HRcop, para=120, alpha=0.4, beta=0.05)
UV   <- simCOP(n, cop=composite1COP, para=para, graphics=FALSE) # set TRUE to view

sapply(v, function(vv) EuvCOP(vv, cop=composite1COP, para=cpara))
# [1] 0.3168844 0.6925082

sapply(v, function(vv) mean( UV$U[UV$V > vv - 50/n & UV$V < vv + 50/n] ) )
# [1] 0.2796518 0.7092755 # general validation is thus shown as n-->large

# If visualized, we see in the lower corner than heuristics suggest a mean further
# to the right of the "singularity" for v=0.2 than v=0.80. For v=0.80, the
# "singularity" appears tighter given the upper tail dependency contrast of the
# coupla in the symmetrical case (alpha=0, beta=0) and changing the parameter to
# a Spearman Rho (say) similar to the para settting in this example. So, 0.70 for
# the mean given v=0.80 makes sense. Further notice that the two estimates of the
# mean are further apparent for v=0.2 relative to v=0.80. Again, this makes sense
# when the copula is visualized even at small n let along large.

# See also Examples under EvuCOP().

\dontrun{
  n <- 5000; Vlo <- rep(0.001, n); Vhi <- 1 - Vlo; Theta <- 3
  Ulo <- simCOPmicro(Vlo, cop=JOcopB5, para=Theta); dlo <- density(Ulo)
  Uhi <- simCOPmicro(Vhi, cop=JOcopB5, para=Theta); dhi <- density(Uhi)
  summary(Ulo)
  Ulomu <- EuvCOP(Vlo[1], cop=JOcopB5, para=Theta)
  summary(Uhi)
  Uhimu <- EuvCOP(Vhi[1], cop=JOcopB5, para=Theta)

  UV <- simCOP(n, cop=JOcopB5, para=Theta, pch=21, bg="lightgreen", col="darkgreen")
  abline(h=Vlo[1], col="salmon") # near the bottom to then form axis for density
  abline(h=Vhi[1], col="purple") # near the   top  to then form axis for density
  lines(dlo$x,   dlo$y/max(dlo$y), col="salmon", lwd=2) # density of the simulation
  # along the line already drawn at the bottom (Vlo)
  lines(dhi$x, 1-dhi$y/max(dhi$y), col="purple", lwd=2) # density of the simulation
  # along the line already drawn at the   top  (Vhi)
  points(Ulomu, Vlo[1], pch=16,  col="salmon", cex=2)
  points(Uhimu, Vhi[1], pch=16,  col="purple", cex=2)
  uv <- seq(0.001, 0.999, by=0.001)
  lines(uv, EuvCOP(u, cop=JOcopB5, para=Theta),     col="blue", lwd=3.5, lty=2)
  lines(    EvuCOP(u, cop=JOcopB5, para=Theta), uv, col="blue", lwd=1.5, lty=2) #}
}
\keyword{copula (utility)}
\keyword{copula operator}
\keyword{copula (conditional distribution)}
\keyword{conditional expectation}
