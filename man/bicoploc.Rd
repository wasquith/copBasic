\encoding{utf8}
\name{bicoploc}
\alias{bicoploc}
\title{Analog to Line of Organic Correlation by Copula Diagonal}
\description{
\emph{HIGHLY EXPERIMENTAL AND SUBJECT TO OVERHAUL OR REMOVAL}---Compute an analog to the \emph{line of organic correlation} (\emph{reduced major axis}) using the \emph{diagonal} of a copula.

Two diagnostic plots can be plotted by the arguments available for this function. The plot for \eqn{(U,V)} coordinate nonexceedance probability domain along with the analyses involving the copula diagonal inversion comes first, which is followed by that for \eqn{(X,Y)} coordinate domain along with the well-known line of organic correlation by the method of L-moments and such a \dQuote{line} (could be a curve) by copula diagonal inversion.
}
\usage{
bicoploc(x, y, xout=NA, xpara=NULL, ypara=NULL, dtypex="nor", dtypey="nor",
               ctype=c("weibull", "hazen", "bernstein", "checkerboard"),
               kumaraswamy=TRUE, plotonuv=TRUE, plotonxy=TRUE, adduv=FALSE, addxy=FALSE,
               a=0, ff=pnorm(seq(-5, +5, by=0.1)), ...)
}
\arguments{
  \item{x}{Numeric vector giving data points of \eqn{X}. If this is a matrix or data frame, then the first and second columns are extracted for the \code{x} and \code{y} internall;}
  \item{y}{Optional numeric vector giving data points of \eqn{Y} depending on the composition of \code{x};}
  \item{xout}{An optional set of numeric values specifying where interpolation through the diagonal inversion is to take place;}
  \item{xpara}{An \pkg{lmomco} package parameter object for the \eqn{X} variable, which if not provided will trigger an method of L-moment parameter estimation for the distribution \code{dtypex};}
  \item{ypara}{An \pkg{lmomco} package parameter object for the \eqn{Y} variable, which if not provided will trigger an method of L-moment parameter estimation for the distribution \code{dtypey};}
  \item{dtypex}{The \pkg{lmomco} package distribution abbreviations (see \code{lmomco::dist.list()}) for the \eqn{X} variable;}
  \item{dtypey}{The \pkg{lmomco} package distribution abbreviations (see \code{lmomco::dist.list()}) for the \eqn{Y} variable;}
  \item{ctype}{Argument of the same name for the empirical copula for dispatch to \code{\link{EMPIRcop}}, and note that the \code{1/n} form is disabled for \code{bicoploc} operations based on limited experiments;}
  \item{kumaraswamy}{A logical to trigger Kumaraswamy distribution smoothing of the copula diagonal inversion from the \emph{empricial copula}. The Kumaraswamy distribution is a distribution having support \eqn{[0,1]} with an explicit quantile function and takes the place of a Beta distribution (see \pkg{lmomco} function \code{quakur()} for more details). The smoothing by Kumaraswamy will provide a continuous real number on the interval, which should insure no flat-lining as one rolls on to or off off the discrete interval provided by the empirical copula for expected sample sizes of the operations anticipated for the \code{bicoploc} function;}
  \item{plotonuv}{A logical to trigger plotting of the analyses in the \eqn{(U,V)} coordinate nonexceedance probability domain along with the analyses involving the copula diagonal inversion. If set true, then \code{adduv} is set false internally;}
  \item{plotonxy}{A logical to trigger plotting of the analyses in the \eqn{(X,Y)} coordinate domain along with the well-known line of organic correlation by the method of L-moments and such a \dQuote{line}  (could be a curve) by copula diagonal inversion. If set true, then \code{addxy} is set false internally;}
  \item{adduv}{A logical when set true will not call the \code{plot()} function for \eqn{(U,V)} coordinate nonexceedance probability domain but the other graphical operations of lines and points will be called;}
  \item{addxy}{A logical when set true will not call the \code{plot()} function for \eqn{(X,Y)} coordinate domain but the other graphical operations of lines and points will be called;}
  \item{a}{Value for the plotting-position formula for \code{lmomco::pp()}, default is \code{a=0}, which specifies \emph{Weibull plotting positions};}
  \item{ff}{The nonexceedance joint probability of of the copula diagonal from which inversion computes the marginal nonexceedance probability values for \eqn{u=v=t}; and}
  \item{...}{Additional arguments to pass.}
}
\value{
  List(s) for the computations and predictions are returned.
}
\note{
The use of a copula diagonal inversion for purposes of a line of organic correlation analog within the text-book literature and elsewhere is unknown to the developers (December 2023).
}
\references{
Kruskal, W.H., 1953, On the uniqueness of the line of organic correlation:
Biometrics, vol.~9, no.~1, pp. 47--58, \doi{10.2307/3001632}.
}
\author{W.H. Asquith}
\seealso{\code{\link{diagCOPatf}}}
\examples{
set.seed(1); nsim <- 50
X  <- rnorm(nsim, mean=3, sd=0.6)
Y  <- rnorm(nsim, mean=0, sd=0.2)
zz <- bicoploc(X,Y, xout=c(2.5, 3.5, 4), dtypex="nor", dtypey="nor")

set.seed(1); nsim <- 50
X  <-           rnorm(nsim, mean=3, sd=0.6)
Y  <- 0.843*X + rnorm(nsim, mean=0, sd=0.2)
zz <- bicoploc(X,Y, xout=c(2.5, 3.5, 4), dtypex="nor", dtypey="nor")
# cor(X,Y, method="spearman") # for nsim=1E6 RHO=0.92367

set.seed(1); nsim <- 50
X  <-           lmomco::rlmomco(nsim, lmomco::vec2par(c(3, 0.6, +0.5), type="pe3"))
Y  <- 0.3*X^2 + lmomco::rlmomco(nsim, lmomco::vec2par(c(0, 0.4, 0),    type="pe3"))
zz <- bicoploc(X,Y, xout=c(2.5, 3.5, 4.5), dtypex="nor", dtypey="nor")
# cor(X,Y, method="spearman") # for nsim=1E6 RHO=0.92366

set.seed(1); nsim <- 50
X  <-           lmomco::rlmomco(nsim, lmomco::vec2par(c(3, 0.6, +0.5), type="pe3"))
Y  <- 0.3*X^2 + lmomco::rlmomco(nsim, lmomco::vec2par(c(0, 0.4, 0),    type="pe3"))
zz <- bicoploc(X,Y, xout=c(2.5, 3.5, 4.5), dtypex="gev", dtypey="gev")

###################################################################################
# COUNTER EXAMPLE? Slope of the LOC is negative. What should or should not be done?
set.seed(1); nsim <- 50
UV <- rCOP(nsim, cop=breveCOP, para=list(cop=W, alpha=0.4, beta=0.3))
X  <- qnorm(UV[,1], mean=3, sd=0.6)
Y  <- qnorm(UV[,2], mean=2, sd=0.2)
zz <- bicoploc(X,Y, xout=c(1.5, 2.5, 3.5, 4), dtypex="gno", dtypey="gno")
}
\keyword{bivariate prediction}
\keyword{copula (diagonal)}
\keyword{copula (diagonal inversion)}
\keyword{line of organic correlation}
\concept{organic correlation}
\concept{reduced major axis}
