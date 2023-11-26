\encoding{utf8}
\name{bicoploc}
\alias{bicoploc}
\title{Analog to Line of Organic Correlation by Copula Diagonal}
\description{
Compute an analog to the \emph{line of organic correlation} (\emph{reduced major axis}) using the \emph{diagonal} of a copula.
}
\usage{
bicoploc(x, y, xout=NA, xpara=NULL, ypara=NULL, dtypex="nor", dtypey="nor",
               ctype=ctype=c("weibull", "hazen", "bernstein", "checkerboard"),
               kumaraswamy=TRUE, plotonxy=TRUE, plotondiag=TRUE, a=0,
               tt=pnorm(seq(-5, +5, by=0.1)), ...)
}
\arguments{
  \item{x}{Numeric vector giving data points of \eqn{X};}
  \item{y}{Numeric vector giving data points of \eqn{Y};}
  \item{xout}{An optional set of numeric values specifying where interpolation through the diagonal inversion is to take place;}
  \item{xpara}{An \pkg{lmomco} package parameter object for the \eqn{X} variable, which if not provided will trigger an method of L-moment parameter estimation for the distribution \code{dtypex};}
  \item{ypara}{An \pkg{lmomco} package parameter object for the \eqn{Y} variable, which if not provided will trigger an method of L-moment parameter estimation for the distribution \code{dtypey};}
  \item{dtypex}{The \pkg{lmomco} package distribution abbreviations (see \code{lmomco::dist.list()}) for the \eqn{X} variable;}
  \item{dtypey}{The \pkg{lmomco} package distribution abbreviations (see \code{lmomco::dist.list()}) for the \eqn{Y} variable;}
  \item{breve}{An alternative way from \code{para} to set the \eqn{\beta} for this function; and}
  \item{...}{Additional arguments to pass to the copula.}
}
\value{
  List(s) for the computations and predictions are returned.
}
\references{
Chang, B., and Joe, H., 2020, Copula diagnostics for asymmetries and conditional
dependence: Journal of Applied Statistics, v. 47, no. 9, pp. 1587--1615, \doi{10.1080/02664763.2019.1685080}.
}
\author{A copBasic User}
\seealso{\code{\link{diagCOPatf}}}
\examples{
set.seed(1); nsim <- 50
X <- rnorm(nsim, mean=3, sd=0.6)
Y <- 0.843*X + rnorm(nsim, mean=0, sd=0.2)
zz <- bicoploc(X,Y, xout=c(2.5, 3.5, 4), dtypex="nor", dtypey="nor")
# cor(X,Y, method="spearman") # for nsim=1E6 RHO=0.9237038

set.seed(1); nsim <- 50
X <- rlmomco(nsim, vec2par(c(3, 0.6, +0.5), type="pe3"))
Y <- 0.3*X^2 + rlmomco(nsim, vec2par(c(0, 0.4, 0), type="pe3"))
zz <- bicoploc(X,Y, xout=c(2.5, 3.5, 4.5), dtypex="nor", dtypey="nor")
# cor(X,Y, method="spearman") # for nsim=1E6 RHO=0.9238148

set.seed(1); nsim <- 50
X <- rlmomco(nsim, vec2par(c(3, 0.6, +0.5), type="pe3"))
Y <- 0.3*X^2 + rlmomco(nsim, vec2par(c(0, 0.4, 0), type="pe3"))
zz <- bicoploc(X,Y, xout=c(2.5, 3.5, 4.5), dtypex="gev", dtypey="gev")
}
\keyword{bivariate prediction}
\keyword{copula (diagonal)}
\keyword{copula (diagonal inversion)}
