\encoding{utf8}
\name{wolfCOPtest}
\alias{wolfCOPtest}
\title{Empirically Created Schweizer and Wolff Sigma Test Against an Independence Copula}
\description{
\emph{EXPERIMENTAL}---Perform a \emph{Schweizer--Wolff Sigma} test for \eqn{\sigma_\mathbf{C}} stemming from \eqn{\mathrm{NULL} : \mathbf{C} = \mathbf{\Pi}} (\emph{Independence copula}, \eqn{\mathbf{\Pi}}, \code{\link{P}}). The Schweizer--Wolff Sigma is defined under \code{\link{wolfCOP}}.
}
\usage{
wolfCOPtest(x, y, asuv=FALSE, aslist=TRUE, na.rm=TRUE, digits=6,
                  probs=c(0.90, 0.95, 0.98, 0.99, 0.995), ...)
}
\arguments{
  \item{x}{Nonexceedance probabilities \eqn{u} in the \eqn{X} direction if \code{asuv=TRUE}, else data in the \eqn{X} direction for which nonexceedance probabilities will be computed by plotting-position formula. If \code{x} is a single value, it is assumed to be a previously computed \eqn{\sigma_\mathbf{C}} and \code{y} will be considered as the sample size. Note, the statistical estimation is based on samples of \eqn{n \ge 9}, if the sample size is less than this, then internally \eqn{n = 9} is set for the computations;}
  \item{y}{Nonexceedance probabilities \eqn{v} in the \eqn{Y} direction if \code{asuv=TRUE}, else data in the \eqn{Y} direction for which nonexceedance probabilities will be computed by plotting-position formula (\code{lmomco::pp()} function);}
  \item{asuv}{Logical identifying whether the \code{x} and \code{y} are already in nonexceedance probabilities or rather such need to be computed internally by plotting-position formula (\code{lmomco::pp()} function);}
  \item{aslist}{Logical to return a list of the test results rather than a vector;}
  \item{na.rm}{A logical to trigger a \code{complete.cases} operation on internal matrix of the data as probabilities;}
  \item{digits}{Digits for rounding before return;}
  \item{probs}{Quantiles by nonexceedance probability to extract from the fitted distribution---these are the critical values in the upper tail. If \code{probs} is empty, then 0.95 nonexceedance probability is silently used; and}
  \item{...}{Additional arguments to pass to \code{lmomco::pp}.}
}
\details{
The algorithm is based on simulation of \eqn{\mathbf{\Pi}} and computation of the L-moments of a logit transform of the computed \eqn{\sigma_\mathbf{\Pi}} for which nonlinear regression of the L-moments to the sample size have been created. Development scripts based on Monte Carlo simulation are located at \code{inst/make_}\code{wolfCOPtest/} of \pkg{copBasic} sources.

A logit transformation and its inversion are defined as
\deqn{\eta = \mathrm{log}\biggl(\frac{\sigma_\mathbf{C}}{1-\sigma_\mathbf{C}}\biggr)\mbox{\  and}
}
\deqn{
\eta^{(-1)} = \frac{1}{1 + \mathrm{exp}(-\sigma_\mathbf{C})}\mbox{.}
}

The predicted L-moments are used to drive either generalized normal distribution model (GNO) (\code{lmomco::pargno}, \code{lmomco::cdfgno}, \code{lmomco::quagno}) or Pearson type III (PE3) (\code{lmomco::parpe3}, \code{lmomco::cdfpe3}, \code{lmomco::quape3}) of the \emph{logit transform} of \eqn{\sigma_\mathbf{\Pi}} to compute the p-value and supporting information.

The switch between approximating distributions occurs at relatively small sample sizes (see internals of \code{wolfCOPtest}). Evidence exists of the decision between these two parametric model based on L-moment ratio diagram evidence and inspection of upper tail from exact sampling distributions (\eqn{n \in 3, \cdots, 25}). The small-sample algorithms are located at \code{inst/} \code{make_wolfCOPtest/}\code{smlsam/} \code{make_}\code{wolfCOPtestP_}\code{smlsam.R}).
}
\value{
  The value for \eqn{\sigma_\mathbf{C}} is returned.
}
\author{W.H. Asquith}
\seealso{\code{\link{wolfCOP}}, \code{\link{wolfCOPsamc}}
}
\examples{
uv <- simCOP(n=62, cop=PSP, graphics=FALSE, seed=6)
wolfCOPtest(uv[,1], uv[,2], asuv=TRUE, aslist=FALSE)

uv <- simCOP(n=62, cop=P,   graphics=FALSE, seed=6)
wolfCOPtest(uv[,1], uv[,2], asuv=TRUE, aslist=FALSE)

\dontrun{
  alpha <- 0.05; n <- 50; nsim <- 60000
  set.seed(1)
  wolf   <- replicate(nsim, { wolfCOP(para=simCOP(n=n, P, graphics=FALSE),
            as.sample=TRUE) })
  wolfpv <- sapply(wolf, function(w) {
                   wolfCOPtest(w, n, asuv=TRUE, aslist=FALSE)[5] })
  message("Schweizer–Wolff Sigma(n=", n, ") type II error rate = ",
           1 - sum(wolfpv < alpha)/nsim)
  set.seed(1)
  rhopv <- replicate(nsim, { uv <- simCOP(n=n, cop=P, graphics=FALSE);
                   cor.test(uv[,1], uv[,2], method="spearman")$p.value })
  message("           Spearman Rho(n=", n, ") type II error rate = ",
           1 - sum(rhopv < alpha)/nsim) #}

\dontrun{
  W <- seq(0, 1, by=0.005)
  for(i in 3:35) {
    Pd <- sapply(W, function(w) wolfCOPtest(w, i)$p.value[1])
    Pe <- sapply(W, function(w) wolfCOPtest(w, i)$p.value[2])
    plot( W, Pd, type="l", main=i) # pval by parametric distribution
    lines(W, Pe, col="red") # small sample from wolfCOP_data_smlsam
  } #}
}
\keyword{copula (hypothesis test)}
\keyword{Schweizer and Wolff Sigma}
\keyword{generalized normal distribution}
\keyword{Pearson type III distribution}
\keyword{GNO (see generalized normal distribution)}
\keyword{PE3 (see Pearson type III distribution)}
