\encoding{utf8}
\name{wolfCOPtest}
\alias{wolfCOPtest}
\title{Empirically Created Schweizer and Wolff Sigma Test Against an Independence Copula}
\description{
\emph{EXPERIMENTAL}---Perform a \emph{Schweizer--Wolff Sigma} test for \eqn{\sigma_\mathbf{C}} stemming from \eqn{\mathrm{NULL} : \mathbf{C} = \mathbf{\Pi}} (\emph{Independence copula}, \eqn{\mathbf{\Pi}}, \code{\link{P}}). The Schweizer--Wolff Sigma is defined under \code{\link{wolfCOP}}. This test has been highly tuned from small samples to large by extensive simulation; however, it is intimately tied to the version of sample computation by \code{\link{wolfCOP}} using a particular form of the empirical copula (refer to that function for details). The \code{wolfCOPtest} function internally loads the \code{\link{wolfCOPtest_data_smlsam}} tables, and most users will not need to access this dataset themselves unless for curiosity.
}
\usage{
wolfCOPtest(x, y, asuv=FALSE, aslist=TRUE, na.rm=TRUE, digits=6,
                  probs=c(0.90, 0.95, 0.98, 0.99, 0.995), ...)
}
\arguments{
  \item{x}{Nonexceedance probabilities \eqn{u} in the \eqn{X} direction if \code{asuv=TRUE}, else data in the \eqn{X} direction for which nonexceedance probabilities will be computed by plotting-position formula. If \code{x} is a single value, it is assumed to be a previously computed \eqn{\sigma_\mathbf{C}} and \code{y} will be considered as the sample size. Note, the statistical estimation is based on samples of \eqn{n \ge 3}, if the sample size is less than this, then internally \eqn{n = 3} is set for the computations;}
  \item{y}{Nonexceedance probabilities \eqn{v} in the \eqn{Y} direction if \code{asuv=TRUE}, else data in the \eqn{Y} direction for which nonexceedance probabilities will be computed by plotting-position formula (\code{lmomco::pp()} function);}
  \item{asuv}{Logical identifying whether the \code{x} and \code{y} are already in nonexceedance probabilities or rather such need to be computed internally by plotting-position formula (\code{lmomco::pp()} function);}
  \item{aslist}{Logical to return a list of the test results rather than a vector;}
  \item{na.rm}{A logical to trigger a \code{complete.cases} operation on internal matrix of the data as probabilities;}
  \item{digits}{Digits for rounding before return;}
  \item{probs}{Quantiles by nonexceedance probability to extract from the fitted distribution---these are the critical values in the upper tail. If \code{probs} is empty, then 0.95 nonexceedance probability is silently used; and}
  \item{...}{Additional arguments to pass to \code{lmomco::pp()}.}
}
\details{
The algorithm is based on simulation of \eqn{\mathbf{\Pi}} and computation of the L-moments of a \emph{logit transform} of the computed \eqn{\sigma_\mathbf{\Pi}} for which nonlinear regression of the L-moments to the sample size have been created. Development scripts based on Monte Carlo simulation are located in the directory \code{inst/make_}\code{wolfCOPtest/} of \pkg{copBasic} sources.

A logit transformation and its inversion are defined as
\deqn{\eta = \mathrm{log}\biggl(\frac{\sigma_\mathbf{C}}{1-\sigma_\mathbf{C}}\biggr)\mbox{\  and}
}
\deqn{
\eta^{(-1)} = \frac{1}{1 + \mathrm{exp}(-\sigma_\mathbf{C})}\mbox{.}
}

The predicted L-moments are used to drive either distributions from the \pkg{lmomco} package generalized normal distribution model (GNO) (\code{pargno}, \code{cdfgno}, \code{quagno}) or Pearson type III (PE3) (\code{parpe3}, \code{cdfpe3}, \code{quape3}) of the \emph{logit transform} of \eqn{\sigma_\mathbf{\Pi}} to compute the p-value and supporting information. Evaluation of an L-moment ratio diagram shows a heuristic change from GNO to PE3 at about \eqn{n = 40}, but inspection of the upper tail of the empirical distributions (and hence p-values in the range of those used for rejection in practice) suggest the PE3 tail is better. Therefore, the function exclusively uses the PE3. Sample-sample empirical distributions (\eqn{n \in 3, \cdots, 40}) will have p-values reported on output and missing for \eqn{n \ge 41} (refer to \code{\link{wolfCOPtest_data_smlsam}} for more details).
}
\value{
  Extensive list or named vector output the test is returned.
}
\author{W.H. Asquith}
\seealso{\code{\link{wolfCOP}}, \code{\link{wolfCOPsamc}}, \code{\link{wolfCOPtest_data_smlsam}}
}
\examples{
uv <- simCOP(n=62, cop=PSP, graphics=FALSE, seed=8)  # approx pvalue 0.001
wolfCOPtest(uv[,1], uv[,2], asuv=TRUE, aslist=FALSE)

uv <- simCOP(n=62, cop=P,   graphics=FALSE, seed=8)
wolfCOPtest(uv[,1], uv[,2], asuv=TRUE, aslist=FALSE) # approx pvalue 0.3
# Note, the small sample p-values are NA for the two previous computations because the
# sample size has become "large", but in the next example, we can have a look at the two
# p-values by the Pearson type III distribution to the predicted L-moments for the sample
# size and then the small-sample empirical tables in object wolfCOPtest_data_smlsam.
uv <- simCOP(n=13, cop=P,   graphics=FALSE, seed=8)
wolfCOPtest(uv[,1], uv[,2], asuv=TRUE, aslist=FALSE)
#   p.value(dist_pe3) p.value(sample_le40)
#           0.2529151            0.2493159

\dontrun{
  # This example has a long run time but will show that the Sigma is performing close to
  # Spearman Rho, which helps establish the veracity of the algorithms. If the user were
  # to switch out the copula P and replace with CIRCcop (circular copula), then one will
  # see Sigma performing similarly but Rho failing because the CIRCcop is an example of
  # perfect dependence in a setting of independence. Use of much larger nsim for this
  # example will hone the apparent rejection to about 0.05 for P but clearly we see
  # divergence in performance for CIRCcop.
  alpha <- 0.05; n <- 42; nsim <- 3000; cop <- P
  set.seed(1)
  wolf   <- replicate(nsim, { wolfCOP(para=simCOP(n=n, cop,
                              graphics=FALSE),  as.sample=TRUE) })
  wolfpv <- sapply(wolf, function(w) {
                   wolfCOPtest(w, n, asuv=TRUE, aslist=FALSE)[5] })
  message("Schweizer-Wolff Sigma(n=", n, ") type I error rate = ",
           sum(wolfpv < alpha) / nsim)
  set.seed(1)
  rhopv <- replicate(nsim, { uv <- simCOP(n=n, cop=cop, graphics=FALSE);
                   cor.test(uv[,1], uv[,2], method="spearman")$p.value })
  message("           Spearman Rho(n=", n, ") type I error rate = ",
           sum(rhopv  < alpha) / nsim)
  #       P : Schweizer-Wolff Sigma(n=42) type I error rate = 0.045
  #       P :          Spearman Rho(n=42) type I error rate = 0.043
  # CIRCcop : Schweizer-Wolff Sigma(n=42) type I error rate = 0.277
  # CIRCcop :          Spearman Rho(n=42) type I error rate = 0.0013
  # These example computations shown that Rho does does since the perfect dependence of
  # CIRCcop because its "correlation" is zero but Sigma has much high rejection rate.}
}
\keyword{copula (hypothesis test)}
\keyword{Schweizer and Wolff Sigma}
\keyword{generalized normal distribution}
\keyword{Pearson type III distribution}
\keyword{GNO (see generalized normal distribution)}
\keyword{PE3 (see Pearson type III distribution)}
\keyword{logit transform}
