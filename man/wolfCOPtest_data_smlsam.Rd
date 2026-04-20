\encoding{utf8}
\name{wolfCOPtest_data_smlsam}
\docType{data}
\alias{wolfCOPtest_data_smlsam}
\title{Small Sample Components for Schweizer--Wolff Sigma Test}
\description{
Data for the \emph{ test} components for small samples. Empirical sampling distributions of \emph{Schweizer--Wolff Sigma} (\eqn{\hat\sigma_\mathbf{C}}, \code{\link{wolfCOP}}) for \eqn{n \in 3, \cdots, 40} of the \eqn{\mathrm{NULL} : \mathbf{C} = \mathbf{\Pi}} using simulation sizes of 900,000 for each \eqn{n} on \eqn{\mathbf{C} = \mathbf{\Pi}}. So, total simulation size was 34,200,000. Continuing, this data set is a concatenated table of empirical distributions formed from the unique number of \eqn{\hat\sigma_\mathbf{C}} (discrete variables in small samples), which large-simulation attainable up to \eqn{\pm4.26} standard deviations, were computed and tabulated in the \code{wolfCOPtest_data_smlsam} table. Algorithms for the small-sample components are located at \code{inst/}\code{make_wolfCOPtest/}\code{smlsam/} \code{make_}\code{wolfCOPtestP_}\code{smlsam.R}, and this script is the source of the \eqn{n = 40} in sources of function \code{\link{wolfCOPtest}}. Finally, for the smallest of sample sizes, the \emph{empirical distribution} is a step function, and for each plateau, a nudge using \code{.Machine$double.eps} has been used to make \code{approx()} for sensible in \code{\link{wolfCOPtest}}. In the \bold{Examples} herein, the jaggedness of the empirical distribution can readily be seen and visually smoothing by about \eqn{n > 10}. Finally, the \code{\link{wolfCOPtest}} function internally loads the \code{wolfCOPtest_data_smlsam} tables, and most users will not need to access this dataset themselves unless for curiosity.
}
\usage{
data(wolfCOPtest_data_smlsam)
}
\format{
An \R \code{data.frame} titled \code{wolfCOPtest_data_smlsam} with
\describe{
  \item{n}{Sample size;}
  \item{probs}{Nonexceedance probability of empirical distribution of \eqn{\hat\sigma_\mathbf{C}} (\code{\link{wolfCOP}}). Note, \code{quantile(..., type=4)} was used for computations because of structural similarity to sample \eqn{\hat\sigma_\mathbf{C}} as implemented in \pkg{copBasic};}
  \item{wolfemp}{Value for \eqn{\hat\sigma_\mathbf{C}} (\code{\link{wolfCOP}});}
}
}
\examples{
data(wolfCOPtest_data_smlsam)

\dontrun{
  W <- seq(0, 1, by=0.005) # Span all Schweizer-Wolff Sigma. If the first argument to
  # wolfCOPtest is scalar and so is the second, these are treated as sigma and the
  # sample size. The parametric dist. is the generalized normal for all these sample
  # sizes and switches to Pearson type III at n > 41 (refer to inst/make_wolfCOPtest/).
  for(i in 3:40) { # samples sizes 3 to 40 in wolfCOP_data_smlsam
    Pd <- sapply(W, function(w) wolfCOPtest(w, i)$p.value[1])
    Pe <- sapply(W, function(w) wolfCOPtest(w, i)$p.value[2])
    # p-value by the generalized normal (parametric) distribution
    plot( W, Pd, type="l", col="grey50", lwd=7, main=i)
    # p-value derived from wolfCOP_data_smlsam tables:
    lines(W, Pe, col="darkgreen", lwd=2) # small sample p-values
  } #}
}
\keyword{datasets}
\keyword{Schweizer and Wolff Sigma}
