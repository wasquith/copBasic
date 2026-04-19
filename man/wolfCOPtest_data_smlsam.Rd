\encoding{utf8}
\name{wolfCOPtest_data_smlsam}
\docType{data}
\alias{wolfCOPtest_data_smlsam}
\title{Small Sample Components for Schweizer--Wolff Sigma Test}
\description{
Data for the \emph{ test} components for small samples. Empirical sampling distributions of \emph{Schweizer--Wolff Sigma} (\eqn{\hat\sigma_\mathbf{C}}, \code{\link{wolfCOP}}) for \eqn{n \in 3, \cdots, 30} of the \eqn{\mathrm{NULL} : \mathbf{C} = \mathbf{\Pi}} using simulation sizes of 900,000 on the \eqn{\mathbf{C} = \mathbf{\Pi}} for each of the \eqn{n} concatenated table of empirical distributions were formed from the unique number of \eqn{\hat\sigma_\mathbf{C}} attainable up to \eqn{\pm4.26} standard deviations were computed and tabulated in the the \code{wolfCOPtest_data_smlsam} table. Algorithms for the small-sample components are located at \code{inst/}\code{make_wolfCOPtest/}\code{smlsam/} \code{make_}\code{wolfCOPtestP_}\code{smlsam.R}, and this script is the source of the \eqn{n = 30} in sources of function \code{\link{wolfCOPtest}}.
}
\usage{
data(wolfCOPtest_data_smlsam)
}
\format{
An \R \code{data.frame} titled \code{wolfCOPtest_data_smlsam} with
\describe{
  \item{n}{Sample size;}
  \item{probs}{Nonexceedance probability of empirical distribution of \eqn{\hat\sigma_\mathbf{C}} (\code{\link{wolfCOP}}) (note, the \code{quantile(..., type=4)} is used for the computation);}
  \item{wolfemp}{Value for \eqn{\hat\sigma_\mathbf{C}} (\code{\link{wolfCOP}});}
}
}
\examples{
\dontrun{
 data(wolfCOPtest_data_smlsam) # }
}
\keyword{datasets}
\keyword{Schweizer and Wolff Sigma}
