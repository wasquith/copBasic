\encoding{utf8}
\name{Wshuf}
\alias{Wshuf}
\title{Shuffles of \enc{Fréchet}{Frechet}--Hoeffding Lower-Bound Copula}
\description{
Compute shuffles of \emph{\enc{Fréchet}{Frechet}--Hoeffding lower-bound copula} (Nelsen, 2006, p. 173), which is defined by partitioning \eqn{\mathbf{W}} within \eqn{\mathcal{I}^2} into \eqn{n} subintervals:
\deqn{\mathbf{W}_n(u,v) = \mathrm{max}\biggl(\frac{k-1}{n}, u+v-\frac{k}{n} \biggr)}
for points within the partitions
\deqn{(u,v) \in \biggl[\frac{k-1}{n}, \frac{k}{n}\biggr]\times \biggl[ \frac{k-1}{n}, \frac{k}{n}\biggr]\mbox{,\ }k = 1,2,\cdots,n}
and for points otherwise out side the partitions
\deqn{\mathbf{W}_n(u,v) = \mathrm{min}(u,v)\mbox{.}}
The support of \eqn{\mathbf{W}_n} consists of \eqn{n} line segments connecting coordinate pairs \eqn{\{(k-1)/n,\, k/n\}} and \eqn{\{k/n,\, (k-1)/n\}} as stated by Nelsen (2006). The \emph{Spearman Rho} (\code{\link{rhoCOP}}) is defined by \eqn{\rho_\mathbf{C} = 1 - (2/n^2)}, and the \emph{Kendall Tau} (\code{\link{tauCOP}}) by \eqn{\tau_\mathbf{C} = 1 - (2/n)}.
}
\usage{
Wshuf(u, v, para=1, ...)
}
\arguments{
  \item{u}{Nonexceedance probability \eqn{u} in the \eqn{X} direction;}
  \item{v}{Nonexceedance probability \eqn{v} in the \eqn{Y} direction;}
  \item{para}{A positive integer \eqn{n \in 1, 2, \cdots}; and}
  \item{...}{Additional arguments to pass.}
}
\value{
  Value(s) for the copula are returned.
}
\references{
Nelsen, R.B., 2006, An introduction to copulas: New York, Springer, 269 p.
}
\author{W.H. Asquith}
\seealso{\code{\link{W}}, \code{\link{Mshuf}}}
\examples{
Wshuf(0.4, 0.6, para=5)

\dontrun{
  UV <- simCOP(1000, cop=Wshuf, para=4) # Nelsen (2006, exer. 5.12a, p. 172, fig. 5.3a)}
}
\keyword{copula (formulas)}
\keyword{copula (limits/bounds)}
\concept{Shuffle of Frechet lower-bound copula}
\concept{Shuffle of Frechet-Hoeffding lower-bound copula}
\concept{Shuffle of Frechet--Hoeffding lower-bound copula}
\concept{Shuffle of Frechet lower bound}
\concept{Shuffle of Frechet-Hoeffding lower-bound}
\concept{Shuffle of Frechet--Hoeffding lower-bound}
\keyword{copula}
\keyword{shuffle}
\keyword{Nelsen (2006) Examples and Exercises}

