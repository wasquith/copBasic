\encoding{utf8}
\name{EMPIRgrid_fast}
\alias{EMPIRgrid_fast}
\title{Fast Grid of the Bivariate Empirical Copula}
\description{
Generate a fast gridded representation of the \emph{bivariate empirical copula} (\enc{Hernández-Maldonado}{Hernandez-Maldonado} \emph{et al.}, 2024). Mathematics for the \code{"1/n"} and \code{"bernstein"} copula are shown under \code{\link{EMPIRcop}}. This function is a much more recent addition to the \pkg{copBasic} than its historical ancestor of \code{\link{EMPIRgrid}} because of publication of the algorithm by \enc{Hernández-Maldonado}{Hernandez-Maldonado} \emph{et al.} (2024).
}
\usage{
EMPIRgrid_fast(uv=NULL, ctype=c("1/n", "bernstein"), para=NULL,
               gridonly=FALSE, verbose=FALSE, ...)
}
\arguments{
  \item{uv}{A data frame or matrix of \eqn{u} and \eqn{v} vectors of nonexceedance probabilities (see also \code{para});}
  \item{ctype}{A means to trigger the definition of \eqn{\mathbf{C}_n(u,v)} or \eqn{\mathbf{C}^\mathcal{B}_n(u,v)};}
  \item{para}{A matrix of \eqn{u} and \eqn{v} vectors as an alternative to \code{uv} for syntax compatibility for other situations in \pkg{copBasic};}
  \item{verbose}{A logical controlling whether the progress during grid building is to be shown;}
  \item{gridonly}{A logical to return the \code{empcop} component (\eqn{\mathbf{C}_{n}(u,v)} or \eqn{\mathbf{C}^\mathcal{B}_n(u,v)}); and}
  \item{...}{Additional arguments to pass to \code{\link{EMPIRcop}}.}
}
\value{
  An \R \code{list} of the gridded values of \eqn{u}, \eqn{v}, and \eqn{\mathbf{C}_{n}(u,v)} values of the bivariate empirical copula is returned. (Well only \eqn{\mathbf{C}_{n}(u,v)} is in the form of a grid as an \R \code{matrix}.) The \code{deluv} used to generated the grid also is returned. Note, there is no \code{deluv} argument to this function as there is for the historical \code{\link{EMPIRgrid}} of the package. The \code{deluv} is set as \eqn{1/n} (inversion of sample size).
}
\references{
\enc{Hernández-Maldonado}{Hernandez-Maldonado}, V.M., Erdely, A., \enc{Díaz-Viera}{Diaz-Viera}, M., and Rios, L., 2024, Fast procedure to compute empirical and Bernstein copulas: Applied Mathematics and Computation, v. 477, article 128827, 14 p., \doi{10.1016/j.amc.2024.128827}.
}
\author{W.H. Asquith}
\seealso{\code{\link{EMPIRgrid}}, \code{\link{EMPIRgrid_lkup}}}
\examples{
  n <- 5; uv <- as.matrix(matrix(runif(n *2), ncol=2))
  empCOP <- EMPIRgrid_fast(uv, ctype="1/n",       gridonly=TRUE)
            image(empCOP, main="1/n")
  empCOP <- EMPIRgrid_fast(uv, ctype="bernstein", gridonly=TRUE)
            image(empCOP, main="Bernstein")

\dontrun{
  # This permutation asymmetric copula is shown as a means to cross check
  # the orientation of the matrix computations within the function.
  n <- 500
  para <- list(cop=GHcop, para=c(12,15), breve=0.5)
  uv <- copBasic::simCOP(n, cop=breveCOP, para=para, graphics=FALSE)
  empCOP <- EMPIRgrid_fast(uv, gridonly=TRUE)
  image(empCOP); points(uv) #}
}
\keyword{empirical copula}
\keyword{empirical copula (utility)}
\keyword{copulatic surface}
\concept{Bernstein copula}
\concept{Bernstein empirical copula}
\keyword{Bernstein polynomial}
\keyword{gridded empirical copula}
