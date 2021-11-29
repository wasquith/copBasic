\encoding{utf8}
\name{blomatrixCOP}
\alias{blomatrixCOP}
\alias{blomatrixCOPdec}
\alias{blomatrixCOPiqr}
\title{A Matrix of Blomqvist-like Betas of a Copula}
\description{
Compute the \emph{Blomqvist-like Betas matrix} \eqn{\beta^\circ_\mathbf{C}}-matrix of a copula, which is defined at presumably strategic points within \eqn{\mathcal{I}^2} as (for \code{as.blomCOPss=FALSE} argument)

\deqn{\beta^\circ_\mathbf{C} = \frac{\mathbf{C}(u^\circ,v^\circ)}{\mathbf{\Pi}(1/2, 1/2)} - 1\mbox{,}}

where the \eqn{u^\circ} and \eqn{v^\circ} are of two types of gridded locations in \eqn{\mathcal{I}^2} space and if \eqn{u^\circ = 1/2} and \eqn{v^\circ = 1/2}, then central location of the matrix is \emph{Blomqvist Beta} (\code{\link{blomCOP}}). The definition of \eqn{\beta^\circ_\mathbf{C}} is such that the matrix is entirely zero for the \emph{independence copula} (\eqn{\mathbf{\Pi}(u,v)}) (\code{\link{P}}) when \eqn{\mathbf{C}(u^\circ,v^\circ) = \mathbf{\Pi}(u,v)} at the medial location \eqn{u,v=1/2}. Also, the definition here might be unique to the \pkg{copBasic} package. The decile version (\code{blomatrixCOPdec}) of this function uses \eqn{u^\circ \in (1,5,9)/10} and \eqn{v^\circ \in (1,5,9)/10}. Whereas, the quartile version (\code{blomatrixCOPiqr}) of this function uses  \eqn{u^\circ \in (25,50,75)/100} and \eqn{v^\circ \in (25,50,75)/100}. If \code{as.blomCOPss=TRUE} argument is set (\bold{default operation}), then the coordinate locations in the matrix become the \eqn{\beta^\diamond_\mathbf{C}} of \code{\link{blomCOPss}}. As a rule \eqn{\beta^\circ_\mathbf{C} \ne \beta^\diamond_\mathbf{C}}.
}
\usage{
blomatrixCOPdec(cop=NULL, para=NULL, as.sample=FALSE, as.blomCOPss=TRUE,
                  ctype=c("weibull", "hazen", "1/n",
                          "bernstein", "checkerboard"), ...)
blomatrixCOPiqr(cop=NULL, para=NULL, as.sample=FALSE, as.blomCOPss=TRUE,
                  ctype=c("weibull", "hazen", "1/n",
                          "bernstein", "checkerboard"), ...)
}
\arguments{
  \item{cop}{A copula function;}
  \item{para}{Vector of parameters or other data structure, if needed, to pass to the copula;}
  \item{as.sample}{A logical controlling whether an optional \R \code{data.frame} in \code{para} is used to compute the  \eqn{\hat\beta^\circ_\mathbf{C}}-matrix at which point the \code{ctype} argument will be passed to multiple calls of \code{\link{EMPIRcop}};}
  \item{as.blomCOPss}{A logical to trigger \code{\link{blomCOPss}} for each of the \eqn{(u,v)} locations where for the \code{\link{blomCOPss}} calls (\eqn{\beta^\diamond_\mathbf{C}(\mathbf{u}, \mathbf{v})}): \eqn{u \mapsto (u,u) \mapsto \mathbf{u}} and  \eqn{v  \mapsto (v,v) \mapsto \mathbf{v}};}
  \item{ctype}{Argument of the same as \code{\link{EMPIRcop}}; and}
  \item{...}{Additional arguments to pass to the copula.}
}
\value{
  The matrix for \eqn{\beta^\circ_\mathbf{C}} is returned depending on whether the decile or quartile version has been called.
}
\author{W.H. Asquith}
\seealso{\code{\link{blomCOP}}, \code{\link{blomCOPss}}
}
\examples{
round(blomatrixCOPdec(cop=P), digits=8);     round(blomatrixCOPiqr(cop=P), digits=8)
#          U|V=0.10 U|V=0.50 U|V=0.90        #          U|V=0.25 U|V=0.50 U|V=0.75
# U|V=0.90        0        0        0        # U|V=0.75        0        0        0
# U|V=0.50        0        0        0        # U|V=0.50        0        0        0
# U|V=0.10        0        0        0        # U|V=0.25        0        0        0

round(blomatrixCOPdec(cop=PSP, as.blomCOPss=TRUE ), digits=8)
#           U|V=0.10   U|V=0.50   U|V=0.90
# U|V=0.90 0.4736842  0.8181818  0.5153268
# U|V=0.50 0.8181818  0.3333333  0.6459330
# U|V=0.10 0.8901099  0.4736842  0.1708292

round(blomatrixCOPdec(cop=PSP, as.blomCOPss=FALSE), digits=8)
#            U|V=0.10   U|V=0.50  U|V=0.90
# U|V=0.90 0.09890110 0.81818182 4.2631579
# U|V=0.50 0.05263158 0.33333333 0.8181818
# U|V=0.10 0.01010101 0.05263158 0.0989011

\dontrun{
set.seed(1)
td <- c(0.10, 0.50, 0.90)
UVn <- simCOP(n=5000, cop=glueCOP, col=8,
          para=list(glue=0.4, cop1 =PLcop,           cop2=PLcop,
                             para1=PLpar(rho=-0.5), para2=PLpar(rho=+0.5)))
points(td, rep(td[3], 3), cex=2, lwd=2, pch=3, col="red")
points(td, rep(td[2], 3), cex=2, lwd=2, pch=3, col="red")
points(td, rep(td[1], 3), cex=2, lwd=2, pch=3, col="red")

print(blomatrixCOPdec(as.sample=TRUE, para=UVn, ctype="weibull"))
#             U|V=0.10 U|V=0.50 U|V=0.90
# U|V=0.90 -0.08222222   -0.580   -0.190
# U|V=0.50  0.30800000    0.112    0.264
# U|V=0.10  0.84000000    0.620    0.262

BMdn <- blomatrixCOPdec(cop=glueCOP,
          para=list(glue=0.4,  cop1=PLcop,            cop2=PLcop,
                              para1=PLpar(rho=-0.5), para2=PLpar(rho=+0.5)))
print(round(BMdn, digits=8))
#             U|V=0.10   U|V=0.50   U|V=0.90
# U|V=0.90 -0.08110464 -0.5569028 -0.2053772
# U|V=0.50  0.33449815  0.1202744  0.2424668
# U|V=0.10  0.75217766  0.6013719  0.2424668 }
}
\keyword{copula (characteristics)}
\keyword{copula (properties)}
\keyword{Blomqvist Beta}
