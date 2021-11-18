\encoding{utf8}
\name{blomatrixCOP}
\alias{blomatrixCOPdec}
\alias{blomatrixCOPiqr}
\title{A Matrix of Blomqvist-like Betas of a Copula}
\description{
Compute the \emph{Blomqvist-like Betas matrix} \eqn{\beta^\circ_\mathbf{C}}-matrix of a copula, which is defined at presumably strategic points within \eqn{\mathcal{I}^2} as

\deqn{\beta^\circ_\mathbf{C} = \frac{\mathbf{C}(u^\circ,v^\circ)}{\mathbf{\Pi}(1/2, 1/2)} - 1\mbox{,}}

where the \eqn{u^\circ} and \eqn{v^\circ} are of two types of gridded locations in \eqn{\mathcal{I}^2} space and if \eqn{u^\circ = 1/2} and \eqn{v^\circ = 1/2}, then central location of the matrix is \emph{Blomqvist Beta} (\code{\link{blomCOP}}). The definition of \eqn{\beta^\circ_\mathbf{C}} is such that the matrix is entirely zero for the \emph{independence copula} (\eqn{\mathbf{\Pi}(u,v)}) (\code{\link{P}}) when \eqn{\mathbf{C}(u^\circ,v^\circ) = \mathbf{\Pi}(u,v)} at the medial location \eqn{u,v=1/2}. Also, the definition here might be unique to the \pkg{copBasic} package. The decile version (\code{blomatrixCOPdec}) of this function uses \eqn{u^\circ \in (1,5,9)/10} and \eqn{v^\circ \in (1,5,9)/10}. Whereas, the quartile version (\code{blomatrixCOPiqr}) of this function uses  \eqn{u^\circ \in (2,5,7)/10} and \eqn{v^\circ \in (2,5,7)/10}.
}
\usage{
blomatrixCOPdec(cop=NULL, para=NULL, as.sample=FALSE,
                  ctype=c("na", "weibull", "hazen",
                                "1/n", "bernstein", "checkerboard"), ...)
blomatrixCOPiqr(cop=NULL, para=NULL, as.sample=FALSE,
                  ctype=c("na", "weibull", "hazen",
                                "1/n", "bernstein", "checkerboard"), ...)
}
\arguments{
  \item{cop}{A copula function;}
  \item{para}{Vector of parameters or other data structure, if needed, to pass to the copula;}
  \item{as.sample}{A logical controlling whether an optional \R \code{data.frame} in \code{para} is used to compute the \eqn{\hat\beta^\circ_\mathbf{C}}-matrix; and}
\item{as.sample}{A logical controlling whether an optional \R \code{data.frame} in \code{para} is used to compute the  \eqn{\hat\beta^\circ_\mathbf{C}}-matrix. The argument is automatically set to true if \code{ctype != "na"} and then the \code{ctype} triggers either the Weibull sample estimator through the \emph{empirical copula}. The combination \code{as.sample=TRUE} and \cr \code{ctype="na"} will trigger the Weibull sample estimator as if \code{ctype="weibull"};}
  \item{ctype}{Argument of the same as \code{\link{EMPIRcop}} with the exception of the \code{"na"} support; and}
  \item{...}{Additional arguments to pass to the copula.}
}
\value{
  The matrix for \eqn{\beta^\circ_\mathbf{C}} is returned depending on whether the decile or quartile version has been called.
}
\author{W.H. Asquith}
\seealso{\code{\link{blomCOP}}
}
\examples{
round(blomatrixCOPdec(cop=P), digits=8);     round(blomatrixCOPiqr(cop=P), digits=8)
#          U|V=0.10 U|V=0.50 U|V=0.90        #          U|V=0.25 U|V=0.50 U|V=0.75
# U|V=0.90        0        0        0        # U|V=0.75        0        0        0
# U|V=0.50        0        0        0        # U|V=0.50        0        0        0
# U|V=0.10        0        0        0        # U|V=0.25        0        0        0

round(blomatrixCOPdec(cop=PSP), digits=8);   round(blomatrixCOPiqr(cop=PSP), digits=8)
#            U|V=0.10   U|V=0.50  U|V=0.90   #            U|V=0.25  U|V=0.50  U|V=0.75
# U|V=0.90 0.09890110 0.81818182 4.2631579   # U|V=0.75 0.23076923 0.6000000 1.2857143
# U|V=0.50 0.05263158 0.33333333 0.8181818   # U|V=0.50 0.14285714 0.3333333 0.6000000
# U|V=0.10 0.01010101 0.05263158 0.0989011   # U|V=0.25 0.06666667 0.1428571 0.2307692

\dontrun{
set.seed(1)
td <- c(0.10, 0.50, 0.90)
UVn <- simCOP(n=5000, cop=glueCOP, col=8,
          para=list(glue=0.4, cop1 =PLcop,           cop2=PLcop,
                             para1=PLpar(rho=-0.5), para2=PLpar(rho=+0.5)))
points(td, rep(td[3], 3), cex=2, lwd=2, pch=3, col="red")
points(td, rep(td[2], 3), cex=2, lwd=2, pch=3, col="red")
points(td, rep(td[1], 3), cex=2, lwd=2, pch=3, col="red")

blomatrixCOPdec(para=UVn, ctype="weibull")
#             U|V=0.10    U|V=0.50    U|V=0.90
# U|V=0.90 -0.12888889 -0.58000000 -0.74000000
# U|V=0.50  0.02444444  0.11200000  0.30800000
# U|V=0.10  0.02197531  0.06888889  0.09333333

BMdn <- blomatrixCOPdec(cop=glueCOP,
          para=list(glue=0.4,  cop1=PLcop,            cop2=PLcop,
                              para1=PLpar(rho=-0.5), para2=PLpar(rho=+0.5)))
print(round(BMdn, digits=8))
#             U|V=0.10    U|V=0.50   U|V=0.90
# U|V=0.90 -0.14709225 -0.55690283 -0.7299417
# U|V=0.50  0.01671506  0.12027439  0.3344981
# U|V=0.10  0.02064803  0.06681911  0.0835753 }

\dontrun{

# }
}
\keyword{copula (characteristics)}
\keyword{copula (properties)}
