\encoding{utf8}
\name{MOcop}
\alias{MOcop}
\title{The Marshall--Olkin Copula}
\description{
The \emph{Marshall--Olkin copula} (Dobrowolski and Kumar, 2014, eq. 2.6; Joe, 2014, p. 183) having parameters \eqn{\alpha, \beta \in (0,1)} for \eqn{u^\alpha \ge v^\beta} is

\deqn{\mathbf{C}_{(\alpha,\beta)}(u,v) = \mathbf{MO}(u,v) = u^{1-\alpha}v\mbox{,}}

and for \eqn{u^\alpha < v^\beta} is

\deqn{\mathbf{C}_{(\alpha,\beta)}(u,v) = \mathbf{MO}(u,v) = uv^{1-\beta}\mbox{.}}

As \eqn{\alpha, \beta \rightarrow 1}, the copula limits to the \emph{comonotonicity coupla} (\eqn{\mathbf{M}(u,v)}; \code{\link{M}}), and as \eqn{\alpha, \beta \rightarrow 0}, the copula limits to the \emph{independence copula} (\eqn{\mathbf{\Pi}(u,v)}; \code{\link{P}}). The copula can be expected to have a visible singularity as the parameters increase and (or) a simulation size becomes large. The copula with \eqn{\alpha = \beta} is \emph{permutation symmetric} (\code{\link{isCOP.permsym}}, \code{\link{breveCOP}}) and is known as the \emph{\enc{Cuadras--Augé}{Cuadras--Auge}} copula. \emph{Kendall tau} (\code{\link{tauCOP}}) of the copula is

\deqn{\tau_\mathbf{C} = \alpha\beta / (\alpha + \beta - \alpha\beta)\mbox{,\ }}

and \emph{Spearman Rho} (\code{\link{rhoCOP}}) is

\deqn{\rho_\mathbf{C} = 3\alpha\beta / (2\alpha + 2\beta - \alpha\beta\mbox{.}}
}
\usage{
MOcop(u, v, para=NULL, ...)
}
\arguments{
  \item{u}{Nonexceedance probability \eqn{u} in the \eqn{X} direction;}
  \item{v}{Nonexceedance probability \eqn{v} in the \eqn{Y} direction;}
  \item{para}{A vector (single element) of parameters---the \eqn{\alpha, \beta} parameters of the copula; and}
  \item{...}{Additional arguments to pass.}
}
\value{
  Value(s) for the copula are returned.
}
\references{
Dobrowolski, E., and Kumar, P., 2014, Some properties of the Marshall--Olkin and generalized Cuadras--Augé familes of copula: Australian Journal of Mathematical Analysis and Applications: v. 11, no. 1, art. 2, pp. 1--13, accessed on August 10, 2025, at \url{https://ajmaa.org/searchroot/files/pdf/v11n1/v11i1p2.pdf}.

Joe, H., 2014, Dependence modeling with copulas: Boca Raton, CRC Press, 462 p.
}
\author{W.H. Asquith}
\seealso{\code{\link{P}}, \code{\link{taildepCOP}}
}
\examples{
UV <- simCOP(500, cop=MOcop, para=c(0.3, 0.75))

\dontrun{
  # This example developed as a demonstrator of the dissarray that the
  # Spearman Footrule measures as a nonunique relation to Spearman Rho.
  # The permutation asymmetry causes the relation between the two to
  # "paint" for positive association as shown in the graphic produced.
  cex <- c(1.3, 1.0, 0.7, 0.4)
  col <- c("#7D0112", "#67C7B1", "#1F28A2", "#CCC86C")
  pch <- c(1, 16, 1, 16)
  plot(c(-0.5, 1), c(0, 1), type="n", xlab="Spearman Footrule", ylab="")
  mtext(paste0("Absolute value of Spearman rho\n",
               "Actual sign is the same as for the Footrule"), side=2, line=2)
  mtext(paste0("Attainable Marshall-Olkin copula, MOcop() with permutation\n",
               "inclusion (breveCOP) and reflections (COP)"), line=1)
  for(r in 1:4) { # reflections, COP
    for(k in seq(-1, 1, by=0.1)) { # permutation asymmetry insertion, breveCOP
      for(a in seq(0, 1, by=0.11)) {   # increment slightly different from beta
        for(b in seq(0, 1, by=0.09)) { # increment slightly different from alpha
          para <- list(cop=MOcop, para=c(a,b), breve=k, reflect=r)
          fotT <- footCOP(cop=breveCOP, para=para)
          rhoT <-  rhoCOP(cop=breveCOP, para=para)
          points(fotT, abs(rhoT), lwd=0.8, pch=pch[r], cex=cex[r], col=col[r])
          #UV <- simCOP(100, cop=breveCOP, para=para)
        }
      }
    }
  } #}
}
\keyword{Marshall-Olkin copula}
\keyword{copula (formulas)}
\keyword{copula}
\concept{Marshall-Olkin copula}
