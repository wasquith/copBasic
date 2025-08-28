\encoding{utf8}
\name{gEVcop}
\alias{gEVcop}
\title{The Gaussian-based (Extreme Value) Copula}
\description{
The \emph{g-EV copula} (Joe, 2014, p. 105) is a limiting form of the \emph{Gaussian copula}:
\deqn{
\mathbf{C}_{\rho}(u,v) = \mathbf{gEV}(u,v; \rho) =
\mathrm{exp}\bigl(-A(x,y; \rho)\bigr)\mbox{,}
}
where \eqn{x = -\log(u)}, \eqn{y = -\log(v)}, and
\deqn{
A(x,y; \rho) = y\mbox{,}
}
for \eqn{0 \le x/(x+y) \le \rho^2/(1+\rho^2)},
\deqn{
A(x,y; \rho) = (x+y - 2\rho\sqrt{xy})/(1-\rho^2)\mbox{,}
}
for \eqn{\rho^2/(1+\rho^2) \le x/(x+y) \le 1/(1+\rho^2)},
\deqn{
A(x,y; \rho) = x\mbox{,}
}
for \eqn{1/(1+\rho^2) \le x/(x+y) \le 1} and where \eqn{\rho \in [0,1]}. A somewhat curious observation is that this copula has relative hard boundaries into the upper-left and lower-right corners when compared to the other copulas supported by the \pkg{copBasic} package. In other words, the hull defined by the copula has a near hard (not fuzzy) curvilinear boundaries that adjust with the parameter \eqn{\rho}.
}
\usage{
gEVcop(u, v, para=NULL, ...)
}
\arguments{
  \item{u}{Nonexceedance probability \eqn{u} in the \eqn{X} direction;}
  \item{v}{Nonexceedance probability \eqn{v} in the \eqn{Y} direction;}
  \item{para}{The parameter \eqn{\rho}; and}
  \item{...}{Additional arguments to pass.}
}
\value{
  Value(s) for the copula are returned.
}
\references{
Joe, H., 2014, Dependence modeling with copulas: Boca Raton, CRC Press, 462 p.
}
\author{W.H. Asquith}
\seealso{\code{\link{tEVcop}}
}
\examples{
\dontrun{
  UV <- simCOP(200, cop=gEVcop, para=0.8) #}

\dontrun{
  # Joe (2014, p. 105) has brief detail indicating rho = [0,1] and though it seems
  # Rho would be a Pearson correlation, this does not seem to be the case. The Rho
  # seems to start with that of the Gaussian and then through the extreme-value
  # transform, it just assumes the role of a parameter called Rho.
  rho <- 0.8
  UV <- simCOP(2000, cop=gEVcop, para=rho)
  P <- cor(UV[,1], UV[,2], method="pearson")
  if(abs(P - rho) < 0.001) {
    print("Yes same")
  } else { print("nope not") } # }

\dontrun{
  gEVparameter <- seq(0.02, 1, by=0.02)
  SpearmanRho <- sapply(gEVparameter, function(k) rhoCOP(cop=gEVcop, para=k))
  plot(gEVparameter, SpearmanRho)

  for(rho in gEVparameter) UV <- simCOP(n=1000, cop=gEVcop, para=rho) #}

\dontrun{
  # Demonstration of L-comoment ratio diagram as a function of Spearman Rho
  # for the gEVcop(), GHcop(), GLcop(), and HRcop() where each copula is
  # fit to the Spearman Rho and the L-comoments computed. The copulas here
  # are permutation symmetric, which means the T*.12 and T*.21 are identical,
  # but plotting for the T*.21 are shown here to make it easier for the user
  # to insert permutation asymmetry say by blending in the breveCOP() or other
  # composition frameworks. L-cokurtosis is negated to have it simply plot
  # away from the L-coskew.
  gefunc <- function(k, rho=NA) rho - rhoCOP(cop=gEVcop,para=pnorm(k))
  ghfunc <- function(k, rho=NA) rho - rhoCOP(cop=GHcop, para=exp(k)+1)
  glfunc <- function(k, rho=NA) rho - rhoCOP(cop=GLcop, para=exp(k)  )
  hrfunc <- function(k, rho=NA) rho - rhoCOP(cop=HRcop, para=exp(k)  )
  D <- NULL; rs <- c(0.002, 0.005, seq(0.04, 0.98, by=0.02), 0.995, 0.998)
  for(r in rs) { # Spearman Rhos
    ge   <- pnorm(uniroot(gefunc, c(-4,+4), rho=r)$root)
    GE   <- lcomCOP(cop=gEVcop, para=ge)
    gh   <- exp(  uniroot(ghfunc, c(-9,4), rho=r)$root)+1
    GH   <- lcomCOP(cop=GHcop, para=gh)
    gl   <- exp(  uniroot(glfunc, c(-9,4), rho=r)$root)
    GL   <- lcomCOP(cop=GLcop, para=gl)
    hr   <- exp(  uniroot(hrfunc, c(-9,4), rho=r)$root)
    HR   <- lcomCOP(cop=HRcop, para=hr)
    x <- data.frame(rho=r,
                    GEpara=ge, GHpara=gh, GLpara=gl, HRpara=hr,
                    GET2.12=GE$lcomUV[2], GET2.21=GE$lcomVU[2],
                    GET3.12=GE$lcomUV[3], GET3.21=GE$lcomVU[3],
                    GET4.12=GE$lcomUV[4], GET4.21=GE$lcomVU[4],
                    GHT2.12=GH$lcomUV[2], GHT2.21=GH$lcomVU[2],
                    GHT3.12=GH$lcomUV[3], GHT3.21=GH$lcomVU[3],
                    GHT4.12=GH$lcomUV[4], GHT4.21=GH$lcomVU[4],
                    GLT2.12=GL$lcomUV[2], GLT2.21=GL$lcomVU[2],
                    GLT3.12=GL$lcomUV[3], GLT3.21=GL$lcomVU[3],
                    GLT4.12=GL$lcomUV[4], GLT4.21=GL$lcomVU[4],
                    HRT2.12=HR$lcomUV[2], HRT2.21=HR$lcomVU[2],
                    HRT3.12=HR$lcomUV[3], HRT3.21=HR$lcomVU[3],
                    HRT4.12=HR$lcomUV[4], HRT4.21=HR$lcomVU[4] )
    D <- rbind(D, x)
  }
  row.names(D) <- NULL

  cols <- c("darkorange3", "seagreen3", "darkorchid3", "skyblue3")
  plot(rs, D$GHT3.12, ylim=c(-0.08,0.08), type="n",
       xlab="SPEARMAN RHO", ylab="L-COSKEW OR NEGATED L-COKURTOSIS")
  lines(rs, D$GET3.12, col=cols[1]); lines(rs, -D$GET4.12, col=cols[1], lty=2)
  lines(rs, D$GET3.21, col=cols[1]); lines(rs, -D$GET4.21, col=cols[1], lty=2)
  lines(rs, D$GHT3.12, col=cols[2]); lines(rs, -D$GHT4.12, col=cols[2], lty=2)
  lines(rs, D$GHT3.21, col=cols[2]); lines(rs, -D$GHT4.21, col=cols[2], lty=2)
  lines(rs, D$GLT3.12, col=cols[3]); lines(rs, -D$GLT4.12, col=cols[3], lty=2)
  lines(rs, D$GLT3.21, col=cols[3]); lines(rs, -D$GLT4.21, col=cols[3], lty=2)
  lines(rs, D$HRT3.12, col=cols[4]); lines(rs, -D$HRT4.12, col=cols[4], lty=2)
  lines(rs, D$HRT3.21, col=cols[4]); lines(rs, -D$HRT4.21, col=cols[4], lty=2)
  txt <- c("L-coskew gEVcop()", "Negated L-cokurtosis gEVcop()",
           "L-coskew GHcop()",  "Negated L-cokurtosis GHcop()",
           "L-coskew GLcop()",  "Negated L-cokurtosis GLcop()",
           "L-coskew HRcop()",  "Negated L-cokurtosis HRcop()")
  legend(0.3, 0.02, txt, col=c(cols[c(1,1,2,2,3,3,4,4)]), lty=c(1,2)) # }
}
\keyword{g-EV copula}
\keyword{copula (formulas)}
\keyword{copula}
\keyword{package copula (comparison to)}
\concept{g-EV copula}
\keyword{copula (extreme value)}
\keyword{extreme value copula}
\keyword{Spearman Rho}
\keyword{Pearson correlation coefficient}
\keyword{L-comoments (ratio diagram)}
