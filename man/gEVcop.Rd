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
  clfunc <- function(k, rho=NA) rho - rhoCOP(cop=COP,
                           para=list(reflect=2, cop=CLcop,  para=exp(k)-1))
  gefunc <- function(k, rho=NA) rho - rhoCOP(cop=gEVcop,    para=pnorm(k))
  ghfunc <- function(k, rho=NA) rho - rhoCOP(cop=GHcop,     para=exp(k)+1)
  glfunc <- function(k, rho=NA) rho - rhoCOP(cop=GLcop,     para=exp(k)  )
  hrfunc <- function(k, rho=NA) rho - rhoCOP(cop=HRcop,     para=exp(k)  )
  nnfunc <- function(k, rho=NA) rho - rhoCOP(cop=rN4220cop, para=exp(k)  )
  D <- NULL; rs <- c(0.002, 0.005, seq(0.04, 0.98, by=0.02), 0.995, 0.998)
  for(r in rs) { # Spearman Rhos
    cl   <- exp(  uniroot(clfunc, c(-9,+5), rho=r)$root)-1
    CL   <- lcomCOP(cop=COP, para=list(cop=CLcop, para=cl, reflect=2))
    fr   <- FRcop(rhotau=r, cortype="rho")$para
    FR   <- lcomCOP(cop=FRcop, para=fr)
    ge   <- pnorm(uniroot(gefunc, c(-4,+4), rho=r)$root)
    GE   <- lcomCOP(cop=gEVcop, para=ge)
    gh   <- exp(  uniroot(ghfunc, c(-9,+4), rho=r)$root)+1
    GH   <- lcomCOP(cop=GHcop, para=gh)
    gl   <- exp(  uniroot(glfunc, c(-9,+4), rho=r)$root)
    GL   <- lcomCOP(cop=GLcop, para=gl)
    hr   <- exp(  uniroot(hrfunc, c(-9,+4), rho=r)$root)
    HR   <- lcomCOP(cop=HRcop, para=hr)
    nn   <- exp(  uniroot(nnfunc, c(-9,+4), rho=r)$root)
    NN   <- lcomCOP(cop=rN4220cop, para=nn)
    x <- data.frame(rho=r, CLpara=cl, FRpara=fr, GEpara=ge,
                GHpara=gh, GLpara=gl, HRpara=hr, NNpara=nn,
                 CLT2.12=CL$lcomUV[2],   CLT2.21=CL$lcomVU[2],
                 CLT3.12=CL$lcomUV[3],   CLT3.21=CL$lcomVU[3],
                 CLT4.12=CL$lcomUV[4],   CLT4.21=CL$lcomVU[4],
                 FRT2.12=FR$lcomUV[2],   FRT2.21=FR$lcomVU[2],
                 FRT3.12=FR$lcomUV[3],   FRT3.21=FR$lcomVU[3],
                 FRT4.12=FR$lcomUV[4],   FRT4.21=FR$lcomVU[4],
                 GET2.12=GE$lcomUV[2],   GET2.21=GE$lcomVU[2],
                 GET3.12=GE$lcomUV[3],   GET3.21=GE$lcomVU[3],
                 GET4.12=GE$lcomUV[4],   GET4.21=GE$lcomVU[4],
                 GHT2.12=GH$lcomUV[2],   GHT2.21=GH$lcomVU[2],
                 GHT3.12=GH$lcomUV[3],   GHT3.21=GH$lcomVU[3],
                 GHT4.12=GH$lcomUV[4],   GHT4.21=GH$lcomVU[4],
                 GLT2.12=GL$lcomUV[2],   GLT2.21=GL$lcomVU[2],
                 GLT3.12=GL$lcomUV[3],   GLT3.21=GL$lcomVU[3],
                 GLT4.12=GL$lcomUV[4],   GLT4.21=GL$lcomVU[4],
                 HRT2.12=HR$lcomUV[2],   HRT2.21=HR$lcomVU[2],
                 HRT3.12=HR$lcomUV[3],   HRT3.21=HR$lcomVU[3],
                 HRT4.12=HR$lcomUV[4],   HRT4.21=HR$lcomVU[4],
                 NNT2.12=NN$lcomUV[2],   NNT2.21=NN$lcomVU[2],
                 NNT3.12=NN$lcomUV[3],   NNT3.21=NN$lcomVU[3],
                 NNT4.12=NN$lcomUV[4],   NNT4.21=NN$lcomVU[4] )
    D <- rbind(D, x)
  }
  row.names(D) <- NULL

  cols <- c("turquoise", "firebrick3", "darkorange1",
            "forestgreen", "darkorchid3", "skyblue3", "black")
  plot(rs, D$GHT3.12, ylim=c(-0.15, 0.20), type="n",
       xlab="SPEARMAN RHO", ylab="L-COSKEW (cT3) OR NEGATED L-COKURTOSIS (cT4)")
  lines(rs, D$CLT3.12, col=cols[1]); lines(rs, -D$CLT4.12, col=cols[1], lty=2)
  lines(rs, D$CLT3.21, col=cols[1]); lines(rs, -D$CLT4.21, col=cols[1], lty=2)
  lines(rs, D$FRT3.12, col=cols[2]); lines(rs, -D$FRT4.12, col=cols[2], lty=2)
  lines(rs, D$FRT3.21, col=cols[2]); lines(rs, -D$FRT4.21, col=cols[2], lty=2)
  lines(rs, D$GET3.12, col=cols[3]); lines(rs, -D$GET4.12, col=cols[3], lty=2)
  lines(rs, D$GET3.21, col=cols[3]); lines(rs, -D$GET4.21, col=cols[3], lty=2)
  lines(rs, D$GHT3.12, col=cols[4]); lines(rs, -D$GHT4.12, col=cols[4], lty=2)
  lines(rs, D$GHT3.21, col=cols[4]); lines(rs, -D$GHT4.21, col=cols[4], lty=2)
  lines(rs, D$GLT3.12, col=cols[5]); lines(rs, -D$GLT4.12, col=cols[5], lty=2)
  lines(rs, D$GLT3.21, col=cols[5]); lines(rs, -D$GLT4.21, col=cols[5], lty=2)
  lines(rs, D$HRT3.12, col=cols[6]); lines(rs, -D$HRT4.12, col=cols[6], lty=2)
  lines(rs, D$HRT3.21, col=cols[6]); lines(rs, -D$HRT4.21, col=cols[6], lty=2)
  lines(rs, D$NNT3.12, col=cols[7]); lines(rs, -D$NNT4.12, col=cols[7], lty=2)
  lines(rs, D$NNT3.21, col=cols[7]); lines(rs, -D$NNT4.21, col=cols[7], lty=2)
  txt <- c("cT3 CLcop(reflected)",  "Negated cT4 CLcop(reflected)",
           "cT3 FRcop",  "Negated cT4 FRcop",
           "cT3 gEVcop", "Negated cT4 gEVcop",
           "cT3 GHcop",  "Negated cT4 GHcop",
           "cT3 GLcop",  "Negated cT4 GLcop",
           "cT3 HRcop",  "Negated cT4 HRcop",
           "cT3 rN4220cop", "Negated cT4 rN4220cop")
   legend("bottomleft", txt, col=c(cols[c(1,1,2,2,3,3,4,4,5,5,6,6,7,7)]),
           cex=0.8, ncol=2, lty=c(1,2))
  mtext("L-comoment ratio diagram for various copula", line=0.8)
  # Now, simulation, note that in samples T[2+].[12|21] are not perfectly equal
  # even for known permutation asymmetry. These could be average or used in
  # isolation. Here, we just plot the two separately. So, there are double the
  # data points. It takes a lot of samples for single sample to look authoritative.
  T2 <- T3 <- T4 <- NULL
  for(i in 1:60) {
    UV <- simCOP(1500, cop=GLcop, para=1.4, graphics=FALSE)
    lcm <- lmomco::lcomoms2(UV, nmom=4); srho <- (lcm$T2[1,2] + lcm$T2[2,1]) / 2
    points(rep(srho, 2),  c(lcm$T3[1,2], lcm$T3[2,1]), cex=0.6, lwd=0.8, pch=24)
    points(rep(srho, 2), -c(lcm$T4[1,2], lcm$T4[2,1]), cex=0.6, lwd=0.8, pch=25)
    T2 <- c(T2, srho)
    T3 <- c(T3, mean(c(lcm$T3[1,2], lcm$T3[2,1])))
    T4 <- c(T4, mean(c(lcm$T4[1,2], lcm$T4[2,1])))
  } # BTW : compare the srho to srho <- cor(UV, method="spearman")
  mT2 <- mean(T2); mT3 <- mean(T3); mT4 <- mean(T4)
  points(mT2,  mT3, col=cols[5], lwd=2, bg="darkorchid1", pch=24, cex=2)
  points(mT2, -mT4, col=cols[5], lwd=2, bg="darkorchid1", pch=25, cex=2)
  # ----------------------------------------------------
  # L-comoment ratio diagram of L-coskew and L-cokurtosis.
  # Again both 12 and 21 (UV and VU) are plotted but for permutation
  # symmetry these "with respect to's" have the same L-comoments.
  plot(D$GHT3.12, xlim=c(0, 0.20), ylim=c(-0.05, 0.10), type="n",
       xlab="L-COSKEW (cT3)", ylab="L-COKURTOSIS (cT4)")
  lines(D$CLT3.12, D$CLT4.12, col=cols[1])
  lines(D$CLT3.12, D$CLT4.12, col=cols[1], lty=2)
  lines(D$CLT3.21, D$CLT4.21, col=cols[1])
  lines(D$CLT3.21, D$CLT4.21, col=cols[1], lty=2)
  lines(D$FRT3.12, D$FRT4.12, col=cols[2])
  lines(D$FRT3.12, D$FRT4.12, col=cols[2], lty=2)
  lines(D$FRT3.21, D$FRT4.21, col=cols[2])
  lines(D$FRT3.21, D$FRT4.21, col=cols[2], lty=2)
  lines(D$GET3.12, D$GET4.12, col=cols[3])
  lines(D$GET3.12, D$GET4.12, col=cols[3], lty=2)
  lines(D$GET3.21, D$GET4.21, col=cols[3])
  lines(D$GET3.21, D$GET4.21, col=cols[3], lty=2)
  lines(D$GHT3.12, D$GHT4.12, col=cols[4])
  lines(D$GHT3.12, D$GHT4.12, col=cols[4], lty=2)
  lines(D$GHT3.21, D$GHT4.21, col=cols[4])
  lines(D$GHT3.21, D$GHT4.21, col=cols[4], lty=2)
  lines(D$GLT3.12, D$GLT4.12, col=cols[5])
  lines(D$GLT3.12, D$GLT4.12, col=cols[5], lty=2)
  lines(D$GLT3.21, D$GLT4.21, col=cols[5])
  lines(D$GLT3.21, D$GLT4.21, col=cols[5], lty=2)
  lines(D$HRT3.12, D$HRT4.12, col=cols[6])
  lines(D$HRT3.12, D$HRT4.12, col=cols[6], lty=2)
  lines(D$HRT3.21, D$HRT4.21, col=cols[6])
  lines(D$HRT3.21, D$HRT4.21, col=cols[6], lty=2)
  lines(D$NNT3.12, D$NNT4.12, col=cols[7])
  lines(D$NNT3.12, D$NNT4.12, col=cols[7], lty=2)
  lines(D$NNT3.21, D$NNT4.21, col=cols[7])
  lines(D$NNT3.21, D$NNT4.21, col=cols[7], lty=2)
  points(T3, T4, lwd=0.8, cex=0.8)
  points(mT3, mT4, cex=2, pch=21, col=cols[6], bg="darkorchid1")
  txt <- c(    "CLcop(reflected)[UV]",           "CLcop(reflected)[VU]",
            "FRcop[UV]",  "FRcop[VU]",    "gEVcop[UV]",    "gEVcop[VU]",
            "GHcop[UV]",  "GHcop[VU]",     "GLcop[UV]",     "GLcop[VU]",
            "HRcop[UV]",  "HRcop[VU]", "rN4220cop[UV]", "rN4220cop[VU]")
   legend("bottomright", txt, col=c(cols[c(1,1,2,2,3,3,4,4,5,5,6,6,7,7)]),
           cex=0.8, ncol=2, lty=c(1,2))
   mtext("L-comoment ratio diagram for various copula", line=0.8) # }
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
\keyword{L-comoment ratio diagram}
