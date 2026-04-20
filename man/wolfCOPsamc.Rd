\encoding{utf8}
\name{wolfCOPsamc}
\alias{wolfCOPsamc}
\title{Monte Carlo Integration of Schweizer--Wolff Sigma}
\description{
Perform \emph{Monte Carlo integration} for the measure of association known as \emph{Schweizer--Wolff Sigma} \eqn{\sigma_\mathbf{C} \in (0,1)} of an \emph{empirical copula} (\eqn{\mathbf{C}_n(u,v)}). This function can be useful as the \code{as.sample} computation provided by \code{\link{wolfCOP}} scales very poorly as sample size increases. Please refer to \code{\link{wolfCOP}} for more details, including theory, about \eqn{\sigma_\mathbf{C}}. Sample estimation can be made by

\deqn{\hat\sigma_\mathbf{C} = 12\int\!\!\int_{\mathcal{I}^2} \bigl|\mathbf{C}_n(u,v) - uv\bigr|\,\mathrm{d}u\mathrm{d}v\mbox{,}}

for which the integrals in the above formula are replaced with Monte Carlo integration. Allowance is made internally to swap out \eqn{\mathbf{C}_n(u,v)} for a copula family \eqn{\mathbf{C}(u,v)}, which means this function can be used for integration for families that might have convergence problems in \code{\link{wolfCOP}}.
}
\usage{
wolfCOPsamc(u=NULL, v=NULL, cop=EMPIRcop, para=NULL, para_has_paras=FALSE,
            tol=0.005, minit=10, maxit=100, subdivisions=200, large_n=5E3,
            trace=FALSE, verbose=FALSE, ...)
}
\arguments{
  \item{u}{Optional nonexceedance probability \eqn{u} in the \eqn{X} direction. This argument can also by a two-column data frame from which the \code{para} argument to \code{\link{EMPIRcop}} (if \code{cop=EMPIRcop}) is internally constructed;}
  \item{v}{Optional nonexceedance probability \eqn{v} in the \eqn{Y} direction from which the \code{para} argument to \code{\link{EMPIRcop}} (if \code{cop=EMPIRcop}) is internally constructed;}
  \item{cop}{The empirical copula, \code{\link{EMPIRcop}}, or a copula family;}
  \item{para}{The optional parameter content for the \code{para} argument to the specified copula in \code{cop}. For an empirical copula application (default), the \code{para} must be a two-column data frame, but if \code{cop} is set to a copula family, then this function can be used as an alternative to Monte Carlo integration in lieu of numerical integration in \code{\link{wolfCOP}} (refer to \bold{Examples} how to use the function in this way);}
  \item{para_has_paras}{A logical to inform processing whether the \code{para} argument contains parameters for a copula family in lieu of \code{para} containing the sample for the empirical copula (\code{cop=EMPIRcop}). This function supports the alternative use of the \code{wolfCOPsamc} function for Monte Carlo integration in lieu of numerical integration in \code{\link{wolfCOP}} (refer to \bold{Examples} how to use the function in this way). Internally, \code{para_has_paras} is used to disable random subsetting by setting the sample size variable internally to \code{large_n}. Technically, it is possible internally to detection whether \code{para_has_paras} needed to be true (even if false) as the number of random variates for the \code{halton()} function of the \pkg{randtoolbox} poackage is computed, and the code will silently proceed as if \code{para_has_paras=TRUE};}
  \item{tol}{A tolerance in absolute percent change in the global mean as the iterations progress. This tolerance should be good enough for many practical large sample applications, but can be reduced considerable if the function is used as the backdoor to integration of a given copula family instead of a sample estimator;}
  \item{minit}{Minimum number of iterations regardless of \code{tol} setting;}
  \item{maxit}{Maximum number of iterations;}
  \item{subdivisions}{Number of subdivisions and somewhat an analog to the term used in numerical integration. This subdivisions will be the size of the Monte Carlo sample for each iteration;}
  \item{large_n}{For an empirical copula, as the sample size gets large, not only does the \eqn{\sigma_\mathbf{C}} computation slow down because a double summation is involved, the counting that that empirical copula must also perform gets large, therefore, the \code{para} ultimately fed to the empirical copula itself is randomly subsetting by \code{large_n} to speed up iteration time (at expense of increasing time that \code{tol} presumably is reached);}
  \item{trace}{A logical, if true, that will produce a scatter plot of the sample (or the current iteration subset [refer to \code{large_n}]) along with the random sample points used for the current Monte Carlo iteration;}
  \item{verbose}{A logical providing percent change results by iteration to the console; and}
  \item{...}{Additional arguments to pass to the \code{halton()} function of the \pkg{randtoolbox} and the copula specified by \code{cop} (\code{\link{EMPIRcop}}) or a copula family.}
}
\value{
  A list of the Monte Carlo integration results is returned containing
\item{wolves}{A vector of the \eqn{\sigma_\mathbf{C}} by iteration;}
\item{estimates}{The minimum, mean (the formal estimate \eqn{\hat\sigma_\mathbf{C}}), maximum, standard deviation (computed from second L-moment), L-scale (second L-moment), L-skew, and L-kurtosis of the \code{wolves};}
\item{its}{Number of iterations before termination by \code{maxits} or upon reaching \code{tol};}
\item{lastpctchg}{The percent change on the last iteration.}
}
\note{
Beyond about 1,000 sample sizes \code{wolfCOP(as.sample=TRUE)} scales in time poorly, so we see in the following recipe the giant increase in speed at a small cost in precision on the Schweizer--Wolff Sigma estimate.
\preformatted{
  UV <- simCOP(3000, cop=CLcop, para=3, graphics=FALSE, seed=473)
  theowolf <- wolfCOP( cop=CLcop, para=3) # 0.7864391
  system.time(  smwolf <- wolfCOP(    para=UV, as.sample=TRUE) )
  system.time(  mcwolf <- wolfCOPsamc(para=UV,   usetime=TRUE) )
  #   user  system elapsed   WOLF(sample) = 0.8002883
  # 146.273  13.643 159.894
  #   user  system elapsed   WOLF(mean)   = 0.799871285
  #  14.893   1.225  16.118
}


If one wanted to make some comparisons of computed results between \code{\link{wolfCOP}} and \code{wolfCOPsamc}, the following recipe, with a suitably large sample size, that we see a general spread around the equal value line. The use of Monte Carlo integration approach is not advised when computational speed of \code{\link{wolfCOP}} is sufficient.
\preformatted{
  nsam <- 1000; nsim <- 20; wolf <- mcwolf <- rep(NA, length=nsim)
  for(i in seq_len(nsim)) {
    uv <- as.data.frame(matrix(runif(nsam*2), ncol=2))
      wolf[i] <- wolfCOP(    para=uv, as.sample=TRUE)
    mcwolf[i] <- wolfCOPsamc(para=uv)$estimates[2]
    print(c(i, wolf[i], mcwolf[i]))
  }
  df <- data.frame(sigma=wolf, sigma_mc=mcwolf)
  df <- df[complete.cases(df),]
  xylim <- range( c(range(df$sigma, na.rm=TRUE), range(df$sigma_mc)) )
  plot(df$sigma, df$sigma_mc, xlim=xylim, ylim=xylim,
       xlab="Schweizer-Wolff sigma by whole-sample counting",
       ylab="Schweizer-Wolff sigma by Monte Carlo integration")
  lines(par()$usr[1:2], par()$usr[3:4]) # equal value line }
}
\references{
Nelsen, R.B., 2006, An introduction to copulas: New York, Springer, 269 p.
}
\author{W.H. Asquith}
\seealso{\code{\link{wolfCOP}}, \code{\link{EMPIRcop}}
}
\examples{
n <- 40; para=c(1.3, 2.4)              # Small sample, for speed.
UV <- simCOP(40, para=para, cop=GHcop, seed=1) # Simulate a sample
wolfCOP(    para=para, cop=GHcop)      # Theo. value + next 2 lines
wolfCOPsamc(para=para, cop=GHcop)$estimate[2] # Same as next line
wolfCOPsamc(para=para, cop=GHcop, para_has_paras=TRUE)$estimate[2]
wolfCOP(    para=UV, as.sample=TRUE) # original package implementation
wolfCOPsamc(para=UV)$estimate[2] # Monte Carlo estimate

\dontrun{
  # Test the trace and verbose arguments
  para <- list(cop=c(P, M), part=c(0,0.5,1))
  UV   <- simCOP(1000, cop=ORDSUWcop, para=para)
  SW   <- wolfCOPsamc(para=UV, tol=0.0005, trace=TRUE, verbose=TRUE)
}

\dontrun{
  # Test the various interfaces for a sample computation along with
  # setting the usetime to the random number subsystem.
  UV <- simCOP(1000, cop=CLcop, para=3, graphics=FALSE)
  mc1 <- wolfCOPsamc(para=UV, usetime=TRUE)$estimates[2]
  mc2 <- wolfCOPsamc(u=UV,    usetime=TRUE)$estimates[2]
  mc3 <- wolfCOPsamc(u=UV[,1], v=UV[,2], usetime=TRUE)$estimates[2] #}

\dontrun{
  nsam <- 2000; nsim <- 400; mcwolves <- rep(NA, length=nsim)
  x <- y <- NULL
  lmomco::plotlmrdia(lmomco::lmrdia(), xlim=c(0,0.4), ylim=c(0,0.3),
             autolegend=TRUE, xleg="topleft", ncol=2, nolimits=TRUE,
             nogov=TRUE, nogpa=TRUE, noexp=TRUE, nonor=TRUE,
             noray=TRUE, nouni=TRUE, noaep4=TRUE)
  for(i in seq_len(nsim)) {
    uv <- as.data.frame(matrix(runif(nsam*2), ncol=2))
    mcwolves[i] <- wolfCOPsamc(para=uv)$estimates[2]
    print(c(i, mcwolves[i]))
    lmr <- lmomco::lmoms(mcwolves[! is.na(mcwolves)], no.stop=TRUE)
    if(is.null(lmr) | i < 20) next
    x1 <- lmr$ratios[3]; y1 <- lmr$ratios[4]
    if((par()$usr[1] <= x1 & x1 <= par()$usr[2]) &
       (par()$usr[3] <= y1 & y1 <= par()$usr[4])) {
      if(! is.null(x)) segments(x,y,x1=x1,y1=y1, lwd=0.5, col="grey50")
      points(x1, y1, pch=16, cex=0.9, col="grey50")
    }
   x <- x1; y <- y1
  }
  points(x1, y1, pch=21, col="black", bg="white", lwd=1.3, cex=1.2) # }
}
\keyword{copula (characteristics)}
\keyword{copula (properties)}
\keyword{Monte Carlo integration}
\keyword{Schweizer and Wolff Sigma}
\keyword{Package randtoolbox}
