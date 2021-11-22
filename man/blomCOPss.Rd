\encoding{utf8}
\name{blomCOPss}
\alias{blomCOPss}
\title{Blomqvist (Schmid--Schmidt) Betas of a Copula}
\description{
Compute the \emph{Blomqvist (Schmid--Schmidt) Betas} \eqn{\beta^\diamond_\mathbf{C}} (Schmid and Schmidt, 2007) defined for arbitrary dimension \eqn{d} of a copula \eqn{\mathbf{C}_(u_1, \cdots, u_d; \Theta)} (\code{\link{COP}}) and survival function of a copula \eqn{\overline\mathbf{C}(u_1, \cdots, u_d; \Theta)} (\code{\link{surfuncCOP}}) for parameters \eqn{\Theta}, though the \pkg{copBasic} package is built around bivariate copula only, as
\deqn{\beta^\diamond_\mathbf{C} = h_d(\mathbf{u}, \mathbf{v})\bigl[
\bigl(\mathbf{C}(\mathbf{u}) + \overline\mathbf{C}(\mathbf{v})\bigr) - g_d(\mathbf{u}, \mathbf{v})
\bigr]\mbox{,}}
where \eqn{h_d} and \eqn{g_d} are norming constants defined below. The superscript \eqn{\diamond} (diamond) is chosen for this package because of the alliteration to \dQuote{dimension}. The bold face for \eqn{\mathbf{u}} and \eqn{\mathbf{v}} are indications that these arguments are vectors of length \eqn{d} reflecting \dQuote{cutting points} on nonexceedance probabilities in each of the dimensions. The \eqn{\mathbf{u}} functions as the arguments \eqn{(u,v)} pair used in copula of this package and represents the first cutting point for a \eqn{\mathrm{Pr}[U \le u, V \le v] = \mathbf{C}(u,v)}, and \eqn{\mathbf{v}} functions as the arguments \eqn{u,v} pair for this package and represents the second cutting point for a \eqn{\mathrm{Pr}[U > u, V > v] = 1 - u - v + \mathbf{C}(u,v) = \overline\mathbf{C}(u,v)}. This notation of vectored (bold face) and nonvectored \dQuote{u} and \dQuote{v} is a little obtuse but as the properties of \eqn{\beta^\diamond_\mathbf{C}} are summarized clarity for the reader is anticipated.

The norming constant \eqn{h_d} is defined as
\deqn{h_d(\mathbf{u}, \mathbf{v}) = \frac{1}{\bigl( \mathrm{min}(u_1,   \cdots,   u_d) + \mathrm{min}(1-v_1, \cdots, 1-v_d) - g_d(\mathbf{u}, \mathbf{v})\bigr)}\mbox{,}}
and \eqn{g_d} is defined as
\deqn{g_d(\mathbf{u}, \mathbf{v}) = \prod^d_{i=1}u_i + \prod^d_{i=1}(1-v_i)\mbox{,}}
where the cutting points \eqn{\mathbf{u}} and \eqn{\mathbf{v}} are in a domain \eqn{D : \{(\mathbf{u}, \mathbf{v})\} \in [0,1]^{2d}} given \eqn{\mathbf{u} \le \mathbf{v}} and \eqn{\mathbf{u} > 0} or \eqn{\mathbf{v} < \mathbf{1}}. The reader must careful remember that these \eqn{\mathbf{u}} and \eqn{\mathbf{v}} are vectors of probabilities. Finally, the norming constants provide for
\eqn{-1 \le \beta^\diamond_\mathbf{C} \le +1}, \eqn{\beta^\diamond_\mathbf{C} = 1} if \eqn{\mathbf{C} = \mathbf{M}} \emph{comonotonicity copula} (\code{\link{M}}) (\code{blomCOPss(cop=M) == 1}), and \eqn{\beta^\diamond_\mathbf{C} = 0} if \eqn{\mathbf{C} = \mathbf{P}} \emph{independence copula} (\code{\link{P}}) (\code{blomCOPss(cop=P)}) using the function argument defaults for \eqn{d=2} dimensions \eqn{\mathbf{u} = (1,1)/2} for \code{uu} and \eqn{\mathbf{v} = (1,1)/2} for \code{vv}.

Schmid and Schmidt (2007) list three important cases extending the \eqn{\mathbf{M}} and \eqn{\mathbf{P}} examples. First, \eqn{\beta^\diamond_\mathbf{C}(\mathbf{1/2}, \mathbf{1/2}) = \beta_\mathbf{C}(1/2, 1/2)}, which is \emph{Blomqvist Beta} (\code{\link{blomCOP}}) and measures overall dependence.

Second, \eqn{\beta^\diamond_\mathbf{C}(\mathbf{u}, \mathbf{v})} with \eqn{\mathbf{u} < 1/2 < \mathbf{v}}, which measures dependence in the tail regions. (Note, the author things \dQuote{regions} in plural is need in the previous sentence; Schmid and Schmidt use the singular. This is potentially important as seemingly simultaneous tail dependency in the lower and upper perspectives would be provided. More discussion is provided in \bold{Examples}.)

Third and presumably very important in practical applications, \eqn{\mathrm{lim}_{p\downarrow 0}\, \beta^\diamond_\mathbf{C}(\mathbf{p}, \mathbf{1}) = \lambda^L_{\beta^\diamond_\mathbf{C}}} for \eqn{\mathbf{u} = \mathbf{p} = (p,\cdots,p)} measures lower-tail dependence. This measure is equal to the \emph{lower-tail dependence parameter} \eqn{\lambda^L_\mathbf{C} = \lambda^L_{\beta^\diamond_\mathbf{C}}} without some of the computational nuances required as defined at \code{\link{taildepCOP}}.

Schmid and Schmidt (2007) do not list how the \emph{upper-tail dependence parameter} \eqn{\lambda^U_\mathbf{C}} could be computed in terms of \eqn{\beta^\diamond}.

\eqn{\lambda^U_{\beta^\diamond_\mathbf{C}}}
}
\usage{
blomCOPss(cop=NULL, para=NULL, as.sample=FALSE,
          uu=rep(0.5, 2), vv=rep(0.5, 2), trap.nan=TRUE, ...)
}
\arguments{
  \item{cop}{A copula function;}
  \item{para}{Vector of parameters or other data structure, if needed, to pass to the copula;}
  \item{as.sample}{A logical controlling whether an optional \R \code{data.frame} in \code{para} is used to compute the \eqn{\hat\beta^\diamond_\mathbf{C}}; and}
\item{as.sample}{A logical controlling whether an optional \R \code{data.frame} in \code{para} is used to compute the  \eqn{\hat\beta^\diamond_\mathbf{C}}. The argument is automatically set to true if \code{ctype != "na"} and then the \code{ctype} triggers either the Weibull sample estimator through the \emph{empirical copula}. The combination \code{as.sample=TRUE} and \cr \code{ctype="na"} will trigger the Weibull sample estimator as if \code{ctype="weibull"};}
  \item{uu}{;}
  \item{vv}{;}
  \item{trap.nan}{A logical ; and}
  \item{...}{Additional arguments to pass to the copula.}
}
\value{
  The \eqn{\beta^\diamond_\mathbf{C}} is returned.
}
\references{
Schmid, Friedrich, and Schmidt, Rafael, 2007, Nonparametric inference on multivariate versions of Blomqvist's beta and related measures of tail dependence: Metrika, v. 66, pp. 323--354, \doi{10.1007/s00184-006-0114-3}.
}
\author{W.H. Asquith}
\seealso{\code{\link{blomCOP}}
}
\examples{
blomCOP(  cop=PSP) # [1] 0.3333333
blomCOPss(cop=PSP) # [1] 0.3333333

\dontrun{
para <- c(3, 6) # define parameters of two-parameter GHcop
UV <- simCOP(1000, cop=GHcop, para=para) # simulate to show general structure

# compute the tail dependencies from havling into the limits
taildepCOP(cop=GHcop, para=para, plot=TRUE)
# lower tail dependency = 0.96222
# upper tail dependency = 0.74088
# The two parameters influence how strongly the tail dependencies are.

# Schmid and Schmidt (2007, eq. 24) define the lower-tail dependency in
# terms of the Beta and p-->0 Beta(c(p,p), c(1,1)). Lets compute these and
# produce content suitable to show on the tail-dependency plot that the
# assertion for the lower dependency by Beta() is correct, which it is and then
# extend to the upper-tail dependency parameter that the authors seem to not
# have defined.
usr <- par()$usr[1:2]         # grab horizontal edges of the plot, and set up the
uuLO <- rep(pnorm(usr[1]), 2) # the uu for the lower tail and the vv for the upper
vvUP <- rep(pnorm(usr[2]), 2) # tail and then plot both with overplotting symbols
# lower-tail estimate and see how it plots along the value from taildepCOP()
SchmidsL <- blomCOPss(cop=GHcop, para=para, uu=uuLO, vv=c(1,1))
points(usr[1], SchmidsL, col="darkgreen", cex=2, pch=1, lwd=2)
points(usr[1], SchmidsL, col="darkgreen", cex=2, pch=3, lwd=2)
points(usr[1], SchmidsL, col="darkgreen", cex=2, pch=4, lwd=2)
# upper-tail estimate and see how it plots along the value from taildepCOP()
SchmidsU <- blomCOPss(cop=GHcop, para=para, uu=c(0,0), vv=vvUP)
points(usr[2], SchmidsU, col="darkgreen", cex=2, pch=1, lwd=2)
points(usr[2], SchmidsU, col="darkgreen", cex=2, pch=3, lwd=2)
points(usr[2], SchmidsU, col="darkgreen", cex=2, pch=4, lwd=2)

# now for sake of curiosity, let us see how the trajectory of the Blomqvist
# (Schmid--Schmidt) Betas at arriving at the tail dependencies as p-->0|1.
psl <- pnorm(seq(0, usr[1], by=-diff(range(c(0, usr[1])))/1000))
lines(qnorm(psl), sapply(psl, function(p) {
              blomCOPss(cop=GHcop, para=para, uu=rep(p, 2), vv=c(1,1)) }),
      col="darkgreen", lty=2)
psu <- pnorm(seq(0, usr[2], by= diff(range(c(0, usr[2])))/1000))
lines(qnorm(psu), sapply(psu, function(p) {
              blomCOPss(cop=GHcop, para=para, uu=c(0,0), vv=rep(p, 2)) }),
      col="darkgreen", lty=2) # }
}
\keyword{copula (characteristics)}
\keyword{copula (properties)}
\keyword{Blomqvist beta}
\keyword{Blomqvist (Schmid-Schmidt) beta}
