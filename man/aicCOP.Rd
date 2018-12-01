\encoding{utf8}
\name{aicCOP}
\alias{aicCOP}
\title{ The Akaike Information Criterion between a Coupla and an Empirical Copula }
\description{
Compute the \emph{Akaike information criterion} \eqn{\mathrm{AIC}_\mathbf{C}} (Chen and Shenglian, 2019, p. 29), which is computed using \emph{mean square error} \eqn{\mathrm{MSE}_\mathbf{C}} as

\deqn{\mathrm{MSE}_\mathbf{C} = \frac{1}{n}\sum_{i=1}^n (P_{e_i} - P_i)^2\mbox{ and}}
\deqn{\mathrm{AIC}_\mathbf{C} = 2m + n\log(\mathrm{MSE}_\mathbf{C})\mbox{,}}

where \eqn{P_{e_i}} is the empirical joint probability for the \eqn{i}th observation, \eqn{P_i} is the joint probability from a fitted copula (itself having \eqn{m}) parameters. The empirical joint probability comes from \code{\link{EMPIRcop}}. The \eqn{\mathrm{AIC}_\mathbf{C}} is in effect saying that the best copula will have its joint probabilities plotting on a 1:1 line with the empirical joint probabilities, which is an \eqn{\mathrm{AIC}_\mathbf{C} = -\infty}. From the \eqn{\mathrm{MSE}_\mathbf{C}} shown above, the root mean square error \code{\link{rmseCOP}} and Bayesian information criterion \code{\link{bicCOP}} can be computed. These goodnesses of fit can assist in deciding one copula favorability over another.
}
\usage{
aicCOP(u,v, cop=NULL, para=NULL, ...)
}
\arguments{
  \item{u}{Nonexceedance probability \eqn{u} in the \eqn{X} direction;}
  \item{v}{Nonexceedance probability \eqn{v} in the \eqn{Y} direction;}
  \item{cop}{A copula function;}
  \item{para}{Vector of parameters or other data structure, if needed, to pass to the copula;}
  \item{...}{Additional arguments to pass to either copula (likely most commonly to the empirical copula).}
}
\value{
  The value for \eqn{\mathrm{AIC}_\mathbf{C}} is returned.
}
\references{
Chen, Lu, and Guo, Shenglian, 2019, Copulas and its application in hydrology and water resources: Springer Nature, Singapore, ISBN 978--981--13--0574--0.
}
\author{ W.H. Asquith}
\seealso{\code{\link{EMPIRcop}}, \code{\link{bicCOP}}, \code{\link{rmseCOP}}
}
\examples{
S <- simCOP(80, cop=GHcop, para=5) # Simulate some probabilities, but we
# must then treat these as data and recompute empirical probabilities. 
U <- lmomco::pp(S$U, sort=FALSE); V <- lmomco::pp(S$V, sort=FALSE)
# The parent distribution is Gumbel-Hougaard extreme value copula.
# But in practical application we don't know that but say we speculate that
# perhaps the Galambos extreme value might be the parent. Then maximum
# likelihood is used to fit the single parameter.
pGL <- mleCOP(U,V, cop=GLcop, interval=c(0,20))$par

aics <- c(aicCOP(U,V, cop=GLcop, para=pGL),
          aicCOP(U,V, cop=P),
          aicCOP(U,V, cop=PSP))
print(aics) # We will see that the first AIC is the smallest as the
# Galambos has the nearest overall behavior than the P and PSP copulas.
}
\keyword{copula (goodness of fit)}
