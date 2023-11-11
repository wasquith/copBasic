\encoding{utf8}
\name{ORDSUMcop}
\alias{ORDSUMcop}
\title{Ordinal Sums of Copula}
\description{
Compute ordinal sums of a copula (Nelsen, 2006, p. 63) within \eqn{\mathcal{I}^2} into \eqn{n} partitions (possibly infinite) within \eqn{\mathcal{I}^2}. According to Nelsen, letting \eqn{\mathcal{J}} denote a \emph{partition} \eqn{\mathcal{I}^2} and \eqn{\mathcal{J}_i = [a_i,\, b_i]} be the \eqn{i}th partition that does not overlap with others, let also \eqn{\mathbf{C}_i} be a copula for the \eqn{i}th partition, then the \emph{ordinal sum} of these \eqn{\mathbf{C}_i} with parameters \eqn{\Theta_i} \emph{with respect to} \eqn{\mathcal{J}_i} is the copula \eqn{\mathbf{C}} given by
\deqn{\mathbf{C}\bigl(u,v; \mathcal{J}_i, \mathbf{C}_i, \Theta_i, i \in 1,2,\cdots,n\bigr) = a_i + (b_i-a_i)\mathbf{C}_i\biggl(\frac{u-a_i}{b_i-a_i},\, \frac{u-a_i}{b_i-a_i};  \Theta_i,\biggr)\  \mathrm{for}\ (u,v) \in \mathcal{J}^2\mbox{,}}
for points within the partitions and otherwise outside the partitions as
\deqn{\mathbf{C}\bigl(u,v; \mathcal{J}_i, \mathbf{C}_i, i \in 1,2,\cdots,n\bigr) = \mathbf{M}(u,v)\  \mathrm{for}\ (u,v) \ni \mathcal{J}^2\mbox{.}} Finally, Nelsen (2006, theorem 3.2.1) states that a copula is an ordinal sum if and only if for a \eqn{t} if \eqn{\mathbf{C}(t,t)=t} for \eqn{t \in (0,1)}. The \emph{diagonal of a coupla} can be useful for quick assessment (see \bold{Examples}).
}
\usage{
ORDSUMcop(u,v, para=list(cop=W, para=NA, part=c(0,1)), ...)
}
\arguments{
  \item{u}{Nonexceedance probability \eqn{u} in the \eqn{X} direction;}
  \item{v}{Nonexceedance probability \eqn{v} in the \eqn{Y} direction;}
  \item{para}{A list of sublists for the coupla, parameters, and partitions (see \bold{Examples}) and some allowance for intelligent in fill of \code{para} is made within the sources; and}
  \item{...}{Additional arguments to pass.}
}
\value{
  Value(s) for the copula are returned.
}
\references{
Nelsen, R.B., 2006, An introduction to copulas: New York, Springer, 269 p.
}
\author{W.H. Asquith}
\seealso{\code{\link{copBasic-package}}}
\examples{
\dontrun{
  para = list(cop=c(CLcop, M, PLcop, GHcop), para=list(4, NA, 0.1, c(3,4)),
              part=list(c(0,0.25), c(0.25,0.35), c(0.35,0.85), c(0.85,1)))
  UV <- simCOP(n=100, cop=ORDSUMcop, para=para, ploton=FALSE)
  plot(c(0,1), c(0,1), xlab="U, NONEXCEEDANCE PROBABILITY",
                       ylab="V, NONEXCEEDANCE PROBABILITY")
  for(k in seq_len(length(para$part))) {
    a <- para$part[[k]][1]; b <- para$part[[k]][2]
    polygon(c(a, b, b, a, a), c(a,a,b,b,a), lty=2, lwd=0.8, col="lightgreen")
    text((a+b)/2, (a+b)/2, k, cex=3, col="blue")
  }
  points(UV, pch=21, cex=0.8, col=grey(0.1), bg="white") #}

\dontrun{
  para = list(cop=c(GHcop), para=list(c(2,3)), # internally replicated
              part=list(c(0,0.2), c(0.2,0.3), c(0.3,0.5), c(0.5,0.7), c(0.7,1)))
  UV <- simCOP(n=100, cop=ORDSUMcop, para=para, ploton=FALSE)
  plot(c(0,1), c(0,1), xlab="U, NONEXCEEDANCE PROBABILITY",
                       ylab="V, NONEXCEEDANCE PROBABILITY")
  for(k in seq_len(length(para$part))) {
    a <- para$part[[k]][1]; b <- para$part[[k]][2]
    polygon(c(a, b, b, a, a), c(a,a,b,b,a), lty=2, lwd=0.8, col="lightgreen")
    text((a+b)/2, (a+b)/2, k, cex=3, col="blue")
  }
 points(UV, pch=21, cex=0.8, col=grey(0.1), bg="white") #}

\dontrun{
  # In this example, it is important that the delt is of the resolution
  # matching the  edges of the partitions.
  para = list(cop=P, para=list(NULL),
              part=list(c(0,0.257), c(0.257,0.358), c(0.358,1)))
  DI <- diagCOP(cop=ORDSUMcop, para=para, delt=0.001)
  if(sum(DI$diagcop == DI$t) >= 1) {
    message("The ORDSUMcop() operation is an ordinal sum if there exists\n",
            "a t=(0,1) exists such that C(t,t)=t by Nelsen (2006, theorem 3.2.1).")
  }
  abline(0,1, col="red") #}
}
\keyword{copula (formulas)}
\concept{Ordinal Sums}
\keyword{copula}
\keyword{ordinal sums}

