"bicoploc" <-
function(x, y, xout=NA, xpara=NULL, ypara=NULL, dtypex="nor", dtypey="nor",
               ctype=ctype=c("weibull", "hazen", "bernstein", "checkerboard"),
               kumaraswamy=TRUE, plotonxy=TRUE, plotondiag=TRUE, a=0,
               tt=pnorm(seq(-5, +5, by=0.1)), ...) {

  ctype <- match.arg(ctype)

  FFedge <- c(0.001, 0.999) # reasonably wide range when plotting marginal distributions, these
  qFFedge <- qnorm(FFedge)
  qtt <- qnorm(tt)

  u <- lmomco::pp(x, a=a, sort=FALSE)
  v <- lmomco::pp(y, a=a, sort=FALSE)
  UV <- data.frame(U=u, V=v)
  sortX <- sort(x) # sort to make other diagnostic plots to potentially come easier to assemble
  if(is.null(xpara)) xpara <- lmomco::lmom2par(lmomco::lmoms(x), type=dtypex)
  if(is.null(ypara)) ypara <- lmomco::lmom2par(lmomco::lmoms(y), type=dtypey)
  loc   <- lmomco::lmrloc(data.frame(X=x, Y=y))
  # compute the line of organic correlation (reduced major axis)
  Yloc     <- loc$loc_lmr[2] * sortX + loc$loc_lmr
  Yloc_out <- loc$loc_lmr[2] * xout  + loc$loc_lmr[1]
  names(Yloc)     <- NULL
  names(Yloc_out) <- NULL
  # LOC ONLY PRESERVES THE MEAN AND STANDARD DEVIATION.

  # Invert the diagonals of copula. The tt is the joint probability Pr[X <= x & Y <= y] and we
  # want to solve on the diagonal for U=V=tt on the supposition that the diagonal of the copula
  # is the mathematical structure useful for the problem at hand
  ec <- eco <- diagCOPatf(tt, cop=EMPIRcop, para=UV, ctype=ctype)
  if(kumaraswamy) { # https://en.wikipedia.org/wiki/Kumaraswamy_distribution
    eckurpara <- disfitqua(ec, tt, type="kur", init.para=vec2par(c(1.5, 2), type="kur"))
    ec <- qlmomco(tt, eckurpara) # the replacement of the "smoothed diagonal of the empirical
    # copula through the use of the Kumaraswamy distribution fit to RMSE of method of percentiles.
  }

  if(plotondiag) {
    xylim <- qFFedge
    plot(qtt, qnorm(eco), xlim=xylim, ylim=xylim, type="p", las=1, pch="E", cex=0.8, col="brown",
         xlab="SNV OF COPULA DIAGONAL (Vector of 't' joint probability)",                 xaxs="i",
         ylab="SNV INVERSE COPULA DIAGONAL (t = C(U,V)) U = V NONEXCEEDANCE PROBABILITY", yaxs="i")
    points(qnorm(UV$U), qnorm(UV$V), pch=21, col=grey(0), bg=grey(0.95)) # UV current iteration
    axis(3, axTicks(1), labels=FALSE, lwd=NA, lwd.ticks=1)
    axis(4, axTicks(2), labels=FALSE, lwd=NA, lwd.ticks=1)
    if(kumaraswamy) lines( qtt, qnorm(qlmomco(tt, eckurpara)), col="forestgreen")
    points(qtt, qnorm(ec), pch=16, cex=0.7, col="forestgreen")
  }

  # solve for Y, like missing record estimation using the best available information on the Y
  # distribution's parameters (ie, the MOMENTS)
  #Ysdd0 <- lmomco::par2qua(dd,  ypara,  paracheck=FALSE)
  Ysec0 <- lmomco::par2qua(ec,  ypara,  paracheck=FALSE)
  #suppressWarnings( Ysdd <- approx(par2qua(dd, xpara), Ysdd0, xout=sortX, rule=2)$y )
  suppressWarnings( Ysec     <- approx(par2qua(ec, xpara), Ysec0, xout=sortX, rule=2)$y )
  suppressWarnings( Ysec_out <- approx(par2qua(ec, xpara), Ysec0, xout=xout,  rule=2)$y )
  # Warning message:   # What needs to be done? Certainly rule=2 is not quite right.
  # In regularize.values(x, y, ties, missing(ties), na.rm = na.rm):collapsing to unique 'x' values
  names(Ysec)     <- NULL
  names(Ysec_out) <- NULL
  xlim  <- range(c(x, qlmomco(FFedge, xpara)))
  ylim  <- range(c(y, qlmomco(FFedge, ypara)))
  if(plotonxy) { # intermediate plotting
    plot(x,y, type="p", xlim=xlim, ylim=ylim, cex=1.1, pch=21, col="black", bg="white", las=1,
               xlab="X OF DISTRIBUTION", ylab="Y OF DISTRIBUTION")
    axis(3, axTicks(1), labels=FALSE, lwd=NA, lwd.ticks=1)
    axis(4, axTicks(2), labels=FALSE, lwd=NA, lwd.ticks=1)
    abline(loc$loc_lmr[1], loc$loc_lmr[2], lty=2, lwd=3, col=grey(0.5))
    #lines(lmomco::par2qua(dd,  xparas), Ysdd0, lwd=5, col="chocolate1"   ) # disabled
    # disabled the dd because in reality we could not know the parent copula itself
    lines(lmomco::par2qua(ec,  xpara), Ysec0, lwd=2, col="forestgreen"  )
    #points(sortX, Ysdd, pch=16, cex=0.7, col="chocolate1")
    points(sortX, Ysec,     pch=16, cex=0.7, col="forestgreen")
    points(xout,  Yloc_out, pch=22, cex=1.5, lwd=1.5, col="red", bg="white")
    points(xout,  Ysec_out, pch=16, cex=1.2, col="red")
  }
  zz <- list(loc=data.frame(xout=xout, loc=Yloc_out, bicoploc=Ysec_out),
             diag=data.frame(jointprob=tt, uvempir=ec))
  return(zz)
}
