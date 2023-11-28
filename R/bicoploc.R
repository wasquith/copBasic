"bicoploc" <-
function(x, y, xout=NA, xpara=NULL, ypara=NULL, dtypex="nor", dtypey="nor",
               ctype=c("weibull", "hazen", "bernstein", "checkerboard"), kumaraswamy=TRUE,
               plotonuv=TRUE, plotonxy=TRUE, adduv=FALSE, addxy=FALSE, snv=FALSE,
               a=0, ff=pnorm(seq(-5, +5, by=0.05)), verbose=TRUE, ...) {

  lo <- .Machine$double.eps; hi <- 1 - lo
  # being at and close to the edges of probability helps to ensure that we by default stress the
  # algorithm in how we handle or trap issues on the edges.
  ff <- sort( unique( c(0, lo, lo^0.5, ff, 1-lo^0.5, hi, 1) ) ) # assurance policy

  if(is.matrix(x) | is.data.frame(x)) {
    x <- x[,1]; y <- y[,2]
  }

  if(plotonuv) adduv <- FALSE
  if(plotonxy) addxy <- FALSE

  ctype    <- match.arg(ctype)
  ctypePCH <- toupper(unlist(strsplit(ctype, ""))[1]) # extract first letter and upper case it
  PCHcex   <- 1.1 # size of the letter for plotting
  PCHcol   <- "turquoise"        # symbol plotting color for the empirical copula diagonal
  DIAcol   <- "forestgreen"  # color for plotting of kumaraswamy points and line smooth of diagonal
  OUTcol   <- "red"    #            symbol plotting color for the Xout, if given
  OUTbg    <- "white"  # background symbol plotting color for the Xout, if given
  DATcol   <- "black"
  DATbg    <- "white"  # background symbol plotting color for the Xs, if given

  FFedge  <- c(0.001, 0.999) # reasonably wide range when plotting marginal distributions, these
  qFFedge <- qnorm(FFedge)   # standard normal variates of reasonably deep in the tails
  qff     <- qnorm(ff)       # standard normal variates of the ff vector
  # remember, the ff vector is the nonexceedance probability of the diagonal of the copula and hence
  # is that joint probability for which we are playing with.

  if(verbose) message("(1) plotting positions for U and V, recall lmomco::pp() ties.method=first")
  u  <- lmomco::pp(x, a=a, sort=FALSE) # plotting positions
  v  <- lmomco::pp(y, a=a, sort=FALSE) # plotting positions
  UV <- data.frame(U=u, V=v) # to become the para argument for the empirical copula

  sortX <- sort(x) # making other diagnostic plots to potentially come easier to assemble

  # parametric marginal distribution of X and Y (if parameter estimation is otherwise required)
  if(is.null(xpara)) {
    if(verbose) message("(2) starting parameter estimation for X by L-moments for dtypex")
        xpara <- NULL
    try(xpara <- lmomco::lmom2par(lmomco::lmoms(x), type=dtypex), silent=TRUE)
    if(is.null(xpara)) {
      message("parameter estimation for X by method of L-moments has failed, try another dtypex?")
      return(NULL)
    }
  }
  if(is.null(ypara)) {
    if(verbose) message("(3) starting parameter estimation for Y by L-moments for dtypex")
        ypara <- NULL
    try(ypara <- lmomco::lmom2par(lmomco::lmoms(y), type=dtypey), silent=TRUE)
    if(is.null(ypara)) {
      message("parameter estimation for Y by method of L-moments has failed, try another dtypey?")
      return(NULL)
    }
  }

  if(verbose) message("(4) computing line of organic correlation by lmomco::lmrloc()")
  loc      <- lmomco::lmrloc( data.frame(X=x, Y=y) )
  # compute the line of organic correlation (reduced major axis), using the method of L-moments
  if(verbose) message("(5) estimating xout by line of organic correlation")
  Yloc     <- loc$loc_lmr[2] * sortX + loc$loc_lmr[1]
  Yloc_out <- loc$loc_lmr[2] * xout  + loc$loc_lmr[1]
  names(Yloc)     <- NULL
  names(Yloc_out) <- NULL
  # LOC ONLY PRESERVES THE MEAN AND STANDARD DEVIATION.

  # Invert the diagonals of copula. The ff is the joint probability Pr[X <= x & Y <= y] and we
  # want to solve on the diagonal for U=V=tt on the supposition that the diagonal of the copula
  # is the mathematical structure useful for the problem at hand
  if(verbose) message("(6) computing vector of primary diagonal of empirical copula")
  ec <- eco <- COP(ff, ff, cop=EMPIRcop, para=UV, ctype=ctype, ...)
  # Nonexceedances of the U and V and a copy forming the "original empirical", we will optionally
  # smooth the ec (empirical copula) but never the eco (empirical copula original)

  # Here are two other ways to think about this diagonal and ultimately its inverse. Seems just hit
  # it at mass as in above COP() call is sufficient with fine enough resolution for approx()-like
  # calls to come or for setup of the Kumaraswamy [0,1]-bounded quantile function smooth estimation
  # to ensure that we are on the real number line.
  #ed <- diagCOP(cop=EMPIRcop, para=UV, ctype=ctype, ploton=FALSE, lines=FALSE, delt=0.001)
  #ed$t <- c(0, ed$t, 1); ed$diagcop <- c(0, ed$diagcop, 1)
  #ff <- ed$t; ec <- eco <- ed$diagcop
  #ec <- eco <- diagCOPatf(ff, cop=EMPIRcop, para=UV, ctype=ctype)
  if(kumaraswamy) { # https://en.wikipedia.org/wiki/Kumaraswamy_distribution
    kur.init.para <- lmomco::vec2par(c(1.5, 2), type="kur")
    eckurpara     <- lmomco::disfitqua(ec, ff,  type="kur", init.para=kur.init.para)
    ec <- lmomco::qlmomco(ff, eckurpara) # the replacement of the "smoothed diagonal of the empirical
    # copula through the use of the Kumaraswamy distribution fit to RMSE of method of percentiles.
    # The ec are the u=v=t in the general copula literature nomenclature.
  }

  if(verbose) message("(7) determing global sign of association direction by Spearman Rho")
  # recall for the usual line of organic correlation that the sign of Rho sets the sign of the
  # slope of the line, it appears that we can do the same as part of the extension of the logic
  # with copulas to support negatively associated data sets
  rhosign <- sign( cor(x,y, method="spearman") )
  if(rhosign < 0) { # negative association
    ecy  <- 1 - ec
    ecoy <- 1 - eco
  } else {          # positive association
    ecy  <-     ec
    ecoy <-     eco
  }

  # An open question, do we 1-v the plotting of the diagonal on the Y-axis when Rho < 0
  # and also a 1 - lmomco::quakur() also, or does all of this reflect things in the wrong way?
  if(plotonuv) {
    if(! adduv) {
      if(snv) {
        xylim <- qFFedge
        plot(qff, qnorm(eco), xlim=xylim, ylim=xylim, type="n", las=1, xaxs="i", yaxs="i",
             xlab="SNV OF COPULA DIAGONAL (Vector of 'f' joint probability)",
             ylab="SNV INVERSE COPULA DIAGONAL (f = C(U,V)) U = V NONEXCEEDANCE PROBABILITY")
      } else {
        xylim <- c(0,1)
        plot( ff,       eco,  xlim=xylim, ylim=xylim, type="n", las=1,
             xlab="COPULA DIAGONAL (Vector of 'f' joint probability)", xaxs="r", yaxs="r",
             ylab="INVERSE COPULA DIAGONAL (f = C(U,V)) U = V NONEXCEEDANCE PROBABILITY")
      }
      axis(3, axTicks(1), labels=FALSE, lwd=NA, lwd.ticks=1)
      axis(4, axTicks(2), labels=FALSE, lwd=NA, lwd.ticks=1)
    }
    if(snv) {
      points(    qff,     qnorm(eco),  pch=ctypePCH, cex=PCHcex, col=PCHcol)
      points(qnorm(UV$U), qnorm(UV$V), pch=21, col=grey(0), bg=grey(0.95))
      if(kumaraswamy) lines( qff, qnorm(lmomco::qlmomco(ff, eckurpara)), col=DIAcol)
      points(    qff,     qnorm(ecy),  pch=16, cex=0.7, col=DIAcol)
    } else {
      points(     ff,           eco,   pch=ctypePCH, cex=PCHcex, col=PCHcol)
      points(UV$U,      UV$V, pch=21,  col=grey(0), bg=grey(0.95))
      if(kumaraswamy) lines(  ff,       lmomco::qlmomco(ff, eckurpara),  col=DIAcol)
      points(     ff,           ec,    pch=16, cex=0.7, col=DIAcol)
    }
  }

  if(verbose) message("(8) solving for Y predictions along the diagonal inversion")
  # solve for Y, like missing record estimation using the best available information on the Y
  # distribution's parameters (ie, the MOMENTS)
  #Ysdd0 <- lmomco::par2qua(dd,  ypara,  paracheck=FALSE)
  # Back us off from -Inf and +Inf to NaN results when a probability = 0 | 1, to ensure numerical
  tmp <- ecy;  tmp[tmp <= lo] <- lo; tmp[tmp >= hi] <- hi
  Ysec0  <- lmomco::par2qua(tmp,  ypara,  paracheck=FALSE)
  tmp <- ecoy; tmp[tmp <= lo] <- lo; tmp[tmp >= hi] <- hi
  Ysec0o <- lmomco::par2qua(tmp,  ypara,  paracheck=FALSE)

  # The suppressWarnings() because of this
  #   # Warning message:  # when ties exist in the incoming to approx().
  # In regularize.values(x, y, ties, missing(ties), na.rm = na.rm):collapsing to unique 'x' values
  #suppressWarnings( Ysdd <- approx(lmomco::par2qua(dd, xpara), Ysdd0, xout=sortX, rule=2)$y )

  # Back us off from -Inf and +Inf to NaN results when a probability = 0 | 1, to ensure numerical
  # return on the approx() application.
  if(verbose) message("(9) solving for Y predictions along given Xs and given (if any) Xouts")
  tmp <- ec;  tmp[tmp <= lo] <- lo; tmp[tmp >= hi] <- hi
  suppressWarnings( Ysec      <- approx(lmomco::par2qua(tmp, xpara), Ysec0,  xout=sortX, rule=2)$y )
  suppressWarnings( Ysec_out  <- approx(lmomco::par2qua(tmp, xpara), Ysec0,  xout=xout,  rule=2)$y )
  tmp <- eco; tmp[tmp <= lo] <- lo; tmp[tmp >= hi] <- hi
  suppressWarnings( Yseco     <- approx(lmomco::par2qua(tmp, xpara), Ysec0o, xout=sortX, rule=2)$y )
  suppressWarnings( Ysec_outo <- approx(lmomco::par2qua(tmp, xpara), Ysec0o, xout=xout,  rule=2)$y )

  names(Ysec)      <- NULL
  names(Ysec_out)  <- NULL
  names(Ysec_outo) <- NULL

  if(plotonxy) { # real-world unit plotting
    xlim <- range(c(x, lmomco::qlmomco(FFedge, xpara)))
    ylim <- range(c(y, lmomco::qlmomco(FFedge, ypara)))
    if(! addxy) {
      plot(x,y, type="p", xlim=xlim, ylim=ylim, las=1, cex=1.1, pch=21, col=DATcol, bg=DATbg,
                xlab="X OF DISTRIBUTION", ylab="Y OF DISTRIBUTION")
      axis(3, axTicks(1), labels=FALSE, lwd=NA, lwd.ticks=1)
      axis(4, axTicks(2), labels=FALSE, lwd=NA, lwd.ticks=1)
    }
    points(x,y, pch=21, col=DATcol, bg=DATbg)
    abline(loc$loc_lmr[1], loc$loc_lmr[2], lty=2, lwd=3, col=grey(0.5))
    #lines(lmomco::par2qua(dd,  xparas), Ysdd0, lwd=5, col="chocolate1"   ) # disabled
    # disabled the dd because in reality we could not know the parent copula itself
    lines(lmomco::par2qua(ec,  xpara), Ysec0, lwd=2,  col=DIAcol)
    #points(sortX, Ysdd, pch=16, cex=0.7, col="chocolate1")
    points(sortX, Yseco,    pch=ctypePCH, cex=PCHcex, col=PCHcol)
    points(sortX, Ysec,     pch=16, cex=0.7, col=DIAcol)
    points(xout,  Yloc_out, pch=22, cex=1.5, lwd=1.5, col=OUTcol, bg=OUTbg)
    points(xout,  Ysec_out, pch=16, cex=1.2, col=OUTcol)
  }
  df <- data.frame(xout=xout, loc=Yloc_out, bicoploc=Ysec_out, bicoploc_emp=Ysec_outo)
  zz <- list(loc=df,
             diag=data.frame(jtprob=round(ff,  digits=16),
                             uv    =round(ec,  digits=16),
                             uv_emp=round(eco, digits=16)))
  return(zz)
}
