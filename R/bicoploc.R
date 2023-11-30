"bicoploc" <-
function(xp, yp=NULL, xout=NA, xpara=NULL, ypara=NULL, dtypex="nor", dtypey="nor",
         ctype=c("weibull", "hazen", "bernstein", "checkerboard"), kumaraswamy=TRUE,
         plotuv=TRUE, plotxy=TRUE, adduv=FALSE, addxy=FALSE, snv=FALSE, limout=TRUE,
         autoleg=TRUE, xleg="topleft", yleg=NULL, rugxy=TRUE, ruglwd=0.5,
         xlim=NULL, ylim=NULL, a=0, ff=pnorm(seq(-5, +5, by=0.1)), locdigits=6,
         verbose=TRUE, x=NULL, y=NULL, ...) {

  use_parametric_copula <- TRUE

  lo <- .Machine$double.eps; hi <- 1 - lo
  # being at and close to the edges of probability helps to ensure that we by default stress the
  # algorithm in how we handle or trap issues on the edges.
  ff <- sort( unique( c(0, lo, lo^0.5, ff, 1-lo^0.5, hi, 1) ) ) # assurance policy

  if(is.matrix(xp) | is.data.frame(xp)) {
    xp <- xp[,1]; yp <- yp[,2]
  }

  if(any(is.na(xp))) {
    message("some of the xp are NA, this function is not configured to handle such circumstances")
    return(NULL)
  }
  if(any(is.na(yp))) {
    message("some of the yp are NA, this function is not configured to handle such circumstances")
    return(NULL)
  }
  if(length(xp) != length(yp)) {
    message("length of xp and yp are not identical")
    return(NULL)
  }
  n <- length(xp)

  if(plotuv) adduv <- FALSE
  if(plotxy) addxy <- FALSE

  if(autoleg) kumaraswamy <- TRUE

  PPtxt    <- "plotting position"
  aa <- as.character(round(as.numeric(a), digits=4)) # to fit into built-in lmomco
  aa <- switch(aa, "0"      = "Weibull",
                   "0.3175" = "Median",
                   "0.375"  = "Blom",
                   "0.40"   = "Cunnane",
                   "0.44"   = "Gringorten",
                   "0.5"    = "Hazen")
  if(! is.null(aa)) PPtxt  <- paste0(aa, " ", PPtxt)

  init.kur <- c(1.5, 2) # initial guesses at Kumaraswamy parameters, these were somewhat guessed at
  # but seem have robustness. Maybe better could be found.

  ctype    <- match.arg(ctype)
  txt      <- unlist(strsplit(ctype, ""))
  ctypeTXT <- paste0(toupper(txt[1]), paste(txt[2:length(txt)], sep="", collapse="")) # extract first letter and upper case it
  ctypeTXTuv <- paste0("Diagonal inversion points in (U,V) domain strictly by ", ctypeTXT, " empirical copula")
  ctypeTXTxy <- paste0("Diagonal inversion points in (X,Y) domain strictly by ", ctypeTXT, " empirical copula")
  PCHctype <- 4   # a times symbol
  PCHcex   <- 1.2 # size of the letter for plotting
  PCHlwd   <- 1.8
  PCHcol   <- "mediumorchid1"        # symbol plotting color for the empirical copula diagonal
  PCHbg    <- "white"  # background symbol plotting color for the empirical copula diagonal
  COPpch   <- 23
  COPbg    <- "chocolate1"
  COPcol   <- "chocolate1"
  COPptlwd <- 0.8
  COPcex   <- 0.7
  COPlty   <- 1
  COPlwd   <- 6
  DIAeol   <- "darkgreen"  # color for plotting of kumaraswamy line smooth of diagonal
  DIAcol   <- "seagreen"  # color for plotting of kumaraswamy line smooth of diagonal
  DIAbg    <- "lightgreen"
  DIAcex   <- 0.91
  DIApch   <- 21       # background symbol plotting color for the diagonal
  DIAlwd   <- 3
  DIAlty   <- 1
  DIAptlwd <- 1.1
  OUTcex   <- 1.2
  OUTcol   <- "darkred"    #            symbol plotting color for the Xout, if given
  OUTbg    <- "white"  # background symbol plotting color for the Xout, if given
  OUTpch   <- 16       # plotting character for the predictions from the inverse copula diagonal
  OUTptlwd <- 1
  LOCpch   <- 22       # plotting character for line of organic correlation predictions
  DATcol   <- "black"  # symbol plotting color for the Xs
  DATbg    <- "white"  # background symbol plotting color for the Xs
  DATcex   <- 1.1
  DATpch   <- 21
  DATptlwd <- 0.9
  LOClwd   <- 3
  LOClty   <- 2
  LOCcol   <- "blue"
  LOCbg    <- "white"
  LOCcex   <- 1.8
  LOCptlwd <- 1
  LEGcex   <- 0.75    # legend() text expansion factor
  SEGlen   <- 3.5      # width of line segments in the legend()
  RUGcol <- "wheat3" # ticking color for rug() plotting

  FFedge  <- c(0.001, 0.999) # reasonably wide range when plotting marginal distributions, these
  qFFedge <- qnorm(FFedge)   # standard normal variates of reasonably deep in the tails
  qff     <- qnorm(ff)       # standard normal variates of the ff vector
  # remember, the ff vector is the nonexceedance probability of the diagonal of the copula and hence
  # is that joint probability for which we are playing with.

  i <- 0

  if(verbose) message("(", i <- i + 1, ") plotting positions for U and V, recall lmomco::pp() ties.method=first")
  u  <- lmomco::pp(xp, a=a, sort=FALSE) # plotting positions
  v  <- lmomco::pp(yp, a=a, sort=FALSE) # plotting positions
  UV <- data.frame(U=u, V=v) # to become the para argument for the empirical copula

  sortX <- sort(xp) # making other diagnostic plots to potentially come easier to assemble

  # parametric marginal distribution of X and Y (if parameter estimation is otherwise required)
  if(is.null(xpara)) {
    if(verbose) message("(", i <- i + 1,
                  ") starting parameter estimation for X parametric margin by L-moments for dtypex")
        xpara <- NULL
    if(is.null(x)) {
      thex <- xp
    } else {
      if(verbose) message("(", i <- i + 1, ") using the alternative X values in argument x for parameter estimation")
      thex <- x
    }
    try(xpara <- lmomco::lmom2par(lmomco::lmoms(thex[! is.na(thex)]), type=dtypex), silent=TRUE)
    if(is.null(xpara)) {
      message("parameter estimation for X by method of L-moments has failed, try another dtypex?")
      return(NULL)
    }
    if(is.null(dtypex)) {
      if(exists("type", xpara)) dtypex <- xpara$type
    }
  } else {
    if(is.null(x)) { # this secondary look at the x is for support of the rug() call and nothing more
      thex <- xp
    } else {
      thex <- x
    }
    if(! is.null(dtypex) & (dtypex != xpara$type)) {
      dtypex <- "--" # reset to empty to ensure that what was given is used from ypara forevermore
    }
  }

  if(is.null(ypara)) {
    if(verbose) message("(", i <- i + 1,
                  ") starting parameter estimation for Y parametric margin by L-moments for dtypex")
        ypara <- "--"
    if(is.null(y)) {
      they <- yp
    } else {
      if(verbose) message("(", i <- i + 1, ") using the alternative X values in argument x for parameter estimation")
      they <- y
    }
    try(ypara <- lmomco::lmom2par(lmomco::lmoms(they[! is.na(they)]), type=dtypey), silent=TRUE)
    if(is.null(ypara)) {
      message("parameter estimation for Y by method of L-moments has failed, try another dtypey?")
      return(NULL)
    }
    if(is.null(dtypey)) {
      if(exists("type", ypara)) dtypey <- ypara$type
    }
  } else { # this secondary look at the y is for support of the rug() call and nothing more
    if(is.null(y)) {
      they <- yp
    } else {
      they <- y
    }
    if(! is.null(dtypey) & (dtypey != ypara$type)) {
      dtypey <- "--" # reset to empty to ensure that what was given is used from ypara forevermore
    }
  }
  n_thex <- length(thex[! is.na(thex)])
  n_they <- length(they[! is.na(they)])

  if(verbose) message("(", i <- i + 1, ") computing line of organic correlation by lmomco::lmrloc()")
  loc      <- lmomco::lmrloc( data.frame(X=xp, Y=yp) )
  # compute the line of organic correlation (reduced major axis), using the method of L-moments
  if(verbose) message("(", i <- i + 1, ") estimating xout by line of organic correlation")
  Yloc     <- loc$loc_lmr[2] * sortX + loc$loc_lmr[1]
  Yloc_out <- loc$loc_lmr[2] * xout  + loc$loc_lmr[1]
  names(Yloc)     <- NULL
  names(Yloc_out) <- NULL
  #Yloc_unsorted <- loc$loc_lmr[2] * x + loc$loc_lmr[1]
  #Uloc <- lmomco::cdfnor(            x, lmomco::parnor(lmomco::par2lmom(xpara), checklmom=FALSE))
  #Vloc <- lmomco::cdfnor(Yloc_unsorted, lmomco::parnor(lmomco::par2lmom(ypara), checklmom=FALSE))
  # LOC ONLY PRESERVES THE MEAN AND STANDARD DEVIATION.

  # Invert the diagonals of copula. The ff is the joint probability Pr[X <= x & Y <= y] = C(u,v)
  # and we want to solve on the diagonal for U=V=t --> C(t,t) = ff on the supposition that the
  # diagonal of the copula is the mathematical structure useful for the problem at hand
  if(verbose) message("(", i <- i + 1, ") computing vector of primary diagonal of empirical copula")
  ec <- eco <- COP(ff,   ff, cop=EMPIRcop, para=UV, ctype=ctype, ...)
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
    kur.init.para <- lmomco::vec2par(init.kur, type="kur")
    eckurpara     <- lmomco::disfitqua(ec, ff,  type="kur", init.para=kur.init.para)
    ec <- lmomco::qlmomco(ff, eckurpara) # the replacement of the "smoothed diagonal of the empirical
    # copula through the use of the Kumaraswamy distribution fit to RMSE of method of percentiles.
    # The ec are the u=v=t in the general copula literature nomenclature.
  }

  if(verbose) message("(", i <- i + 1, ") determing global sign of association direction by Spearman Rho")
  # recall for the usual line of organic correlation that the sign of Rho sets the sign of the
  # slope of the line, it appears that we can do the same as part of the extension of the logic
  # with copulas to support negatively associated data sets

  ### BEGIN PARAMETRIC COPULA BLOCK
  permsynsim <- 5E4
  rhoS  <- cor(xp, yp, method="spearman")
  infS  <- LzCOPpermsym(cop=EMPIRcop, para=UV, n=permsynsim, type="halton", as.vec=FALSE, ctype=ctype)
  infSv <- LzCOPpermsym(cop=EMPIRcop, para=UV, n=permsynsim, type="halton", as.vec=TRUE,  ctype=ctype)
  tparf <- function(par) c(exp(par[1]), pnorm( par[2] ), pnorm( par[3] ))
  rparf <- function(par) c(log(par[1]), qnorm( par[2] ), qnorm( par[3] ))
  ofunc <- function(par) { # objective function
    mypara <- tparf(par)
    mypara <- list(cop=PLcop, para=mypara[1], alpha=mypara[2], beta=mypara[3])
    rhoT   <- rhoCOP(cop=composite1COP, para=mypara)    # simulated Spearman Rho
    infTv  <- LzCOPpermsym(cop=composite1COP, para=mypara, n=permsynsim, type="halton", as.vec=TRUE)
    (rhoT - rhoS)^2 + mean( (infTv - infSv)^2 )
  }
  init.par <- rparf(c(1, 0.5, 0.5)); rt <- NULL # initial parameter guess
  try( rt <- optim(init.par, fn=ofunc) ) # 3D optimization
  if(is.null(rt)) {
    message("composite1COP parameter estimation returned NULL, disabling parameteric copula")
    use_parametric_copula <- FALSE
  } else {
    para.cop <- tparf(rt$par)
    para.cop <- list(cop=PLcop, para=para.cop[1], alpha=para.cop[2], beta=para.cop[3])
    rhoT <- rhoCOP(      cop=composite1COP, para=para.cop)
    infT <- LzCOPpermsym(cop=composite1COP, para=para.cop, n=permsynsim, type="halton", as.vec=FALSE)

    #JK <- simCOP(1000, cop=composite1COP, para=para.cop)
    faqscop <- c(round(rhoS, digits=6), round(rhoT, digits=6),
                 round(infS, digits=6), round(infT, digits=6),
                 round(para.cop$para,  digits=6),
                 round(para.cop$alpha, digits=6), round(para.cop$beta, digits=6))
    names(faqscop) <- c("SpearmanRhoSample", "SpearmanRhoFittedCopula",
                        "LzCOPpermsymSample", "LzCOPpermsymFittedCopula",
                        "CopulaParameter", "Alpha", "Beta")
    dtt <- COP(ff, ff, cop=composite1COP, para=para.cop, ...)
  }
  ### END PARAMETRIC COPULA BLOCK
  dtt[dtt > 1] <- 1 # assurance in case any numerical issues have caused a nudge outside the
  dtt[dtt < 0] <- 0 # applicable range of probability on the closed interval [0,1]

  rhosign <- sign( rhoS )
  if(rhosign < 0) { # negative association
    ecy  <- 1 - ec
    ecoy <- 1 - eco
    dtty <- 1 - dtt
  } else {          # positive association
    ecy  <-     ec
    ecoy <-     eco
    dtty <-     dtt
  }

  # An open question, do we 1-v the plotting of the diagonal on the Y-axis when Rho < 0
  # and also a 1 - lmomco::quakur() also, or does all of this reflect things in the wrong way?
  if(plotuv) {
    if(! adduv) {
      if(snv) {
        xylim <- qFFedge
        plot(qff, qnorm(eco), xlim=xylim, ylim=xylim, type="n", las=1, xaxs="i", yaxs="i",
             xlab="SNV OF COPULA DIAGONAL (Vector of 'f' joint probability)", mgp=c(2.2, 0.75, 0),
             ylab="SNV INVERSE COPULA DIAGONAL (f = C(U,V))\nTHE U = V NONEXCEEDANCE PROBABILITY")
      } else {
        xylim <- c(0,1)
        plot(ff,       eco,  xlim=xylim, ylim=xylim, type="n", las=1, xaxs="r", yaxs="r",
             xlab="COPULA DIAGONAL (Vector of 'f' joint probability)", mgp=c(2.2, 0.75, 0),
             ylab="INVERSE COPULA DIAGONAL (f = C(U,V))\nTHE U = V NONEXCEEDANCE PROBABILITY")
      }
      axis(3, axTicks(1), labels=FALSE, lwd=NA, lwd.ticks=1)
      axis(4, axTicks(2), labels=FALSE, lwd=NA, lwd.ticks=1)
    }
    if(snv) {
      if(! is.null(para.cop)) lines(qff, qnorm(dtt), col=COPcol, lwd=COPlwd)
      points(    qff,     qnorm(eco),  pch=PCHctype, cex=PCHcex, col=PCHcol, lwd=PCHlwd)
      points(qnorm(UV$U), qnorm(UV$V), pch=DATpch,   cex=DATcex, col=DATcol, lwd=DATptlwd, bg=DATbg)
      if(kumaraswamy) lines( qff, qnorm(lmomco::qlmomco(ff, eckurpara)), col=DIAcol, lwd=DIAlwd, lend=2)
      points(    qff,     qnorm(ecy),  pch=DIApch,   cex=DIAcex, col=DIAeol, lwd=DIAptlwd, bg=DIAbg)
    } else {
      if(! is.null(para.cop)) lines( ff, dtt, col=COPcol, lwd=COPlwd)
      points(     ff,           eco,   pch=PCHctype, cex=PCHcex, col=PCHcol, lwd=PCHlwd)
      points(UV$U,      UV$V,          pch=DATpch,   cex=DATcex, col=DATcol, lwd=DATptlwd, bg=DATbg)
      if(kumaraswamy) lines(  ff, lmomco::qlmomco(ff, eckurpara),col=DIAcol, lwd=DIAlwd, lend=2)
      points(     ff,           ec,    pch=DIApch,   cex=DIAcex, col=DIAeol, lwd=DIAptlwd, bg=DIAbg)
    }
    if(autoleg) {
      lwd    <- c(COPlwd,    DIAlwd,        NA,        NA,        NA, NA)
      lty    <- c(COPlty,    DIAlty,        NA,        NA,        NA, NA)
      col    <- c(COPcol,    DIAcol,    DIAcol,    PCHcol,    DATcol, NA)
      pch    <- c(    NA,        NA,    DIApch,  PCHctype,    DATpch, NA)
      pt.bg  <- c(    NA,        NA,     DIAbg,     PCHbg,     DATbg, NA)
      pt.lwd <- c(    NA,        NA,  DIAptlwd,    PCHlwd,  DATptlwd, NA)
      pt.cex <- c(    NA,        NA,    DIAcex,    PCHcex,    DATcex, NA)
      txt    <- c("Diagonal inverse of parametric asymmetric Plackett copula",
                  "Kumaraswamy smooth for diagonal inverse of the empirical copula coordinates",
                  "Coordinates of empirical copula diagonal inverse from joint probability vector",
                   ctypeTXTuv,
       paste0("Paired data points (n=", prettyNum(n, big.mark=","), ") between U and V by ", PPtxt),
              'Note, interpretation of "best fit" to data points on this plot is not applicable.')
      opts <- par(no.readonly=TRUE)
      par(lend=2)
      legend(xleg, yleg, txt,  bty="o", cex=LEGcex, inset=0.01, seg.len=SEGlen,
             lwd=lwd, lty=lty, col=col, pch=pch, pt.bg=pt.bg, pt.lwd=pt.lwd, pt.cex=pt.cex,
             bg="white", box.col=NA)
      par(opts)
    }
  }

  if(verbose) message("(", i <- i + 1, ") solving for Y predictions along the diagonal inversion")
  # solve for Y, like missing record estimation using the best available information on the Y
  # distribution's parameters (ie, the MOMENTS)
  # Back us off from -Inf and +Inf to NaN results when a probability = 0 | 1, to ensure numerical
  tmp <- ecy;  tmp[tmp <= lo] <- lo; tmp[tmp >= hi] <- hi
  Ysec0  <- lmomco::par2qua(tmp,  ypara,  paracheck=FALSE)
  tmp <- ecoy; tmp[tmp <= lo] <- lo; tmp[tmp >= hi] <- hi
  Ysec0o <- lmomco::par2qua(tmp,  ypara,  paracheck=FALSE)
  tmp <- dtty; tmp[tmp <= lo] <- lo; tmp[tmp >= hi] <- hi
  Ysdtt0 <- lmomco::par2qua(tmp,  ypara,  paracheck=FALSE)

  # The suppressWarnings() because of this
  #   # Warning message:  # when ties exist in the incoming to approx().
  # In regularize.values(x, y, ties, missing(ties), na.rm = na.rm):collapsing to unique 'x' values
  #suppressWarnings( Ysdd <- approx(lmomco::par2qua(dd, xpara), Ysdd0, xout=sortX, rule=2)$y )

  # Back us off from -Inf and +Inf to NaN results when a probability = 0 | 1, to ensure numerical
  # return on the approx() application.
  if(verbose) message("(", i <- i + 1, ") solving for Y predictions along given Xs and given (if any) Xouts")
  tmp <- ec;  tmp[tmp <= lo]  <- lo; tmp[tmp >= hi] <- hi
  suppressWarnings( Ysec      <- approx(lmomco::par2qua(tmp, xpara), Ysec0,  xout=sortX, rule=2)$y )
  suppressWarnings( Ysec_out  <- approx(lmomco::par2qua(tmp, xpara), Ysec0,  xout=xout,  rule=2)$y )
  tmp <- eco; tmp[tmp <= lo]  <- lo; tmp[tmp >= hi] <- hi
  suppressWarnings( Yseco     <- approx(lmomco::par2qua(tmp, xpara), Ysec0o, xout=sortX, rule=2)$y )
  suppressWarnings( Ysec_outo <- approx(lmomco::par2qua(tmp, xpara), Ysec0o, xout=xout,  rule=2)$y )
  tmp <- dtt; tmp[tmp <= lo]  <- lo; tmp[tmp >= hi] <- hi
  suppressWarnings( Ysdtt     <- approx(lmomco::par2qua(tmp, xpara), Ysdtt0, xout=sortX, rule=2)$y )
  suppressWarnings( Ysdtt_out <- approx(lmomco::par2qua(tmp, xpara), Ysdtt0, xout=xout,  rule=2)$y )

  names(Ysec)      <- NULL
  names(Ysec_out)  <- NULL
  names(Ysec_outo) <- NULL

  if(plotxy) { # real-world unit plotting
    myxlim <- range(  c(xp, lmomco::qlmomco(FFedge, xpara)),          na.rm=TRUE)
    myylim <- range(  c(yp, lmomco::qlmomco(FFedge, ypara)),          na.rm=TRUE)
    if(limout) {
      myxlim <- range(c(myxlim, xout),                                na.rm=TRUE)
      myylim <- range(c(myylim, Ysec_out, Ysec_outo, Yloc, Yloc_out), na.rm=TRUE)
    }
    if(! is.null(xlim)) {
      if(length(xlim) == 2 & ! any(is.na(xlim))) {
        myxlim <- xlim
      } else {
        myxlim <- range(c(xlim, myxlim), na.rm=TRUE)
      }
    }
    if(! is.null(ylim)) {
      if(length(ylim) == 2 & ! any(is.na(ylim))) {
        myylim <- ylim
      } else {
        myylim <- range(c(ylim, myylim), na.rm=TRUE)
      }
    }
    if(! addxy) {
      plot(xp, yp, type="n", xlim=myxlim, ylim=myylim, las=1, mgp=c(2.2, 0.75, 0),
                xlab="X OF DISTRIBUTION", ylab="Y OF DISTRIBUTION")
      axis(3, axTicks(1), labels=FALSE, lwd=NA, lwd.ticks=1)
      axis(4, axTicks(2), labels=FALSE, lwd=NA, lwd.ticks=1)
    }
    if(rugxy) {
      rug(thex, ticksize=0.02, side=1, col=RUGcol, lwd=ruglwd)
      rug(they, ticksize=0.02, side=2, col=RUGcol, lwd=ruglwd)
    }
    points(xp, yp, pch=DATpch, cex=DATcex, col=DATcol, bg=DATbg, lwd=DATptlwd)
    abline(loc$loc_lmr[1], loc$loc_lmr[2], lty=LOClty, lwd=LOClwd, col=LOCcol, lend=2)
    lines(lmomco::par2qua(dtt, xpara), Ysdtt0, lty=COPlty, lwd=COPlwd, col=COPcol, lend=2)
    # disabled the dd because in reality we could not know the parent copula itself
    lines(lmomco::par2qua(ec,  xpara), Ysec0,  lwd=DIAlwd,  col=DIAcol, lend=2)
    points(sortX, Yseco,     pch=PCHctype, cex=PCHcex, col=PCHcol, lwd=PCHlwd,   bg=PCHbg)
    points(sortX, Ysec,      pch=DIApch,   cex=DIAcex, col=DIAcol, lwd=DIAptlwd, bg=DIAbg)
    points(xout,  Yloc_out,  pch=LOCpch,   cex=LOCcex, col=OUTcol, lwd=LOCptlwd, bg=LOCbg)
    points(xout,  Ysec_out,  pch=OUTpch,   cex=OUTcex, col=OUTcol, lwd=OUTptlwd, bg=LOCbg)
    points(xout,  Ysdtt_out, pch=COPpch,   cex=COPcex, col=COPcol, lwd=COPptlwd, bg=COPbg)

    if(autoleg) {
      lwd    <- c(LOClwd,    COPlwd,    DIAlwd,        NA,        NA,        NA,        NA,        NA,       NA)
      lty    <- c(LOClty,    COPlty,    DIAlty,        NA,        NA,        NA,        NA,        NA,       NA)
      col    <- c(LOCcol,    COPcol,    DIAcol,    DIAcol,    PCHcol,    DATcol,    OUTcol,    OUTcol,   COPcol)
      pch    <- c(NA,            NA,        NA,    DIApch,  PCHctype,    DATpch,    LOCpch,    OUTpch,   COPpch)
      pt.bg  <- c(NA,            NA,        NA,     DIAbg,     PCHbg,     DATbg,     OUTbg,     OUTbg,    COPbg)
      pt.lwd <- c(NA,             NA,       NA,  DIAptlwd,    PCHlwd,  DATptlwd,  LOCptlwd,  OUTptlwd, COPptlwd)
      pt.cex <- c(NA,            NA,        NA,    DIAcex,    PCHcex,    DATcex,    LOCcex,    OUTcex,   COPcex)
      txt    <- c("Conventional line of organic correlation (LOC) though fit using linear moments",
                  "Organic correlation by parametric asymmetric Plackett copula",
                  "Organic correlation by copula diagonal inverse (Kumaraswamy smooth)",
                  "Coordinates of organic correlation by empirical copula diagonal inverse",
                   ctypeTXTxy,
                  paste0("Paired data points (n=", prettyNum(n, big.mark=","), ") between X and Y"),
                  "LOC predicted Y for requested X",
                  "Organic correlation predicted Y for requested X through empirical copula",
                  "Organic correlation predicted Y for requested X through parametic copula")
      if(is.null(xout) | length(xout[! is.na(xout)]) == 0) {
        nix    <- -c(7, 8)
        txt    <- txt[   nix]
        lwd    <- lwd[   nix]
        lty    <- lty[   nix]
        col    <- col[   nix]
        pch    <- pch[   nix]
        pt.bg  <- pt.bg[ nix]
        pt.lwd <- pt.lwd[nix]
        pt.cex <- pt.cex[nix]
      }
      opts <- par(no.readonly=TRUE)
      par(lend=2)
      legend(xleg, yleg, txt,  bty="o", cex=LEGcex, inset=0.01, seg.len=SEGlen,
             lwd=lwd, lty=lty, col=col, pch=pch, pt.bg=pt.bg, pt.lwd=pt.lwd, pt.cex=pt.cex,
             bg="white", box.col=NA)
      par(opts)
    }
  }

  df <- data.frame(xout=xout, loc=Yloc_out,
                   bicoploc=Ysec_out, bicoploc_emp=Ysec_outo, bicoploc_cop=Ysdtt_out)
  df <- round(df, digits=locdigits)
        dtypes  <- c( dtypex,   dtypey )
  names(dtypes) <- c("dtypex", "dtypey")
  faqs <- c(n, n_thex, n_they)
  names(faqs) <- c("paired_sample_size", "nonmissing_x_sample_size", "nonmissing_y_sample_size")
  zz <- list(loc=df, xpara=xpara, ypara=ypara, dtypes=dtypes, faqs=faqs, faqscop=faqscop,
             diag=data.frame(jtprob=round(ff,  digits=8),
                             uv    =round(ec,  digits=8),
                             uv_emp=round(eco, digits=8),
                             uv_cop=round(dtt, digits=8)))
  return(zz)
}

#set.seed(4); nsim <- 50
#X  <- rnorm(nsim, mean=3, sd=0.6)
#Y  <- rnorm(nsim, mean=0, sd=0.2)
#zz <- bicoploc(X,Y, xout=c(2.5, 3.5, 4), dtypex="nor", dtypey="gov")
#
#set.seed(4); nsim <- 50
#X  <- rexp(nsim, rate=3)
#Y  <- 0.7*X + rnorm(nsim, mean=0, sd=0.2)
#zz <- bicoploc(X,Y, xout=c(2.5, 3.5, 4), dtypex="exp", dtypey="gev")

#set.seed(1); nsim <- 100; npair <- 50
#UV <- rCOP(nsim, cop=composite1COP, para=list(cop=PLcop, para=150, alpha=0.8, beta=0.3))
#X  <- lmomco::qlmomco(UV[,1], lmomco::vec2par(c(3, 0.6, +0.5), type="pe3"))
#Y  <- lmomco::qlmomco(UV[,2], lmomco::vec2par(c(3, 0.4, +0.0), type="pe3"))
#ix <- sample(seq_len(nsim), npair, replace=FALSE)
#Xp <- X[ix]; Yp <- Y[ix]; dtypex <- "gev"; dtypey <- "gev"
#xpara <- lmomco::lmr2par(X, type=dtypex); ypara <- lmomco::lmr2par(Y, type=dtypey)
#plot(10^X, 10^Y, log="xy", las=1,
#     xlab="SOME RISK PHENOMENON IN X-DIRECTION", ylab="SOME RISK PHENOMENON IN Y-DIRECTION")
#zz <- bicoploc(Xp, Yp, xout=c(1.5, 2.5, 3.5, 4), xpara=xpara, ypara=ypara, ylim=c(1.8, 4.5))
#zz <- bicoploc(Xp, Yp, xout=c(1.5, 2.5, 3.5, 4), xpara=xpara, ypara=ypara, ylim=c(1.8, 4.5), x=X, y=Y)
#zz <- bicoploc(X,  Y,  xout=c(1.5, 2.5, 3.5, 4), dtypex="pe3", dtypey="pe3", ylim=c(1.8,  4.5))

