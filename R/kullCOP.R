"kullCOP" <-
function(cop1=NULL, cop2=NULL, para1=NULL, para2=NULL, alpha=0.05,
         del=0, n=1E5, verbose=TRUE, sobol=FALSE, scrambling=0, ...) {

    if(del > 0) {
      lo <- del; hi <- 1 - del # integration limits
    } else {
      lo <- 0; hi <- 1
    }

    UV <- NULL
    if(sobol) {
       if(! exists(".Random.seed")) tmp <- runif(1) # insures definition
       seed <- sample(.Random.seed, 1)
       UV <- randtoolbox::sobol(n=n, dim=2, seed=seed, scrambling=scrambling)
    } else {
       UV <- matrix(data=runif(2*n), ncol=2)
    }
    if(verbose) message("kullCOP: Computing 'f' density values---",
                        appendLF=FALSE)
    f <- densityCOP(UV[,1],UV[,2], cop=cop1, para=para1, ...)
    f[f == 0] <- .Machine$double.eps
    if(verbose) message("done")

    if(verbose) message("kullCOP: Computing 'g' density values---",
                        appendLF=FALSE)
    g <- densityCOP(UV[,1],UV[,2], cop=cop2, para=para2, ...)
    g[g == 0] <- .Machine$double.eps
    if(verbose) message("done")

    h <- g*(log(g/f))
    h <- h[is.finite(h)]
    h <- h[!  is.nan(h)]
    KLdivergence.fg    <- mean(h)
    KLdivergence.fg.sd <-   sd(h)/sqrt(n)

    h <- f*(log(f/g))
    h <- h[is.finite(h)]
    h <- h[!  is.nan(h)]
    KLdivergence.gf    <- mean(h)
    KLdivergence.gf.sd <-   sd(h)/sqrt(n)

    JEFF.divergence <- KLdivergence.fg + KLdivergence.gf
    names(JEFF.divergence) <- "Jeffrey Divergence"

    h <- g*log(g/f)^2
    h <- h[is.finite(h)]
    h <- h[!  is.nan(h)]
    KLvar.fg <- mean(h)
    KLvar.fg.sd <- sd(h)/sqrt(n)

    h <- f*log(f/g)^2
    h <- h[is.finite(h)]
    h <- h[!  is.nan(h)]
    KLvar.gf <- mean(h)
    KLvar.gf.sd <- sd(h)/sqrt(n)

    #    Warning messages:
    #1: In sqrt(KLvar.fg - KLdivergence.fg^2) : NaNs produced
    #2: In kullCOP(cop1 = EMPIRcop, para1 = UVaem, cop2 = EMPIRcop, para2 = UVnoaem) :
    #  NAs introduced by coercion to integer range

    suppressWarnings(sigmaKL.fg <- sqrt(KLvar.fg - KLdivergence.fg^2))
    if(is.nan(sigmaKL.fg)) {
      sigmaKL.fg <- 0
    }
    # message("KLvar.fg=",        KLvar.fg)
    # message("KLdivergence.fg=", KLdivergence.fg)

    sigmaKL.gf <- sqrt(KLvar.gf - KLdivergence.gf^2)
    if(is.nan(sigmaKL.gf)) {
      sigmaKL.gf <- 0
    }
    # message("KLvar.gf=",KLvar.gf)
    # message("KLdivergence.gt=",KLdivergence.gf)

    tmp <- c(sigmaKL.fg/KLdivergence.fg, sigmaKL.gf/KLdivergence.gf)
    tmp <- max(tmp[! is.nan(tmp)]) # even need with other protections above?

    KL.sample.size <- (qnorm(1-alpha) * tmp)^2

    diverge   <- c(KLdivergence.fg, sigmaKL.fg,
                   KLdivergence.gf, sigmaKL.gf)
    divergesd <- c(KLdivergence.fg.sd, KLvar.fg.sd,
                   KLdivergence.gf.sd, KLvar.gf.sd)
    names(diverge)   <- c("KL-diverge.fg", "sigmaKL-diverge.fg",
                          "KL-diverge.gf", "sigmaKL-diverge.gf")
    names(divergesd) <- c("StdDev_KL-diverge.fg", "StdDev_KL-variance.fg",
                          "StdDev_KL-diverge.gf", "StdDev_KL-variance.gf")
    SS <- as.integer(KL.sample.size)
    names(SS) <- "Kullback-Leibler sample size"
    zz <- list(MonteCarlo.sim.size = n,
               divergences         = diverge,
               stdev.divergences   = divergesd,
               Jeffrey.divergence  = JEFF.divergence,
               KL.sample.size      = SS)
    return(zz)
}


"kullCOPint" <-
 function(cop1=NULL, cop2=NULL, para1=NULL, para2=NULL, alpha=0.05,
          del=.Machine$double.eps^0.25, verbose=TRUE, ...) {

    if(del > 0) {
      lo <- del; hi <- 1 - del # integration limits
    } else {
      lo <- 0; hi <- 1
    }

    # Solve equation 5.7 in Joe (2015)
    if(verbose) message("  CPU on Kullback-Leibler, double integrations for divergence ",
                        "(f relative to g)")
    KLfg <- NULL
    try(KLfg <- integrate(function(u) {
        sapply(u, function(u) {
            integrate(function(v) {
              f <- densityCOP(u,v, cop=cop1, para=para1, ...)
              g <- densityCOP(u,v, cop=cop2, para=para2, ...)
              return(g*(log(g/f)))
            }, lo, hi)$value
        })
    }, lo, hi))
    #print(KLfg);stop()
    if(is.null(KLfg)) {
       warning("could not numerically integrate, Kullback-Leibler divergence, ",
               "(f relative to g)")
       KLfg            <- NA
       KLdivergence.fg <- NA
    } else {
       KLdivergence.fg <- KLfg$value
    }

    # Solve equation 5.8 in Joe (2015)
    if(verbose) message("  CPU on Kullback-Leibler, double integrations for divergence ",
                        "(g relative to f)")
    KLgf <- NULL
    try(KLgf <- integrate(function(u) {
        sapply(u, function(u) {
            integrate(function(v) {
              f <- densityCOP(u,v, cop=cop1, para=para1, ...)
              g <- densityCOP(u,v, cop=cop2, para=para2, ...)
              return(f*log(f/g))
            }, lo, hi)$value
        })
    }, lo, hi))
    if(is.null(KLgf)) {
       warning("could not numerically integrate, Kullback-Leibler divergence ",
               "(g relative to f)")
       KLgf            <- NA
       KLdivergence.gf <- NA
    } else {
       KLdivergence.gf <- KLgf$value
    }

    # Solve Jeffrey's Divergence Joe (2015)
    JEFF.divergence <- KLdivergence.fg + KLdivergence.gf
    names(JEFF.divergence) <- "Jeffrey Divergence"
    #JEFF <- NULL
    #try(JEFF <- integrate(function(u) {
    #    sapply(u, function(u) {
    #        integrate(function(v) {
    #          f <- densityCOP(u,v, cop=cop1, para=para1, ...)
    #          g <- densityCOP(u,v, cop=cop2, para=para2, ...)
    #          return((g-f)*log(g/f))
    #        }, lo, hi)$value
    #    })
    #}, lo, hi))
    #print(JEFF)


    # Solve equation 5.9 in Joe (2015)
    if(verbose) message("  CPU on Kullback-Leibler, double integrations for variance ",
                        "(f relative to g)")
    KLvar.fg <- NULL
    try(KLvar.fg <- integrate(function(u) {
        sapply(u, function(u) {
            integrate(function(v) {
              f <- densityCOP(u,v, cop=cop1, para=para1, ...)
              g <- densityCOP(u,v, cop=cop2, para=para2, ...)
              return(g*log(g/f)^2)
            }, lo, hi)$value
        })
    }, lo, hi))
    if(is.null(KLvar.fg)) {
       warning("could not numerically integrate, Kullback-Leibler variance ",
               "(f relative to g)")
       KLvar.fg   <- NA
       sigmaKL.fg <- NA
    } else {
       sigmaKL.fg <- sqrt(KLvar.fg$value - KLdivergence.fg^2)
    }

    if(verbose) message("  CPU on Kullback-Leibler, double integrations for variance ",
                        "(g relative to f)")
    KLvar.gf <- NULL
    try(KLvar.gf <- integrate(function(u) {
        sapply(u, function(u) {
            integrate(function(v) {
              f <- densityCOP(u,v, cop=cop1, para=para1, ...)
              g <- densityCOP(u,v, cop=cop2, para=para2, ...)
              return(f*log(f/g)^2)
            }, lo, hi)$value
        })
    }, lo, hi))
    if(is.null(KLvar.gf)) {
       warning("could not numerically integrate, Kullback-Leibler variance ",
               "(g relative to f)")
       KLvar.gf   <- NA
       sigmaKL.gf <- NA
    } else {
       sigmaKL.gf <- sqrt(KLvar.gf$value - KLdivergence.gf^2)
    }

    alpha <- alpha
    tmp <- max(c(sigmaKL.fg/KLdivergence.fg, sigmaKL.gf/KLdivergence.gf))
    KL.sample.size <- (qnorm(1-alpha) * tmp)^2

    yy <- list(KLintegrate.fg    = KLfg,
               KLintegrate.gf    = KLgf,
               KLvarintegrate.fg = KLvar.fg,
               KLvarintegrate.gf = KLvar.gf)
    diverge <- c(KLdivergence.fg, sigmaKL.fg,
                 KLdivergence.gf, sigmaKL.gf)
    names(diverge) <- c("KL-diverge.fg", "sigmaKL-diverge.fg",
                        "KL-diverge.gf", "sigmaKL-diverge.gf")
    SS <- as.integer(KL.sample.size)
    names(SS) <- "Kullback-Leibler sample size"
    zz <- list(divergences = diverge,
               Jeffrey.divergence=JEFF.divergence,
               KL.sample.size=SS,
               integrations = yy)
    return(zz)
}
