"wolfCOPtest" <-
function(x, y, asuv=FALSE, aslist=TRUE, bylogit=TRUE, dtype="gno",
               na.rm=TRUE, digits=6, probs=c(0.90, 0.95, 0.98, 0.99, 0.995), ...) {
  # The probs are quantile levels of the sigma to report, and these are useful to check against the
  # simulations but also to produce these as critical values should the user be interested in
  # these as well as the p-value.
  if(is.null(probs))      probs <- 0.95 # 95th percentile or rather the 5-percent critical value (upper tail).
  if( length(probs) == 0) probs <- 0.95 # 95th percentile or rather the 5-percent critical value (upper tail).

  # -------------------------------
  dtype <- tolower(dtype)
  suppressWarnings( have_dtype <- lmomco::dist.list(dtype) )
  if(length(grep("does not match", as.character(have_dtype))) == 1) { # "The given type argument does not match a distribution") {
    warning("requested dtype='", dtype, "' does not exist in lmomco::dist.list() function")
    return(NULL)
  }
  if(have_dtype >= 5) {
    warning("function only supports up to 4-parameter distributions ")
  }
  if(have_dtype == "gld") {
    warning("function does support first four Lmoments, but a fifth would be needed for dtype='gld'")
  }
  # -------------------------------

  if(length(x) == 1) { # If x is just one value, then it is treated as the Schweizer-Wolff Sigma
    rwolf <- x[1]; lwolf <- log(rwolf/(1-rwolf)); n <- ngiven <- y[1] # and the sample size is in y[1]
    if(n < 9) {
      warning("sample size is <9; so, will treat it as if it were 9 and preserve original sample by negation")
      ngiven <- -n; n <- 9
    }
  } else {
    if(length(x) != length(y)) {
      warning("length x != length y, returning NULL")
      return(NULL)
    }
    uv <- data.frame(u=x, v=y)
    if(na.rm) uv <- uv[complete.cases(uv),]
    n <- ngiven <- nrow(uv) # sample size
    if(n < 9) {
      warning("sample size is <9; so, will treat it as if it were 9 and preserve original sample by negation")
      ngiven <- -n; n <- 9
    }
    if(! asuv) { # if true, then the user has provided the paired observations of probability
      uv[,1] <- lmomco::pp(uv[,1], sort=FALSE, ...)
      uv[,2] <- lmomco::pp(uv[,2], sort=FALSE, ...)
    }

    rwolf <- copBasic::wolfCOP(para=uv, as.sample=TRUE) # Schweizer-Wolff Sigma : wolf in (0,1)
    lwolf <- log(rwolf / (1 - rwolf)) # logit transform of the Sigma
  }

  if(bylogit) {
    # Nonlinear regression coefficients computed PRESS minimization of residuals for the
    # exponent on log10(sample size) term. The regressions come from simulation of the Sigma
    # distribution (its logit) assuming the Independence copula.
    mucoe <- c(-0.11655885, -1.07964311, 1.1612798, -1.19358)
    l2coe <- c(0.13256181, -0.00178554, 0.09075513, -2.06958)
    t3coe <- c(0.08710575, 0.00329644, 0.15671491, -1.98558)
    t4coe <- c(0.12818454, -0.00139074, 0.04120657, -2.9368)
  } else {
    # Nonlinear regression coefficients computed PRESS minimization of residuals for the
    # exponent on log10(sample size) term. The regressions come from simulation of the Sigma
    # distribution not using its logit assuming the Independence copula.
    mucoe <- c(-2.61652495, 0.14574984, 2.96319764, -0.26006258)
    l2coe <- c(0.13525809, -0.30994223, 0.22869554, 1.1418)
    t3coe <- c(0.22490525, -0.00120998, 0.00490342, -5.5048)
    t4coe <- c(0.16640906, -0.00180767, -0.02055019, -2.5168)
  }

  # Apply the regressions using the hardwired coefficients herein
  mu    <- mucoe[1] + mucoe[2] * log10(n) + mucoe[3] * log10(n)^mucoe[4] # Mean    (Lambda1)
  l2    <- l2coe[1] + l2coe[2] * log10(n) + l2coe[3] * log10(n)^l2coe[4] # Lambda2 (L-scale)
  t3    <- t3coe[1] + t3coe[2] * log10(n) + t3coe[3] * log10(n)^t3coe[4] # Tau3    (L-skew)
  t4    <- t4coe[1] + t4coe[2] * log10(n) + t4coe[3] * log10(n)^t4coe[4] # Tau4    (L-kurtosis)
  lmrs  <- c(mu, l2, t3, t4) # Tidy list of the Lmoments of the logit(Sigma) distribution
  para  <-  lmomco::lmom2par(lmomco::vec2lmom(c(mu, l2, t3, t4)), type=dtype) # Lmoment ratio diagram shows
  # Really close adherence to a generalized normal distribution (3-parameter log-normal)
  # and hence that distribution is chosen here.
  if(bylogit) {
    neps  <-               lmomco::par2cdf(lwolf, para, paracheck=FALSE)   # CDF of logit(SIGMA)
    quans <- 1 / (1 + exp(-lmomco::par2qua(probs, para, paracheck=FALSE))) # QUA of inverse logit --> quans are Sigmas
  } else { # The neps are the nonexceedance probabilities.
    neps  <-               lmomco::par2cdf(rwolf, para, paracheck=FALSE)   # CDF of SIGMA
    quans <-               lmomco::par2qua(probs, para, paracheck=FALSE)   # QUA no retransformation --> are Sigmas
  }

  quans     <- round(quans,     digits=digits) # Estimated upper-tail quantiles of Sigma distribution
  lmrs      <- round(lmrs,      digits=digits) # Lmoments of the logit(Sigma) distribution
  para$para <- round(para$para, digits=digits) # parameters of the fitted distribution
  rwolf     <- round(rwolf,     digits=digits) # Sigma as the test statistic
  lwolf     <- round(lwolf,     digits=digits) # logit(Sigma) as needed for p-value lookup

  # Because of the potential confusion with the logit() transform involved, let us have long names.
  names(lmrs)      <- c("mu_TEXT_sigmas",    "lscale_TEXT_sigmas",  "tau3_TEXT_sigmas", "tau4_TEXT_sigmas" )
  names(para$para) <- paste0("para", seq_len(length(para$para)), "_TEXT_sigmas")
  names(lmrs)      <- gsub("_TEXT_", ifelse(bylogit, "logit", "_"), names(lmrs     ))
  names(para$para) <- gsub("_TEXT_", ifelse(bylogit, "logit", "_"), names(para$para))

  quatxt           <- paste0("fit_f", gsub("\\.", "p", as.character(100 * probs)))
  names(quans)     <- quatxt

  pval <- round(1 - neps, digits=8); names(pval) <- "p.value"

  zz <- c(n, ngiven, rwolf, lwolf, pval, para$para, lmrs, quans)
  names(zz) <- c("sample_size_in_comps", "sample_size_given", "sigma", "logit_sigma", "p.value",
                 names(para$para), names(lmrs), quatxt)
  names(zz) <- gsub("_TEXT_", ifelse(bylogit, "logit", "_"), names(zz))
  if(aslist) {
    wz <- c(rwolf, lwolf); names(wz) <- c("sigma", "logit_sigma")
    zz <- list(sample_size_in_comps=n, sample_size_given=ngiven,
               sigma=wz, p.value=pval, distpara_by_lmoms=para$para)
    if(bylogit) {
      zz$lmoms_logit_sigma <- lmrs # L-moments of the logit(SIGMAS) distribution
    } else {
      zz$lmoms_sigma       <- lmrs # L-moments of the       SIGMAS  distribution
    }
    zz$sigma_quantiles <- quans # Put these last because this length of vector is mutable, and it
    # visually makes these better on the right side of aslist=FALSE (vector return), in particular.
  }
  return(zz)
}



#CPU <- NULL
#for(n in as.integer(10^seq(log10(700), log10(3000), by=.05))) {
#  UV <- simCOP(n=n, cop=CLcop, para=2, graphics=FALSE)
#  THwolf <- wolfCOP(cop=PSP, para=2, )
#  a <- system.time(MC <- wolfCOPsamc(UV[,1], UV[,2]))
#  b <- system.time(SA <- wolfCOP(para=UV, as.sample=TRUE))
#  CPU <- rbind(CPU, data.frame(n=n, a=a[3], b=b[3],
#                               mc=MC$estimates[2], sam=SA))
#}
#plot(  CPU$n, CPU$b, log="xy")
#points(CPU$n, CPU$a, pch=16)
#
#plot(  CPU$n, CPU$mc,  log="xy")
#points(CPU$n, CPU$sam, pch=16)
