"wolfCOPtest2" <-
function(x, y, asuv=FALSE, aslist=TRUE, na.rm=TRUE, digits=6,
               probs=c(0.90, 0.95, 0.98, 0.99, 0.995), ...) {
  # The probs are quantile levels of the sigma to report, and these are useful to check against the
  # simulations but also to produce these as critical values should the user be interested in
  # these as well as the p-value.
  if(is.null(probs)      ) probs <- 0.95 # 95th percentile or rather the 5-percent critical value (upper tail).
  if( length(probs) == 0 ) probs <- 0.95 # 95th percentile or rather the 5-percent critical value (upper tail).
  if( length(probs) == "") probs <- 0.95 # 95th percentile or rather the 5-percent critical value (upper tail).

  if(length(x) == 1) { # If x is just one value, then it is treated as the Schweizer-Wolff Sigma
    rwolf <- x[1]; lwolf <- log(rwolf/(1-rwolf)); n <- y[1] # and the sample size is in y[1]
    if(n < 3) {
      warning("sample size is <3, returning NULL")
      return(NULL)
    }
  } else {
    if(length(x) != length(y)) {
      warning("length x != length y, returning NULL")
      return(NULL)
    }
    uv <- data.frame(u=x, v=y)
    if(na.rm) uv <- uv[complete.cases(uv),]
    n <- nrow(uv) # sample size
    if(n < 3) { # This handling of the sample size dates from an much earlier version of this function
      # that had a lower limit of 9. With the empirical distributions for sample sizes 3-40 now supported,
      # we drop the minimum sample size down to 3 but with the logic here, we effectively permit samples
      # sizes to be incoming down to
      warning("sample size is <3; returning NULL")
      return(NULL)
    }
    if(! asuv) { # if true, then the user has provided the paired observations of probability
      uv[,1] <- lmomco::pp(uv[,1], sort=FALSE, ...)
      uv[,2] <- lmomco::pp(uv[,2], sort=FALSE, ...)
    }

    rwolf <- wolfCOP(para=uv, as.sample=TRUE) # Schweizer-Wolff Sigma : wolf in (0,1)
    lwolf <- log(rwolf / (1 - rwolf)) # logit transform of the Sigma
  }

  dtype <- ifelse(n <= 40, "gno", "pe3") # We can see via inst/make_wolfCOPtest/chck_wolfCOPtestP.R
  # and the L-moment ratio diagram on the logit transform of the sigma, that there is a heuristic
  # change over from GNO to PE3 at about n=40. HOWEVER ------------------------------------------
  # Study of the terminal output of the make_wolfCOPtestP_smlsam.R script and the upper tail,
  # informs us that the PE3 is reaching into the upper tail already better than the GNO by sample
  # size of 6 to 9 at p-values as fine as 0.001, so we move to the pe3 throughout.
  dtype <- "pe3"

  # DISABLED if(bylogit) {
    # Nonlinear regression coefficients computed PRESS minimization of residuals for the
    # exponent on log10(sample size) term. The regressions come from simulation of the Sigma
    # distribution (its logit) assuming the Independence copula.
    mucoe <- c(0.04564608, -1.10684982, 1.02858995, -1.36093758)
    l2coe <- c(0.12772039, -0.00036535, 0.0954984, -2.036914068)
    t3coe <- c(0.08827982, 0.0032132, 0.15446976, -2.001953138)
    t4coe <- c(0.12558378, -0.00054097, 0.04265488, -2.856258)
  # DISABLED } else {
  # DISABLED   # Nonlinear regression coefficients computed PRESS minimization of residuals for the
  # DISABLED   # exponent on log10(sample size) term. The regressions come from simulation of the Sigma
  # DISABLED   # distribution not using its logit assuming the Independence copula.
  # DISABLED   mucoe <- c(-5.13011976, 0.17965821, 5.43973795, -0.14481258)
  # DISABLED   l2coe <- c(0.1117815, -0.07903764, 0.02014877, 1.71258)
  # DISABLED   t3coe <- c(0.22033934, 0.00035634, 0.00707766, -2.8248)
  # DISABLED   t4coe <- c(0.16112678, -0.00011312, -0.01705013, -2.9168)
  # DISABLED }

  # Apply the regressions using the hardwired coefficients herein
  # m <- 3000 # sample sizes for which we declare that Tau3 and Tau4 have become constant, which is
  # is technically close to reality but with the curvilinear regression being used, we eschew the
  # prediction not being monotonic decreasing with sample size want it to have an apparent asymptote.
  m <- ifelse(n > 100000, 100000, n) # This keeps the apparent trajectory of a Tau3 and Tau4 plot
  # having a hook in it as sample sizes increases to infinite.
  mu    <- mucoe[1] + mucoe[2] * log10(n) + mucoe[3] * log10(n)^mucoe[4] # Mean    (Lambda1)
  l2    <- l2coe[1] + l2coe[2] * log10(n) + l2coe[3] * log10(n)^l2coe[4] # Lambda2 (L-scale)
  t3    <- t3coe[1] + t3coe[2] * log10(m) + t3coe[3] * log10(m)^t3coe[4] # Tau3    (L-skew)
  t4    <- t4coe[1] + t4coe[2] * log10(m) + t4coe[3] * log10(m)^t4coe[4] # Tau4    (L-kurtosis)
  if(t4 < (5 * t3^2 - 1)/4) t4 <- (5 * t3^2 - 1)/4
  lmrs  <- c(mu, l2, t3, t4) # Tidy list of the Lmoments of the logit(Sigma) distribution
  if(dtype == "gno") {
    para  <-  lmomco::pargno(lmomco::vec2lmom(c(mu, l2, t3, t4)), useHosking=FALSE)
  } else if(dtype == "pe3") {
    para  <-  lmomco::parpe3(lmomco::vec2lmom(c(mu, l2, t3, t4)), useHosking=FALSE)
  } else {
    stop("should not be here in logic")
  }

  # Lmoment ratio diagram shows
  # Really close adherence to a generalized normal distribution (3-parameter log-normal)
  # and hence that distribution is chosen here.
  # DISABLED if(bylogit) {
    neps  <-               lmomco::par2cdf(lwolf, para, paracheck=FALSE)   # CDF of logit(SIGMA)
    quans <- 1 / (1 + exp(-lmomco::par2qua(probs, para, paracheck=FALSE))) # QUA of inverse logit --> quans are Sigmas
  # DISABLED } else { # The neps are the nonexceedance probabilities.
  # DISABLED   neps  <-               lmomco::par2cdf(rwolf, para, paracheck=FALSE)   # CDF of SIGMA
  # DISABLED   quans <-               lmomco::par2qua(probs, para, paracheck=FALSE)   # QUA no retransformation --> are Sigmas
  # DISABLED }

  quans     <- round(quans,     digits=digits) # Estimated upper-tail quantiles of Sigma distribution
  lmrs      <- round(lmrs,      digits=digits) # Lmoments of the logit(Sigma) distribution
  para$para <- round(para$para, digits=digits) # parameters of the fitted distribution
  rwolf     <- round(rwolf,     digits=digits) # Sigma as the test statistic
  lwolf     <- round(lwolf,     digits=digits) # logit(Sigma) as needed for p-value lookup

  # Because of the potential confusion with the logit() transform involved, let us have long names.
  names(lmrs)      <- c("mu_TEXT_sigmas",    "lscale_TEXT_sigmas",  "tau3_TEXT_sigmas", "tau4_TEXT_sigmas" )
  names(para$para) <- paste0("para", seq_len(length(para$para)), "_TEXT_sigmas")
  names(lmrs)      <- gsub("_TEXT_", "logit", names(lmrs     ))
  names(para$para) <- gsub("_TEXT_", "logit", names(para$para))
  # DISABLED names(lmrs)      <- gsub("_TEXT_", ifelse(bylogit, "logit", "_"), names(lmrs     ))
  # DISABLED names(para$para) <- gsub("_TEXT_", ifelse(bylogit, "logit", "_"), names(para$para))

  quatxt           <- paste0("fit_f", gsub("\\.", "p", as.character(100 * probs)))
  names(quans)     <- quatxt

  # We appear to use less system time on system.time(), when we find the file and load it instead
  # of using the data() call.
  smlsam <- system.file("data/wolfCOPtest_data_smlsam.RData", package="copBasic")
  max_n_in_smlsam <- 40 # max sample size within wolfCOPtest_data_smlsam$n (yes hard wired)
  if(n <= max_n_in_smlsam & file.exists(smlsam)) {
    wolfCOPtest_data_smlsam <- NULL # initialize whether or no so R CMD check --as-cran will pass by visibility
    load(smlsam)
    #data(wolfCOPtest_data_smlsam) # importFrom("utils", "data") # in NAMESPACE required
    sata <- wolfCOPtest_data_smlsam; suppressWarnings( rm(wolfCOPtest_data_smlsam) )
    sata <- sata[sata$n == n,]       # isolate the table to the sample size of interest
    if(nrow(sata) == 0) {
      warning("sample size n=", n, " does not exist in wolfCOPtest_data_smlsam, so set pval_small=NA")
      pval_small <- NA
      names(pval_small) <- paste0("p.value(sample_le", max_n_in_smlsam, ")")
    } else {
      sata <- sata[order(sata$probs),] # should be sorted already but do so again if needing to inspect
      row.names(sata) <- NULL; # print(sata, 16)
      suppressWarnings( nep_small <- approx(sata$wolfemp, y=sata$probs, xout=rwolf)$y )
      pval_small <- round(1 - nep_small, digits=16)
      names(pval_small) <- paste0("p.value(sample_le", max_n_in_smlsam, ")")
    }
  } else {
    pval_small <- NA; names(pval_small) <- paste0("p.value(sample_le", max_n_in_smlsam, ")")
  }

  pval <- round(1 - neps, digits=16); names(pval) <- paste0("p.value(dist_", dtype, ")")
  pval <- c(pval, pval_small)

  zz <- c(n, rwolf, lwolf, pval, para$para, lmrs, quans)
  names(zz) <- c("sample_size", "sigma", "logit_sigma",
                 names(pval), names(para$para), names(lmrs), quatxt)
  names(zz) <- gsub("_TEXT_", "logit", names(zz))
  # DISABLED names(zz) <- gsub("_TEXT_", ifelse(bylogit, "logit", "_"), names(zz))
  if(aslist) {
    wz <- c(rwolf, lwolf); names(wz) <- c("sigma", "logit_sigma")
    zz <- list(sample_size=n, statistic=wz, p.value=pval, distpara_by_lmoms=para$para)
    # DISABLED if(bylogit) {
      zz$lmoms_logit_sigma <- lmrs # L-moments of the logit(SIGMAS) distribution
    # DISABLED } else {
    # DISABLED   zz$lmoms_sigma       <- lmrs # L-moments of the       SIGMAS  distribution
    # DISABLED }
    zz$sigma_quantiles <- quans # Put these last because this length of vector is mutable, and it
    # visually makes these better on the right side of aslist=FALSE (vector return), in particular.
  }
  return(zz)
}



#CPU <- NULL
#for(n in as.integer(10^seq(log10(100), log10(1000), by=.05))) {
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

# probs <- pnorm(seq(-3.72, +3.72, by=0.02))
# para  <- vec2par(wolfCOPtest(0, 2000)$distpara_by_lmoms, type="pe3")
# quans <- 1 / (1 + exp(-lmomco::par2qua(probs, para)))
# plot(qnorm(probs), quans, type="l", ylim=c(0,1))

