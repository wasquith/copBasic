if(Sys.getenv("RSTUDIO") == "1") {
  # Automatic change directory to location of this script when RStudio is running.
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
} else {
  print("RStudio is not running")
}

library(copBasic)

SHOW_PLOT <- FALSE

FCORE <- "mc_wolfPI"
RUNNO  <- 127:127
ns <- unique( as.integer( 10^(c(seq(log10(10), log10(1000), by=0.05)) ) ) )
ns <- rep(c(1200), 1)
#stop()
for(runno in RUNNO) {
  OUTFILE <- paste0(FCORE, "_", runno, ".txt")
  if(file.exists(OUTFILE)) {
    if(runno == 9) stop("all output files are accounted for")
    next
  }
  cat("START\n", file=OUTFILE)
  break
}
message("writing to ", OUTFILE)


d <- 16
nsim <- 100000

x <- seq(min(ns), max(ns), by=1)

par(lend=1, ljoin=1, las=1, mar=c(5, 5, 4, 5)+0.1) # c(bottom, left, top, right)
D <- NULL
for(n in ns) {
  mynsim <- ceiling(10*nsim/n)
  if(mynsim == 1  ) mynsim <- 2
  if(n <= 200) mynsim <-  2000
  if(n <= 100) mynsim <-  3000
  if(n <= 80 ) mynsim <-  5000
  if(n <= 40 ) mynsim <- 10000
  if(n <= 20 ) mynsim <- 20000
  if(mynsim < 1/(1 - 0.995)) mynsim <- 1000
  mynsim <- 1000
  print(c(n, mynsim))
  #wolfP <- replicate(mynsim, { copBasic::wolfCOP(para=simCOP(n=n, cop=P, graphics=FALSE),
  #                                               as.sample=TRUE) })
  wolfP <- rep(NA, mynsim)
  for(i in seq_len(mynsim)) {
    print(system.time(wolfP[i] <- copBasic::wolfCOP(para=as.data.frame(matrix(runif(n*2), ncol=2)), as.sample=TRUE)))
    print(c(mynsim, n, i, wolfP[i]))
  }
  #for(i in seq_len(mynsim)) {
  #  results <- copBasic::wolfCOPsamc(para=as.data.frame(
  #                              matrix(runif(mynsim*2), ncol=2)))
  #  wolfP[i] <- results$estimates[2]
  #  print(c(i, results$its, wolfP[i]))
  #}
  lmr <- Lmoments::Lcoefs( wolfP, rmax=5)
  print(lmr)
  lmc <- lmomco::lmoms.cov(wolfP, se="lmrse")
  print(lmc)
  emp <- round(quantile(wolfP, probs=c(0.9, 0.95, 0.98, 0.99, 0.995), type=6), digits=d)
  df <- data.frame(nsim=mynsim, n=n,
                    mu=round(          lmr[1],    digits=d),
                   var=round((sqrt(pi)*lmr[2])^2, digits=d),
                  lam2=round(          lmr[2],    digits=d),
                  tau3=round(          lmr[3],    digits=d),
                  tau4=round(          lmr[4],    digits=d),
                  tau5=round(          lmr[5],    digits=d),
                  muse=round(          lmc[1],    digits=d),
                  l2se=round(          lmc[2],    digits=d),
                  t3se=round(          lmc[3],    digits=d),
                  t4se=round(          lmc[4],    digits=d),
                  t5se=round(          lmc[5],    digits=d),
                  f90=emp[1], f95=emp[2], f98=emp[3], f99=emp[4], f99p5=emp[5])
  logit <- log(wolfP/(1-wolfP))
  lmr <- Lmoments::Lcoefs( logit, rmax=5)
  print(lmr)
  lmc <- lmomco::lmoms.cov(logit, se="lmrse")
  print(lmc)
  df$logitmu        <- round(          lmr[1],    digits=d)
  df$logitvar       <- round((sqrt(pi)*lmr[2])^2, digits=d)
  df$logitlam2      <- round(          lmr[2],    digits=d)
  df$logittau3      <- round(          lmr[3],    digits=d)
  df$logittau4      <- round(          lmr[4],    digits=d)
  df$logittau5      <- round(          lmr[5],    digits=d)
  df$logitmuse      <- round(          lmc[1],    digits=d)
  df$logitl2se      <- round(          lmc[2],    digits=d)
  df$logitt3se      <- round(          lmc[3],    digits=d)
  df$logitt4se      <- round(          lmc[4],    digits=d)
  df$logitt5se      <- round(          lmc[5],    digits=d)
  D <- rbind(D, df)
  write.table(D, file=OUTFILE, row.names=FALSE, sep="\t")

  if(SHOW_PLOT) {
    mu <- lm(log(D$mu)~log(D$n)); vr <- lm(log(D$var)~log(D$n))
    plot(D$n, D$mu,  type="l", col="seagreen3", log="xy", xlim=range(ns), bty="n",
         xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlab="Sample sizes", ylab="")
    x <- seq(min(ns), max(n), by=1)
    y <- exp(coefficients(mu)[1] + coefficients(mu)[2]*log(x))
    x[x < 10^par()$usr[1]] <- NA; x[x > 10^par()$usr[2]] <- NA
    y[y < 10^par()$usr[3]] <- NA; y[y > 10^par()$usr[4]] <- NA
    lines(x, y, lwd=3, col="seagreen4")
    par(col.axis="black")
    axis(1, at=axTicks(1), labels=TRUE,  lwd=1, lwd.ticks=1)
    axis(3, at=axTicks(1), labels=FALSE, lwd=1, lwd.ticks=1)
    mtext("Mean", side=2, line=2, las=0, col="seagreen4")
    par(col.axis="seagreen4")
    axis(2, at=10^par()$usr[3:4], labels=FALSE, lwd=2, lwd.ticks=0, col="seagreen4")
    axis(2, at=axTicks(2),        labels=TRUE,  lwd=0, lwd.ticks=1, col="seagreen4")
    par(new = TRUE); par(mar=c(5, 5, 4, 5)+0.1)
    par(col.axis="salmon4")
    plot(D$n, D$var, type="l", col="salmon2",   log="xy", xlim=range(ns), bty="n",
         xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlab="", ylab="")
    x <- seq(min(ns), max(n), by=1)
    y <- exp(coefficients(vr)[1] + coefficients(vr)[2]*log(x))
    x[x < 10^par()$usr[1]] <- NA; x[x > 10^par()$usr[2]] <- NA
    y[y < 10^par()$usr[3]] <- NA; y[y > 10^par()$usr[4]] <- NA
    lines(x, y, lwd=3, col="salmon4")
    axis(4, at=10^par()$usr[3:4], labels=FALSE, lwd=2, lwd.ticks=0, col="salmon4")
    axis(4, at=axTicks(2),        labels=TRUE,  lwd=0, lwd.ticks=1, col="salmon4")
    mtext("Variance", side=4, line=2, las=0, col="salmon4")
  }
  message("writing to '", OUTFILE, "'")
  write.table(D, file=OUTFILE, row.names=FALSE, sep="\t")
}
