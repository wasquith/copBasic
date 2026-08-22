setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(copBasic)
library( copula )
library(parallel)

# Here is a cross check of sorts between wolfCOPtest and the Independence test
# available in the copula package. We use a Clayton copula with various bivariate
# associations and will see for these circumstances that wolfCOPtest out performs
# copula::indepTest by just a small amount by rejecting the NULL more often.
OUFILE <- "zz_chck_wolfCOPtestCLsim10kF.txt"
NSIM   <- as.integer( 10000 ); ALPHA <- 0.05
myCOPTXT <- "CLcop"
myCOP    <-  CLcop
myTAUS   <- rev( seq(0, 0.5, by=0.02) )
myPARAf  <- function(tau=NA) (2*tau) / (1-tau)
myPARA   <- NULL
myRTdigits <- 6
myCORES    <- 4
sample_sizes <- c(seq(5, 20, by=1), seq(22, 50, by=2), seq(55, 100, by=5))
sample_sizes <- sample(sample_sizes, length(sample_sizes), replace=FALSE)
# Randomly order the sample sizes by a hunch that it helps how the cluster might
# partition work. Reorder the data frame d is trivial.

"IWfunc" <- function(n) {
  J <- copula::indepTestSim(n, 2, N=100000); WW <- II <- NULL
  for(i in seq_len(NSIM)) {
    uv <- copBasic::rCOP(n, cop=myCOP, para=myPARA); I <- W <- NULL
    try(I <-   copula::indepTest(  uv,         J             )$pvalues[1])
    try(W <- copBasic::wolfCOPtest(uv[,1], uv[,2], asuV=TRUE )$p.value[1])
     if(is.null(I)) I <- NA; if(is.null(W)) W <- NA
    II <- c(II, I);         WW <- c(WW, W)
  }
  II <- II[! is.na(II)]; WW <- WW[! is.na(WW)]
                         Irejrate <-                sum(II < ALPHA) / length(II)
                                        Wrejrate <- sum(WW < ALPHA) / length(WW)
  return( paste0(n, ":", Irejrate, ":", Wrejrate) )
}

D <- NULL; apnd <- FALSE
for(tau in myTAUS) {
  message(paste0(Sys.time(), "; Tau=", tau));             myPARA <- myPARAf( tau )
  suppressWarnings( rhoi <- round( rhoCOP(cop=myCOP, para=myPARA), digits=myRTdigits ) )
  suppressWarnings( taui <- round( tauCOP(cop=myCOP, para=myPARA), digits=myRTdigits ) )
  # -----------------------------------------------------------------------
  cl <- parallel::makeCluster(getOption("cl.cores", myCORES))
        parallel::clusterExport( cl, c("ALPHA", "NSIM", "myCOP", "myPARA"))
   R <- parallel::parSapply(     cl, sample_sizes, IWfunc)
        parallel::stopCluster(   cl   )
  # -----------------------------------------------------------------------
  d <- NULL
  for(i in seq_len(length(R)) ) {
    r <- as.numeric( unlist( strsplit(R[i], ":") ) )
    d <- rbind(d, data.frame(copula=myCOPTXT, nsim=NSIM, n=as.integer(r[1]), alpha=ALPHA,
       rho_given=NA, tau_given=tau, para=myPARA, rho=rhoi, tau=taui, Irejrate=r[2], Wrejrate=r[3]))
  }
  d <- d[order(d$n),]; print( d )
  write.table(d, file=OUFILE, sep="\t", row.names=FALSE, quote=FALSE, append=apnd, col.names=! apnd)
  D <- rbind( D, d)
  print(Sys.time()); apnd <- TRUE
}
message("Done")

write.table(D, file=OUFILE, sep="\t", row.names=FALSE, quote=FALSE)

stop()

ALPHA <- 0.05
D <- read.table("zz_chck_wolfCOPtestCLsim50kA.txt", header=TRUE)
D <- rbind(D, read.table("zz_chck_wolfCOPtestCLsim10kB.txt", header=TRUE))
D <- rbind(D, read.table("zz_chck_wolfCOPtestCLsim10kC.txt", header=TRUE))
D <- rbind(D, read.table("zz_chck_wolfCOPtestCLsim10kD.txt", header=TRUE))
D <- rbind(D, read.table("zz_chck_wolfCOPtestCLsim10kE.txt", header=TRUE))
D <- rbind(D, read.table("zz_chck_wolfCOPtestCLsim10kF.txt", header=TRUE))


H <- NULL
for(t in sort(unique(D$tau_given))) {
  J <- D[D$tau_given == t,]
  for(n in sort(unique(J$n))) {
    G <- J[J$n == n,]
    H <- rbind(H, data.frame(copula="CLcop", nsim=sum(G$nsim), n=G$n[1], alpha=G$alpha[1],
                             rho_given=G$rho_given[1], tau_given=G$tau_given[1], para=G$para[1],
                             rho=weighted.mean(G$rho, w=G$nsim), tau=weighted.mean(G$tau, w=G$nsim),
                             Irejrate=weighted.mean(G$Irejrate, w=G$nsim),
                             Wrejrate=weighted.mean(G$Wrejrate, w=G$nsim)))
  }
}
D <- H
D$Irejrate <- round(D$Irejrate, digits=6)
D$Wrejrate <- round(D$Wrejrate, digits=6)

write.table(D, file="zz_chck_wolfCOPtestCL_data.txt", sep="\t", row.names=FALSE, quote=FALSE)

D <- D[as.character(D$tau_given) %in% as.character(rev( seq(0, 0.5, by=0.04) )),]
cols <- hcl.colors(length(unique(D$tau_given))+5, palette="Roma")
mc   <- floor(median(1:length(cols))); cols <- cols[-c(mc-2, mc-1, mc, mc+1, mc+2)]
cols <- data.frame(tau_given=sort(unique(D$tau_given)), col=cols)

topcol <- cols$col[cols$tau_given == "0.48"]

pdf("zz_chck_wolfCOPtestCL_plots.pdf", useDingbats=FALSE, width=7, height=6)
  par(xpd=NA, bg="white", las=1, lend=2, mgp=c(2.5, 0.8, 0))
  txt <- paste0("Rejection rate by method at alpha = ", ALPHA)
  plot(range(D$n), c(0,1+ALPHA), log="x", type="n", las=1, bty="n",
       xlab="", ylab=txt, xaxs="i", yaxs="i", xaxt="n", yaxt="n")
  mtext("Sample size (logarithmic scale)", side=1, line=1.8)
  polygon(10^c(par()$usr[1], par()$usr[2], par()$usr[2], par()$usr[1], par()$usr[1]),
             c(par()$usr[3], par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]),
          col="white", border=NA)
  px <- c(5, 10, 10, 5, 5); py <- c(0, 0, 0.2, 0.2, 0)
  polygon(px, py, col="grey95", border="grey30", lwd=0.4)
  par(lend=1); lines(10^par()$usr[1:2], rep(ALPHA, 2), lwd=2, col="grey50"); par(lend=2)
  axis(2, at=axTicks(2), labels=TRUE,  lwd=0, lwd.ticks=1)
  axis(4, at=axTicks(2), labels=FALSE, lwd=0, lwd.ticks=1)
  tix <- c(5, 6, 7, 8, 9, 12, 14, 16, 18, 20, 30, 40, 50, 60, 70, 80, 90)
  axis(1, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(3, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  tix <- c(10, 100)
  axis(1, at=tix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.9)
  axis(3, at=tix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.9)
  axis(1, at=c(5, 10, 20, 50, 100), labels=TRUE, lwd=0, lwd.ticks=0)
  tix <- seq(0.1, 0.9, by=0.2)
  axis(2, at=tix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.2)
  axis(4, at=tix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.2)
  for(t in sort(unique(D$tau_given))) {
    J <- D[D$tau_given == t,]; col <- cols$col[cols$tau_given == t]
    lines(J$n, J$Irejrate, col=col, lty=2); lines(J$n, J$Wrejrate, col=col)
    if(t < 0.30) {
      w <- 70; ix <- seq_len(nrow(J));
      wnt <- ix[J$n == w];   xa <- J$n[wnt]; ya <- (J$Irejrate[wnt] + J$Wrejrate[wnt]) / 2
      wnt <- ix[J$n == w]+1; xb <- J$n[wnt]; yb <- (J$Irejrate[wnt] + J$Wrejrate[wnt]) / 2
    } else {
      w <- 18; ix <- seq_len(nrow(J));
      wnt <- ix[J$n == w];   xa <- J$n[wnt]; ya <- (J$Irejrate[wnt] + J$Wrejrate[wnt]) / 2
      wnt <- ix[J$n == w]+1; xb <- J$n[wnt]; yb <- (J$Irejrate[wnt] + J$Wrejrate[wnt]) / 2
    }
    a <- (par()$fin[2]-sum(par()$mai[c(1,3)])) /
         (par()$fin[1]-sum(par()$mai[c(2,4)])) # inches y / inches x  (aspect ratio in real world)
    dx <- (log10(xb)-log10(xa)) / diff(par()$usr[1:2])
    dy <- (      yb -      ya ) / diff(par()$usr[3:4])
    srt <- (180/pi) * atan(a*dy/dx)
    text(w, ya, paste0("T = ", t), col="white", cex=0.7, adj=c(0.5, 0.4), font=2, srt=srt)
    text(w, ya, paste0("T = ", t), col="white", cex=0.7, adj=c(0.5, 0.6), font=2, srt=srt)
    text(w, ya, paste0("T = ", t), col=col,     cex=0.7, adj=c(0.5, 0.5), font=2, srt=srt)
  }
  lxt <- c("Significance level alpha = 0.05",
           "copBasic::wolfCOPtest (colored by Kendall's tau [T])",
           "copula::indepTest (colored by Kendall's tau [T])")
  legend("topleft", lxt, cex=0.7, bty="o", box.lty=0, bg=NA, lwd=c(2, 1, 1),
                         col=c("grey50", topcol, topcol),    lty=c(1, 1, 2))
  polygon(10^c(par()$usr[1], par()$usr[2], par()$usr[2], par()$usr[1], par()$usr[1]),
             c(par()$usr[3], par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), lwd=1)
  text(10^0.849485, 0.025, "Zoom region for next plot", adj=c(0.5, 0.5), cex=0.7)


  par(xpd=NA, lend=2)
  x <- seq(5, 10, by=0.001)
  txt <- paste0("Rejection rate by method at alpha = ", ALPHA)
  plot(range(D$n), c(0,1+ALPHA), log="", type="n", las=1, bty="n", xlim=range(x), ylim=c(0, 0.2),
       xlab="", ylab=txt, xaxs="i", yaxs="i", xaxt="n", yaxt="n")
  mtext("Sample size (arithmetic scale)", side=1, line=1.8)
  polygon(px, py, col="grey95", border=NA)
  par(lend=1); lines(par()$usr[1:2], rep(ALPHA, 2), lwd=2, col="grey50"); par(lend=2)
  axis(2, at=axTicks(2), labels=TRUE,  lwd=0, lwd.ticks=1)
  axis(4, at=axTicks(2), labels=FALSE, lwd=0, lwd.ticks=1)
  tix <- c(5, 6, 7, 8, 9, 10)
  axis(1, at=tix, labels=TRUE,  lwd=0, lwd.ticks=1)
  axis(3, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  tix <- as.character(seq(0,0.2, by=0.01));
  tix <- tix[-grep("[05]$", tix)]; tix <- tix[-grep("\\.[12]$", tix)]
  axis(2, at=tix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.2)
  axis(4, at=tix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.2)
  for(t in rev(sort(unique(D$tau_given)))) {
    J <- D[D$tau_given == t,]; col <- cols$col[cols$tau_given == t]
    y1 <- approx(J$n, J$Irejrate, xout=x, rule=2)$y
    y2 <- approx(J$n, J$Wrejrate, xout=x, rule=2)$y
    wnt1 <- par()$usr[3] <= y1 & y1 <= par()$usr[4]
    wnt2 <- par()$usr[3] <= y2 & y2 <= par()$usr[4]
    x1 <- x2 <- x; x1 <- x1[wnt1]; x2 <- x2[wnt2]
                   y1 <- y1[wnt1]; y2 <- y2[wnt2]
    lines(x1, y1, col=col, lty=2); lines(x2, y2, col=col)
    JJ <-  J[par()$usr[1] <= J$n         & J$n         <= par()$usr[2],]
    JJ <- JJ[par()$usr[3] <= JJ$Irejrate & JJ$Irejrate <= par()$usr[4],]
    JJ <- JJ[par()$usr[3] <= JJ$Wrejrate & JJ$Wrejrate <= par()$usr[4],]
    #points(JJ$n, JJ$Irejrate, col=col, pch=24, bg=col)
    #points(JJ$n, JJ$Wrejrate, col=col, pch=25, bg=col)
    if(t == 0) {
      w <- 7; ix <- seq_len(nrow(J));
      wnt <- ix[J$n == w];   xa <- J$n[wnt]; ya <- (J$Irejrate[wnt] + J$Wrejrate[wnt]) / 2
      wnt <- ix[J$n == w]+1; xb <- J$n[wnt]; yb <- (J$Irejrate[wnt] + J$Wrejrate[wnt]) / 2
    } else if(t < 0.30) {
      w <- 9; ix <- seq_len(nrow(J));
      wnt <- ix[J$n == w];   xa <- J$n[wnt]; ya <- (J$Irejrate[wnt] + J$Wrejrate[wnt]) / 2
      wnt <- ix[J$n == w]+1; xb <- J$n[wnt]; yb <- (J$Irejrate[wnt] + J$Wrejrate[wnt]) / 2
    } else if(t <= 0.4) {
      w <- 6; ix <- seq_len(nrow(J));
      wnt <- ix[J$n == w];   xa <- J$n[wnt]; ya <- (J$Irejrate[wnt] + J$Wrejrate[wnt]) / 2
      wnt <- ix[J$n == w]+1; xb <- J$n[wnt]; yb <- (J$Irejrate[wnt] + J$Wrejrate[wnt]) / 2
    } else {
      w <- 5; ix <- seq_len(nrow(J));
      wnt <- ix[J$n == w];   xa <- J$n[wnt]; ya <- (J$Irejrate[wnt] + J$Wrejrate[wnt]) / 2
      wnt <- ix[J$n == w]+1; xb <- J$n[wnt]; yb <- (J$Irejrate[wnt] + J$Wrejrate[wnt]) / 2
    }
    a <- (par()$fin[2]-sum(par()$mai[c(1,3)])) /
         (par()$fin[1]-sum(par()$mai[c(2,4)])) # inches y / inches x  (aspect ratio in real world)
    dx <- (      xb -      xa ) / diff(par()$usr[1:2])
    dy <- (      yb -      ya ) / diff(par()$usr[3:4])
    srt <- (180/pi) * atan(a*dy/dx)
    if(t <= 0.4) {
      text(w,        ya,        paste0("T = ",t), col="white", cex=0.7, adj=c(0.5,0.4), font=2, srt=srt)
      text(w,        ya,        paste0("T = ",t), col="white", cex=0.7, adj=c(0.5,0.6), font=2, srt=srt)
      text(w,        ya,        paste0("T = ",t), col=col,     cex=0.7, adj=c(0.5,0.5), font=2, srt=srt)
    } else if(t <= 0.44) {
      text(w+0.155,  ya+0.011,  paste0("T = ",t), col="white", cex=0.7, adj=c(0.5,0.4), font=2, srt=srt)
      text(w+0.155,  ya+0.011,  paste0("T = ",t), col="white", cex=0.7, adj=c(0.5,0.6), font=2, srt=srt)
      text(w+0.155,  ya+0.011,  paste0("T = ",t), col=col,     cex=0.7, adj=c(0.5,0.5), font=2, srt=srt)
    } else {
      text(w+0.1455, ya+0.0102, paste0("T = ",t), col="white", cex=0.7, adj=c(0.5,0.4), font=2, srt=srt)
      text(w+0.1455, ya+0.0102, paste0("T = ",t), col="white", cex=0.7, adj=c(0.5,0.6), font=2, srt=srt)
      text(w+0.1455, ya+0.0102, paste0("T = ",t), col=col,     cex=0.7, adj=c(0.5,0.5), font=2, srt=srt)
    }
  }
  lxt <- c("Significance level alpha = 0.05",
           "copBasic::wolfCOPtest (colored by Kendall's tau [T])",
           "copula::indepTest (colored by Kendall's tau [T])")
  legend("bottomright", lxt, cex=0.7, bty="o", box.lty=0, bg=NA, lwd=c(2, 1, 1),
                             col=c("grey50", topcol, topcol),    lty=c(1, 1, 2))
  polygon(   c(par()$usr[1], par()$usr[2], par()$usr[2], par()$usr[1], par()$usr[1]),
             c(par()$usr[3], par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), lwd=1)
dev.off()

