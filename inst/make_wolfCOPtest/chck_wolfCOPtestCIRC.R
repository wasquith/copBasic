setwd(dirname(rstudioapi::getSourceEditorContext()$path))

stop("SAFE STOP AT BEGINNING")

library(copBasic)
library( copula )
library(parallel)

# Here is a cross check of sorts between wolfCOPtest and the Independence test
# available in the copula package. We use a Circular copula with various bivariate
# associations and will see for these circumstances that wolfCOPtest out performs
# copula::indepTest by just a small amount by rejecting the NULL more often.
OUFILE <- "zz_chck_wolfCOPtestCIRCsim10kC.txt"
NSIM   <- as.integer( 150000 ); ALPHA <- 0.05
myCOPTXT <- "CIRCcop"
myCOP    <-  CIRCcop
myPARA   <- NULL
myRTdigits <- 6
myCORES    <- 4
sample_sizes <- c(68, 72, 63, seq(5, 20, by=1), seq(22, 50, by=2), seq(55, 100, by=5))
sample_sizes <- sample(sample_sizes, length(sample_sizes), replace=FALSE)
# Randomly order the sample sizes by a hunch that it helps how the cluster might
# partition work. Reorder the data frame d is trivial.

"IWfunc" <- function(n) {
  J <- copula::indepTestSim(n, 2, N=100000); WW <- II <- RR <- TT <- NULL
  for(i in seq_len(NSIM)) {
    uv <- copBasic::rCOP(n, cop=myCOP, para=myPARA); I <- W <- R <- KT <- NULL
    try(I  <-   copula::indepTest(  uv,         J                    )$pvalues[1])
    try(W  <- copBasic::wolfCOPtest(uv[,1], uv[,2],       asuv=TRUE  )$p.value[1])
    try(R  <-              cor.test(uv[,1], uv[,2], method="spearman")$p.value[1])
    try(KT <-              cor.test(uv[,1], uv[,2], method="kendall" )$p.value[1])
     if(is.null(I)) I <- NA; if(is.null(W)) W <- NA; if(is.null(R)) R <- NA; if(is.null(KT)) KT <- NA
    II <- c(II, I);         WW <- c(WW, W);         RR <- c(RR, R);         TT <- c(TT, KT)
  }
  II <- II[! is.na(II)]; WW <- WW[! is.na(WW)]; RR <- RR[! is.na(RR)]; TT <- TT[! is.na(TT)]
                         Irejrate <-                sum(II < ALPHA) / length(II)
                                        Wrejrate <- sum(WW < ALPHA) / length(WW)
                                        Rrejrate <- sum(RR < ALPHA) / length(RR)
                                        Trejrate <- sum(TT < ALPHA) / length(TT)
  return( paste0(n, ":", Irejrate, ":", Wrejrate, ":", Rrejrate, ":", Trejrate) )
}

D <- NULL; apnd <- FALSE
   cl <- parallel::makeCluster(getOption("cl.cores", myCORES))
        parallel::clusterExport( cl, c("ALPHA", "NSIM", "myCOP", "myPARA"))
   R <- parallel::parSapply(     cl, sample_sizes, IWfunc)
        parallel::stopCluster(   cl   )
  # -----------------------------------------------------------------------
  d <- NULL
  for(i in seq_len(length(R)) ) {
    r <- as.numeric( unlist( strsplit(R[i], ":") ) )
    d <- rbind(d, data.frame(copula=myCOPTXT, nsim=NSIM, n=as.integer(r[1]), alpha=ALPHA,
               rho_given=NA, tau_given=NA, para=NA, rho=0, tau=0,
               Irejrate=r[2], Wrejrate=r[3], Rrejrate=r[4], Trejrate=r[5]))
  }
  d <- d[order(d$n),]; print( d )
  write.table(d, file=OUFILE, sep="\t", row.names=FALSE, quote=FALSE, append=apnd, col.names=! apnd)
  D <- rbind( D, d); rownames(D) <- NULL
  print(Sys.time()); apnd <- TRUE

message("Done")

write.table(D, file=OUFILE, sep="\t", row.names=FALSE, quote=FALSE)

stop()


ALPHA <- 0.05
D <- read.table("zz_chck_wolfCOPtestCIRCsim10kB.txt", header=TRUE)
D <- rbind(D, read.table("zz_chck_wolfCOPtestCIRCsim10kC.txt", header=TRUE))

H <- NULL
  for(n in sort(unique(D$n))) {
    G <- D[D$n == n,]
    H <- rbind(H, data.frame(copula="CIRCcop", nsim=sum(G$nsim), n=G$n[1], alpha=G$alpha[1],
                             rho_given=G$rho_given[1], tau_given=G$tau_given[1], para=G$para[1],
                             rho=weighted.mean(G$rho, w=G$nsim), tau=weighted.mean(G$tau, w=G$nsim),
                             Irejrate=weighted.mean(G$Irejrate, w=G$nsim, na.rm=TRUE),
                             Wrejrate=weighted.mean(G$Wrejrate, w=G$nsim, na.rm=TRUE),
                             Rrejrate=weighted.mean(G$Rrejrate, w=G$nsim, na.rm=TRUE),
                             Trejrate=weighted.mean(G$Trejrate, w=G$nsim, na.rm=TRUE)))
  }
D <- H



D$Irejrate <- round(D$Irejrate, digits=6)
D$Wrejrate <- round(D$Wrejrate, digits=6)
D$Rrejrate <- round(D$Rrejrate, digits=6)
D$Trejrate <- round(D$Trejrate, digits=6)




#pdf("zz_chck_wolfCOPtestCIRC_plot.pdf", useDingbats=FALSE, width=7, height=6)
#  par(xpd=NA, bg="white", las=1, lend=2, mgp=c(2.5, 0.8, 0))
#  txt <- paste0("Null hypothesis rejection rate (or 50x) by method at alpha = ", ALPHA)
#  plot(range(D$n), c(0,1+ALPHA), log="x", type="n", las=1, bty="n",
#       xlab="", ylab=txt, xaxs="i", yaxs="i", xaxt="n", yaxt="n")
#  mtext("Sample size (logarithmic scale)", side=1, line=1.8)
#  polygon(10^c(par()$usr[1], par()$usr[2], par()$usr[2], par()$usr[1], par()$usr[1]),
#             c(par()$usr[3], par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]),
#             col="white", border=NA)
#  par(lend=1); lines(10^par()$usr[1:2], rep(ALPHA, 2), lwd=2, col="grey50"); par(lend=2)
#  axis(2, at=axTicks(2), labels=TRUE,  lwd=0, lwd.ticks=1)
#  axis(4, at=axTicks(2), labels=FALSE, lwd=0, lwd.ticks=1)
#  tix <- c(5, 6, 7, 8, 9, 12, 14, 16, 18, 20, 30, 40, 50, 60, 70, 80, 90)
#  axis(1, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
#  axis(3, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
#  tix <- c(10, 100)
#  axis(1, at=tix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.9)
#  axis(3, at=tix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.9)
#  axis(1, at=c(5, 10, 20, 50, 100), labels=TRUE, lwd=0, lwd.ticks=0)
#  tix <- seq(0.1, 0.9, by=0.2)
#  axis(2, at=tix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.2)
#  axis(4, at=tix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.2)
#  lines(D$n, D$Irejrate,    col="#A55A0B", lwd=2.0, lty=2)
#  lines(D$n, D$Wrejrate,    col="#0070A3", lwd=2.0, lty=1)
#  lines(D$n, 50*D$Trejrate, col="#00A0A9", lwd=1.2, lty=4)
#  lines(D$n, 50*D$Rrejrate, col="red2",    lwd=1.2, lty=4)
#  lxt <- c("Significance level alpha = 0.05",
#           "copBasic::wolfCOPtest on Circular copula (copBasic::CIRCcop)",
#           "copula::indepTest on Circular copula (copBasic::CIRCcop)",
#           '50 x rate Mann-Kendall test, stats::cor.test(method="kendall") on Circular copula',
#           '50 x rate Spearman Rho test, stats::cor.test(method="spearman") on Circular copula')
#  legend("topleft", lxt, cex=0.7, bty="o", box.lty=0, bg=NA, lwd=c(2, 2, 2, 1.2, 1.2), seg.len=3.5,
#         col=c("grey50", "#0070A3", "#A55A0B", "#00A0A9", "red2"), lty=c(1, 1, 2, 4, 4))
#  polygon(10^c(par()$usr[1], par()$usr[2], par()$usr[2], par()$usr[1], par()$usr[1]),
#             c(par()$usr[3], par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), lwd=1)
#dev.off()

hyphen <- "\u00ad"
pvt <- paste0("p", hyphen, "value")
seed <- 2; set.seed(2); n <- 50
J <- copula::indepTestSim(n, 2, N=100000)

pdf("zz_chck_wolfCOPtestCLCIRC_plot.pdf", useDingbats=FALSE, width=7, height=6)
  xlab <- "U, nonexceedance probability"; ylab <- "V, nonexceedance probability"
  tix <- seq(0, 1, by=0.1)
  layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE))
  par(xpd=NA, bg="white", las=1, lend=2, mgp=c(2.5, 0.8, 0), mai=c(0.24, 0.7, 0.24, 0.06))
  uv1 <- simCOP(n, cop=CLcop,   para=CLcop(tau=0.00)$para, seed=seed, graphics=FALSE)
  uv2 <- simCOP(n, cop=CLcop,   para=CLcop(tau=0.24)$para, seed=seed, graphics=FALSE)
  uv3 <- simCOP(n, cop=CLcop,   para=CLcop(tau=0.48)$para, seed=seed, graphics=FALSE)
  uv4 <- simCOP(n, cop=CIRCcop, para=NULL,                 seed=seed, graphics=FALSE)
  i1  <- sprintf("%0.2f", round( copula::indepTest(uv1, J)$pvalues[1], digits=2))
  i2  <- sprintf("%0.2f", round( copula::indepTest(uv2, J)$pvalues[1], digits=2))
  i3  <- sprintf("%0.2f", round( copula::indepTest(uv3, J)$pvalues[1], digits=2))
  i4  <- sprintf("%0.2f", round( copula::indepTest(uv4, J)$pvalues[1], digits=2))
  r1  <- sprintf("%0.2f", round( cor.test(uv1[,1], uv1[,2], method="spearman")$p.value, digits=2))
  r2  <- sprintf("%0.2f", round( cor.test(uv2[,1], uv2[,2], method="spearman")$p.value, digits=2))
  r3  <- sprintf("%0.2f", round( cor.test(uv3[,1], uv3[,2], method="spearman")$p.value, digits=2))
  r4  <- sprintf("%0.2f", round( cor.test(uv4[,1], uv4[,2], method="spearman")$p.value, digits=2))
  t1  <- sprintf("%0.2f", round( cor.test(uv1[,1], uv1[,2], method="kendall" )$p.value, digits=2))
  t2  <- sprintf("%0.2f", round( cor.test(uv2[,1], uv2[,2], method="kendall" )$p.value, digits=2))
  t3  <- sprintf("%0.2f", round( cor.test(uv3[,1], uv3[,2], method="kendall" )$p.value, digits=2))
  t4  <- sprintf("%0.2f", round( cor.test(uv4[,1], uv4[,2], method="kendall" )$p.value, digits=2))
  w1  <- sprintf("%0.2f", round( wolfCOPtest(uv1[,1], uv1[,2])$p.value[1], digits=2))
  w2  <- sprintf("%0.2f", round( wolfCOPtest(uv2[,1], uv2[,2])$p.value[1], digits=2))
  w3  <- sprintf("%0.2f", round( wolfCOPtest(uv3[,1], uv3[,2])$p.value[1], digits=2))
  w4  <- sprintf("%0.2f", round( wolfCOPtest(uv4[,1], uv4[,2])$p.value[1], digits=2))
  plot(uv1[,1], uv1[,2], type="n", las=1, bty="n", xlim=c(0,1), ylim=c(0,1),
       xlab="", ylab=ylab, xaxs="i", yaxs="i", xaxt="n", yaxt="n")
  axis(1, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(3, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(2, at=tix, labels=TRUE,  lwd=0, lwd.ticks=1)
  axis(4, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  polygon(c(par()$usr[1], par()$usr[2], par()$usr[2], par()$usr[1], par()$usr[1]),
          c(par()$usr[3], par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), lwd=1)
  points(uv1[,1], uv1[,2], pch=21, col="#7D0112", bg="white", lwd=1.5)
  legend("topleft", c(paste0("Clayton copula, given Kendall Tau = 0"),
                      paste0("  sample Kendall Tau ",    pvt, " = ",  t1),
                      paste0("  sample Spearman Rho ",   pvt, " = ",  r1),
                      paste0("  copBasic::wolfCOPtest ", pvt, " = ",  w1),
                      paste0("  copula::indepTest ",     pvt, " = ",  i1)),
         bty="o", box.lty=0, bg=grey(1, 0.5), cex=0.7, y.intersp=1, inset=0.01)

  par(xpd=NA, bg="white", las=1, lend=2, mgp=c(2.5, 0.8, 0), mai=c(0.24, 0.4, 0.24, 0.34))
  plot(uv2[,1], uv2[,2], type="n", las=1, bty="n", xlim=c(0,1), ylim=c(0,1),
       xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n")
  axis(1, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(3, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(2, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(4, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  polygon(c(par()$usr[1], par()$usr[2], par()$usr[2], par()$usr[1], par()$usr[1]),
          c(par()$usr[3], par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), lwd=1)
  points(uv2[,1], uv2[,2], pch=21, col="#78CDB3", bg="white", lwd=1.5)
  legend("topleft", c(paste0("Clayton copula, given Kendall Tau = 0.24"),
                      paste0("  sample Kendall Tau ",    pvt, " = ",  t2),
                      paste0("  sample Spearman Rho ",   pvt, " = ",  r2),
                      paste0("  copBasic::wolfCOPtest ", pvt, " = ",  w2),
                      paste0("  copula::indepTest ",     pvt, " = ",  i2)),
         bty="o", box.lty=0, bg=grey(1, 0.5), cex=0.7, y.intersp=1, inset=0.01)

  par(xpd=NA, bg="white", las=1, lend=2, mgp=c(2.5, 0.8, 0), mai=c(0.5, 0.7, 0.14, 0.06))
  plot(uv3[,1], uv3[,2], type="n", las=1, bty="n", xlim=c(0,1), ylim=c(0,1),
       xlab="", ylab=ylab, xaxs="i", yaxs="i", xaxt="n", yaxt="n")
  axis(1, at=tix, labels=TRUE,  lwd=0, lwd.ticks=1)
  axis(3, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(2, at=tix, labels=TRUE,  lwd=0, lwd.ticks=1)
  axis(4, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  mtext(xlab, side=1, line=1.8, cex=0.83)
  polygon(c(par()$usr[1], par()$usr[2], par()$usr[2], par()$usr[1], par()$usr[1]),
          c(par()$usr[3], par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), lwd=1)
  points(uv3[,1], uv3[,2], pch=21, col="#1F28A2", bg="white", lwd=1.5)
  legend("topleft", c(paste0("Clayton copula, given Kendall Tau = 0.48"),
                      paste0("  sample Kendall Tau ",    pvt, " = ",  t3),
                      paste0("  sample Spearman Rho ",   pvt, " = ",  r3),
                      paste0("  copBasic::wolfCOPtest ", pvt, " = ",  w3),
                      paste0("  copula::indepTest ",     pvt, " = ",  i3)),
         bty="o", box.lty=0, bg=grey(1, 0.5), cex=0.7, y.intersp=1, inset=0.01)

  par(xpd=NA, bg="white", las=1, lend=2, mgp=c(2.5, 0.8, 0), mai=c(0.5, 0.4, 0.14, 0.34))
  plot(uv4[,1], uv4[,2], type="n", las=1, bty="n", xlim=c(0,1), ylim=c(0,1),
       xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n")
  mtext(xlab, side=1, line=1.8, cex=0.83)
  axis(1, at=tix, labels=TRUE,  lwd=0, lwd.ticks=1)
  axis(3, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(2, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(4, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  polygon(c(par()$usr[1], par()$usr[2], par()$usr[2], par()$usr[1], par()$usr[1]),
          c(par()$usr[3], par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), lwd=1)
  points(uv4[,1], uv4[,2], pch=21, col="#0070A3", bg="white", lwd=1.5)
  legend("topleft", c(paste0("Circular copula, theoretical Kendall Tau = 0"),
                      paste0("  sample Kendall Tau ",    pvt, " = ",  t4),
                      paste0("  sample Spearman Rho ",   pvt, " = ",  r4),
                      paste0("  copBasic::wolfCOPtest ", pvt, " = ",  w4),
                      paste0("  copula::indepTest ",     pvt, " = ",  i4)),
         bty="o", box.lty=0, bg=grey(1, 0.5), cex=0.7, y.intersp=1, inset=0.01)
dev.off()



pdf("zz_chck_wolfCOPtestCIRC_plot.pdf", useDingbats=FALSE, width=7, height=6)

  par(xpd=NA, bg="white", las=1, lend=2, mgp=c(2.5, 0.8, 0))
  txt <- paste0("Null hypothesis rejection rate (or 50x) by method at alpha = ", ALPHA)
  plot(range(D$n), c(0,1+ALPHA), log="x", type="n", las=1, bty="n",
       xlab="", ylab=txt, xaxs="i", yaxs="i", xaxt="n", yaxt="n")
  mtext("Sample size (logarithmic scale)", side=1, line=1.8)
  polygon(10^c(par()$usr[1], par()$usr[2], par()$usr[2], par()$usr[1], par()$usr[1]),
             c(par()$usr[3], par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]),
             col="white", border=NA)
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
  lines(D$n, D$Irejrate,    col="#A55A0B", lwd=2.0, lty=2)
  lines(D$n, D$Wrejrate,    col="#0070A3", lwd=2.0, lty=1)
  lines(D$n, 50*D$Trejrate, col="#00A0A9", lwd=1.2, lty=4)
  lines(D$n, 50*D$Rrejrate, col="red2",    lwd=1.2, lty=4)
  lxt <- c("Significance level alpha = 0.05",
           "copBasic::wolfCOPtest on Circular copula (copBasic::CIRCcop)",
           "copula::indepTest on Circular copula (copBasic::CIRCcop)",
           '50 x rate Mann-Kendall test, stats::cor.test(method="kendall") on Circular copula',
           '50 x rate Spearman Rho test, stats::cor.test(method="spearman") on Circular copula')
  legend("topleft", lxt, cex=0.7, bty="o", box.lty=0, bg=NA, lwd=c(2, 2, 2, 1.2, 1.2), seg.len=3.5,
         col=c("grey50", "#0070A3", "#A55A0B", "#00A0A9", "red2"), lty=c(1, 1, 2, 4, 4))
  polygon(10^c(par()$usr[1], par()$usr[2], par()$usr[2], par()$usr[1], par()$usr[1]),
             c(par()$usr[3], par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), lwd=1)
dev.off()



system("pdfcrop --margins 5 zz_chck_wolfCOPtestCLCIRC_plot.pdf zz_chck_wolfCOPtestCLCIRC_plot.pdf")

system("pdfcrop --margins 5 zz_chck_wolfCOPtestCIRC_plot.pdf zz_chck_wolfCOPtestCIRC_plot.pdf")

