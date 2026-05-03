setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(copBasic)
library(lmomco)

#n_at_minLam1 <- optim(1000, function(n) wolfCOPtest(0, n)$lmoms_logit_sigma[1])$par
#n_at_minLam2 <- optim(1000, function(n) wolfCOPtest(0, n)$lmoms_logit_sigma[2])$par
n_at_minTau3 <- optim(1000, function(n) wolfCOPtest(0, n)$lmoms_logit_sigma[3])$par
#n_at_minTau4 <- optim(1000, function(n) wolfCOPtest(0, n)$lmoms_logit_sigma[4])$par
n_at_minTau3 <- 10000

D <- NULL
D <- read.table("aa_lapsims.txt", sep="\t", header=TRUE)
files <- sort(list.files(pattern="^mc_wolfPI"))
for(f in list.files(pattern="^mc_wolfPI")) {
  D <- rbind(D, read.table(f, header=TRUE))
}
D <- D[order(D$n),]
D <- D[complete.cases(D),]
write.table(D, file="aa_allsims.txt", sep="\t", row.names=FALSE, quote=FALSE)

h <- aggregate(D, by=list(D$n), length)
print(paste0(h$Group.1[h$nsim < 3], collapse=","))

Z <- NULL # now compute weighted mean columns
for(k in sort(unique(D$n))) {
  y <- D[D$n == k,]
  z <- y[1,]
  #for(i in 3:ncol(y)) {
  #  z[,i] <- weighted.mean(y[,i], y$nsim)
  #}
  strs <- c("mu", "var", "lam2", "tau3", "tau4", "tau5")
  vtrs <- c("muse", "l2se", "l2se", "t3se", "t4se", "t5se")
  for(i in seq_len(length(strs))) {
    z[,strs[i]] <- weighted.mean(y[,strs[i]], 1/y[,vtrs[i]])
  }
  strs <- c("logitmu",   "logitvar",  "logitlam2", "logittau3", "logittau4", "logittau5")
  vtrs <- c("logitmuse", "logitl2se", "logitl2se", "logitt3se", "logitt4se", "logitt5se")
  for(i in seq_len(length(strs))) {
    z[,strs[i]] <- weighted.mean(y[,strs[i]], 1/y[,vtrs[i]])
  }
  z[,1] <- sum(y$nsim)
  Z <- rbind(Z, z)
}
Z <- Z[is.finite(Z$logitmu),]
Z$wgts <- sqrt(Z$nsim); Z$wgts <- Z$wgts/sum(Z$wgts) * length(Z$wgts)

#col <- 2+as.numeric(Z$n %/% 2 == Z$n / 2)
plot(Z$n, Z$logitmu,   log="x", cex=log10(Z$nsim)-2.9, col=1); #stop()
plot(Z$n, Z$logitlam2, log="x", cex=log10(Z$nsim)-2.9, col=1); #stop()
#stop()
plot(Z$n, Z$logittau3, log="x", cex=log10(Z$nsim)-2.9, col=1, ylim=c(0.1,0.35))
plot(Z$n, Z$logittau4, log="x", cex=log10(Z$nsim)-2.9, col=1, ylim=c(0.1,0.20))
#stop()

AUX <- read.table("wolfCOPaux_A.txt", sep="\t", header=TRUE)
AUX$wgts <- sqrt(AUX$nsim); AUX$wgts <- AUX$wgts/sum(AUX$wgts) * length(AUX$wgts)
AUX <- AUX[order(AUX$n),]
LMR <- aggregate(AUX, by=list(AUX$n), mean); jjj <- aggregate(AUX, by=list(AUX$n), sum )
LMR$nsim <- jjj$nsim; rm(jjj); LMR$Group.1 <- NULL # compute means and total up sims by sample size
AUX <- LMR
AUX$wgts <- sqrt(AUX$nsim); AUX$wgts <- AUX$wgts/sum(AUX$wgts) * length(AUX$wgts)
AUX <- AUX[AUX$n >= 5000,]

AUX <- merge(Z, AUX, all=TRUE)
Z <- AUX
Z <- Z[order(Z$n, decreasing=TRUE),]

plot(Z$n, Z$logitmu,   log="x", cex=log10(Z$nsim)-2.9, col=1); #stop()
plot(Z$n, Z$logitlam2, log="x", cex=log10(Z$nsim)-2.9, col=1); #stop()
#stop()
plot(Z$n, Z$logittau3, log="x", cex=log10(Z$nsim)-2.9, col=1, ylim=c(0.1,0.35))
plot(Z$n, Z$logittau4, log="x", cex=log10(Z$nsim)-2.9, col=1, ylim=c(0.1,0.20))
#stop()

plotlmrdia(lmrdia(), xlim=c(0.1,0.35), ylim=c(0.1,0.20), empty=TRUE, autoaxes=FALSE,
           xaxs="i", yaxs="i", lwd.cex=1.3)
par(las=1)
xtix <- c(0.07, seq(0.10, 0.35, by=0.05) )
ytix <- c(0.11, seq(0.12, 0.20, by=0.02) )
axis(1, at=xtix, labels=sprintf("%0.2f", xtix), lwd=0, lwd.ticks=1)
axis(3, at=xtix, labels=FALSE, lwd=0, lwd.ticks=1)
axis(2, at=ytix, labels=sprintf("%0.2f", ytix), lwd=0, lwd.ticks=1)
axis(4, at=ytix, labels=FALSE, lwd=0, lwd.ticks=1)
xtix <- c(      seq(0.11, 0.34, by=0.01) )
ytix <- c(      seq(0.11, 0.19, by=0.01) )
xtix <- as.numeric( xtix[grep("\\d[05]$", sprintf("%0.2f", xtix), invert=TRUE)] )
axis(1, at=xtix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.3)
axis(3, at=xtix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.3)
axis(2, at=ytix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.3)
axis(4, at=ytix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.3)

points(Z$tau3,      Z$tau4,      pch=21, cex=log10(Z$nsim)-2.9, col="turquoise4", bg="turquoise")
plotlmrdia(lmrdia(usrtrim=TRUE), add=TRUE, nopoints=TRUE, autolegend=TRUE, xleg="topleft",
           noaep4=TRUE, nogev=TRUE, nogpa=TRUE, nogov=TRUE, noglo=TRUE, nopdq3=TRUE,
           nolimits=TRUE, lwd.cex=2, expand.names=TRUE, legendcex=0.7)
nevels <- seq(3,14, by=1)
for(i in seq_len(length(nevels))) {
  wnt <- nevels[i] == Z$n
  with(Z[wnt,], points(weighted.mean(tau3, nsim), weighted.mean(tau4, nsim), pch=16, cex=0.5))
  with(Z[wnt,], text(  weighted.mean(tau3, nsim), weighted.mean(tau4, nsim), nevels[i],
                     pos=1, offset=0.3, cex=0.7, col="grey10", font=2))
}


pdf("wolfCOPlogitTau34.pdf", width=7, height=7, useDingbats=FALSE)
  plotlmrdia(lmrdia(), xlim=c(0.10,0.32), ylim=c(0.11,0.20), empty=TRUE, autoaxes=FALSE,
             xaxs="i", yaxs="i", lwd.cex=1.3)
  txt <- paste0("L-moment Ratio Diagram of Logit Transformation of Schweizer-Wolff Sigma\n",
                "from Monte Carlo Simulation Using the Independence Copula")
  mtext(txt, line=1, cex=0.9, font=2)
  par(las=1, xpd=NA)
  xtix <- c(0.07, seq(0.10, 0.35, by=0.05) ); xtix <- xtix[par()$usr[1] <= xtix & xtix <= par()$usr[2]]
  ytix <- c(0.11, seq(0.12, 0.20, by=0.02) ); ytix <- ytix[par()$usr[1] <= ytix & ytix <= par()$usr[2]]
  axis(1, at=xtix, labels=sprintf("%0.2f", xtix), lwd=0, lwd.ticks=1)
  axis(3, at=xtix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(2, at=ytix, labels=sprintf("%0.2f", ytix), lwd=0, lwd.ticks=1)
  axis(4, at=ytix, labels=FALSE, lwd=0, lwd.ticks=1)
  xtix <- c(      seq(0.10, 0.34, by=0.01) ); xtix <- xtix[par()$usr[1] <= xtix & xtix <= par()$usr[2]]
  ytix <- c(      seq(0.11, 0.19, by=0.01) ); ytix <- ytix[par()$usr[1] <= ytix & ytix <= par()$usr[2]]
  xtix <- as.numeric( xtix[grep("\\d[05]$", sprintf("%0.2f", xtix), invert=TRUE)] )
  axis(1, at=xtix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.3)
  axis(3, at=xtix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.3)
  axis(2, at=ytix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.3)
  axis(4, at=ytix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=-0.3)
  points(Z$logittau3, Z$logittau4, pch=21, cex=log10(Z$nsim)-3, col="salmon4",    bg="salmon1", lwd=0.7  )
  xy <- NULL
  ns <- sort(unique(c(10^seq(0,1,0.001), 10^seq(0, 4, by=0.02))))
  for(n in ns) {
    xy <- rbind(xy, data.frame(tau3=wolfCOPtest(0, n)$lmoms_logit_sigma[3],
                               tau4=wolfCOPtest(0, n)$lmoms_logit_sigma[4]))
  }
  row.names(xy) <- NULL
  xy <- xy[par()$usr[1] <= xy[,1] & xy[,1] <= par()$usr[2],]
  xy <- xy[par()$usr[3] <= xy[,2] & xy[,2] <= par()$usr[4],]
  lines( xy[,1], xy[,2], col="wheat2", lwd=2)
  xy <- NULL
  for(n in c(3:14, 17, 19, 25)) {
    xy <- rbind(xy, data.frame(tau3=wolfCOPtest(0, n)$lmoms_logit_sigma[3],
                               tau4=wolfCOPtest(0, n)$lmoms_logit_sigma[4]))
  }
  row.names(xy) <- NULL
  xy <- xy[par()$usr[1] <= xy[,1] & xy[,1] <= par()$usr[2],]
  xy <- xy[par()$usr[3] <= xy[,2] & xy[,2] <= par()$usr[4],]
  points(xy[,1], xy[,2], col="wheat4", pch=17, cex=0.8)

  txt <- c("Predicted Tau3 and Tau4 logits\n  by regressions within wolfCOPtest(...)",
           "Predictions for the near plotting sample\n  sizes (n = 3-14 plotted)",
           "Weighted mean for a sample size with\n  symbol size scaling to log10 count",
           "Mean for sample size or range of sample\n  sizes (not all labeled, by 100s up to 700)",
           "Weighted mean for samples size\n  within the range 3,000-10,000")
  par(lheight=0.85)
  legend("bottomright", txt,
         box.lty=0, inset=0.01, cex=0.7, y.intersp=1.6, adj=c(0, 0.8),
         lty=c(1, NA, NA, NA, NA), lwd=c(2, NA, NA, NA, NA),
         col=c("wheat2", "wheat4", "salmon4", "black", "black"),
         pch=c(NA, 17, 21, 16, 21), pt.lwd=c(1, 1, 1, 1, 1.21),
         pt.bg=c(NA, NA, "salmon1", NA, "white"), pt.cex=c(NA, 0.8, 1, 0.5, 1.8))
  par(lheight=1)

  plotlmrdia(lmrdia(usrtrim=TRUE), add=TRUE, nopoints=TRUE, autolegend=TRUE, xleg="topleft",
             noaep4=TRUE, nogev=TRUE, nogpa=TRUE, nogov=TRUE, noglo=TRUE, nopdq3=TRUE,
             nolimits=TRUE, lwd.cex=2, expand.names=TRUE, inset=0.01, legendcex=0.8)

  gno <- lmrdia()$gno; pe3 <- lmrdia()$pe3
  x <- pe3[,1]; y <- (gno[,2]+pe3[,2])/2; suppressWarnings( rm(gno, pe3))
  x[x < par()$usr[1] | x > par()$usr[2]] <- NA
  y[y < par()$usr[1] | y > par()$usr[2]] <- NA
  lines(x, y, lty=4, col="deepskyblue3")
  text(0.292, 0.1725, "Halfway between the\ntwo distributions", srt=40, cex=0.9, col="deepskyblue3")

  nevels <- seq(3000, 30000, by=27000)
  for(i in seq_len(length(nevels)-1)) {
    wnt <- nevels[i] <= Z$n-1 & Z$n <= nevels[i+1]+1
    txt <- paste0(nevels[i], "-", nevels[i+1])
    with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim),
                         pch=21, cex=1.2, col="black", bg="white", lwd=1.21))
    #with(Z[wnt,], text(  mean(logittau3), mean(logittau4), txt,
    #                   pos=3, offset=0.3, cex=0.7, col="grey10"))
  }

  nevels <- seq(7,14, by=1)
  for(i in seq_len(length(nevels))) {
    wnt <- nevels[i] == Z$n
    with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
    with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), nevels[i],
                       pos=1, offset=0.26, cex=0.7, col="grey95", font=2))
    with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), nevels[i],
                     pos=1, offset=0.3, cex=0.7, col="grey10", font=2))
  }

  nevels <- seq(16,20, by=2)
  for(i in seq_len(length(nevels)-1)) {
    wnt <- nevels[i] <= Z$n & Z$n <= nevels[i+1]
    txt <- paste0(nevels[i],   "-",  nevels[i+1])
    with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
    with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                       adj=c(1.15, 0.38), offset=0.3, cex=0.7, col="grey95", font=2, srt=90))
    with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                       adj=c(1.15, 0.5), offset=0.3, cex=0.7, col="grey10", font=2, srt=90))
  }

  nevels <- seq(20,60, by=10)
  for(i in seq_len(length(nevels)-1)) {
    wnt <- nevels[i] <= Z$n & Z$n <= nevels[i+1]
    txt <- paste0(nevels[i], "-", nevels[i+1])
    with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
    with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                       adj=c(1.15, 0.38), offset=0.3, cex=0.7, col="grey95", font=2, srt=90))
    with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                       adj=c(1.15, 0.5), offset=0.3, cex=0.7, col="grey10", font=2, srt=90))
  }

  nevels <- seq(35,45, by=10)
  for(i in seq_len(length(nevels)-1)) {
    wnt <- nevels[i] <= Z$n & Z$n <= nevels[i+1]
    txt <- paste0(nevels[i], "-", nevels[i+1])
    with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
    with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                       adj=c(1.15, 0.38), offset=0.3, cex=0.7, col="grey95", font=2, srt=90))
    with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                       adj=c(1.15, 0.5), offset=0.3, cex=0.7, col="grey10", font=2, srt=90))
  }

  nevels <- seq(60,100, by=40)
  for(i in seq_len(length(nevels)-1)) {
    wnt <- nevels[i] <= Z$n & Z$n <= nevels[i+1]
    txt <- paste0(nevels[i], "-", nevels[i+1])
    with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
    with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                       adj=c(1.15, 0.38), offset=0.3, cex=0.7, col="grey95", font=2, srt=90))
    with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                      adj=c(1.15, 0.5), offset=0.3, cex=0.7, col="grey10", font=2, srt=90))
  }

  nevels <- seq(100,300, by=100)
  for(i in seq_len(length(nevels)-1)) {
    wnt <- nevels[i] <= Z$n & Z$n <= nevels[i+1]
    txt <- paste0(nevels[i], "-", nevels[i+1])
    with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
    with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                       adj=c(1.15, 0.38), offset=0.3, cex=0.7, col="grey95", font=2, srt=90))
    with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                       adj=c(1.15, 0.5), offset=0.3, cex=0.7, col="grey10", font=2, srt=90))
  }

  nevels <- seq(300,700, by=100)
  for(i in seq_len(length(nevels)-1)) {
    wnt <- nevels[i] <= Z$n & Z$n <= nevels[i+1]
    txt <- paste0(nevels[i], "-", nevels[i+1])
    with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
    #with(Z[wnt,], text(  mean(logittau3), mean(logittau4), txt,
    #                   pos=3, offset=0.3, cex=0.7, col="grey10"))
  }

dev.off()
par(xpd=FALSE)
