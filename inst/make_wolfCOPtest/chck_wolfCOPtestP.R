setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(copBasic)
library(lmomco)

files <- sort(list.files(pattern="^mc_wolfPI"))
D <- read.table("aa_lapsims.txt",  header=TRUE)
for(f in list.files(pattern="^mc_wolfPI")) {
  D <- rbind(D, read.table(f, header=TRUE))
}
D <- D[order(D$n),]
D <- D[complete.cases(D),]
write.table(D, file="aa_macmini.txt", sep="\t", row.names=FALSE, quote=FALSE)

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
#plot(Z$n, Z$logittau4, log="x", cex=log10(Z$nsim)-2.9, col=1, ylim=c(0.1,0.20))
#stop()


plotlmrdia(lmrdia(), xlim=c(0.1,0.35), ylim=c(0.1,0.20), empty=TRUE, autoaxes=FALSE,
           xaxs="i", yaxs="i", lwd.cex=1.3)
par(las=1)#11768000
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
           nolimits=TRUE, lwd.cex=2, expand.names=TRUE)
nevels <- seq(3,14, by=1)
for(i in seq_len(length(nevels))) {
  wnt <- nevels[i] == Z$n
  with(Z[wnt,], points(weighted.mean(tau3, nsim), weighted.mean(tau4, nsim), pch=16, cex=0.5))
  with(Z[wnt,], text(  weighted.mean(tau3, nsim), weighted.mean(tau4, nsim), nevels[i],
                     pos=1, offset=0.3, cex=0.7, col="grey10", font=2))
}


pdf("wolfCOPlogitTau34.pdf", useDingbats=FALSE)
plotlmrdia(lmrdia(), xlim=c(0.1,0.35), ylim=c(0.11,0.20), empty=TRUE, autoaxes=FALSE,
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
points(Z$logittau3, Z$logittau4, pch=21, cex=log10(Z$nsim)-2.9, col="salmon4",    bg="salmon1"  )
xy <- NULL
for(n in c(3:100, seq(100,3000, by=100))) {
  xy <- rbind(xy, data.frame(tau3=wolfCOPtest(0, n)$lmoms_logit_sigma[3],
                             tau4=wolfCOPtest(0, n)$lmoms_logit_sigma[4]))
}
row.names(xy) <- NULL
lines( xy[,1], xy[,2], col="wheat2", lwd=2)
xy <- NULL
for(n in c(3:14)) {
  xy <- rbind(xy, data.frame(tau3=wolfCOPtest(0, n)$lmoms_logit_sigma[3],
                             tau4=wolfCOPtest(0, n)$lmoms_logit_sigma[4]))
}
row.names(xy) <- NULL
points(xy[,1], xy[,2], col="wheat4", pch=17, cex=0.8)

txt <- c("Predicted logits of Tau3 and Tau4\nby regressions in wolfCOPtest()",
         "Predictions for near plotting sample sizes\n(n 3-14 plotted)",
         "Weighted mean for a sample size with\nsymbol size scaling to log10 count",
         "Mean for sample size or range of sample\nsizes (not all labeled, see source code)",
         "Weighted mean for samples size\nwithin the range 1,000-3,000")
par(lheight=0.8)
legend("bottomright", txt,
       box.lty=0, inset=0.01, cex=0.8, y.intersp=1.6, adj=c(0, 0.8),
       lty=c(1, NA, NA, NA, NA), lwd=c(2, NA, NA, NA, NA),
       col=c("wheat2", "wheat4", "salmon4", "black", "black"),
       pch=c(NA, 17, 21, 16, 21), pt.lwd=c(1, 1, 1, 1, 1.21),
       pt.bg=c(NA, NA, "salmon1", NA, "white"), pt.cex=c(NA, 0.8, 1, 0.5, 1.8))
par(lheight=1)

plotlmrdia(lmrdia(usrtrim=TRUE), add=TRUE, nopoints=TRUE, autolegend=TRUE, xleg="topleft",
           noaep4=TRUE, nogev=TRUE, nogpa=TRUE, nogov=TRUE, noglo=TRUE, nopdq3=TRUE,
           nolimits=TRUE, lwd.cex=2, expand.names=TRUE, inset=0.01, legendcex=0.8)
gno <- lmrdia()$gno; pe3 <- lmrdia()$pe3
x <- pe3[,1]; y <- (gno[,2]+pe3[,2])/2
x[x < par()$usr[1] | x > par()$usr[2]] <- NA
y[y < par()$usr[1] | y > par()$usr[2]] <- NA
lines(x, y, lty=3)



nevels <- seq(1000, 3000, by=2000)
for(i in seq_len(length(nevels)-1)) {
  wnt <- nevels[i] <= Z$n & Z$n <= nevels[i+1]
  txt <- paste0(nevels[i], "-", nevels[i+1])
  with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim),
                       pch=21, cex=1.8, col="black", bg="white", lwd=1.21))
  #with(Z[wnt,], text(  mean(logittau3), mean(logittau4), txt,
  #                   pos=3, offset=0.3, cex=0.7, col="grey10"))
}


nevels <- seq(7,14, by=1)
for(i in seq_len(length(nevels))) {
  wnt <- nevels[i] == Z$n
  with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
  with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), nevels[i],
                     pos=1, offset=0.3, cex=0.7, col="grey10", font=2))
}

nevels <- seq(16,20, by=2)
for(i in seq_len(length(nevels)-1)) {
  wnt <- nevels[i] <= Z$n & Z$n <= nevels[i+1]
  txt <- paste0(nevels[i], "-", nevels[i+1])
  with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
  with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                     adj=c(1.15, 0.5), offset=0.3, cex=0.7, col="grey10", font=2, srt=90))
}

nevels <- seq(20,60, by=10)
for(i in seq_len(length(nevels)-1)) {
  wnt <- nevels[i] <= Z$n & Z$n <= nevels[i+1]
  txt <- paste0(nevels[i], "-", nevels[i+1])
  with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
  with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                     adj=c(1.15, 0.5), offset=0.3, cex=0.7, col="grey10", font=2, srt=90))
}


nevels <- seq(35,45, by=10)
for(i in seq_len(length(nevels)-1)) {
  wnt <- nevels[i] <= Z$n & Z$n <= nevels[i+1]
  txt <- paste0(nevels[i], "-", nevels[i+1])
  with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
  with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                     adj=c(1.15, 0.5), offset=0.3, cex=0.7, col="grey10", font=2, srt=90))
}

nevels <- seq(60,100, by=40)
for(i in seq_len(length(nevels)-1)) {
  wnt <- nevels[i] <= Z$n & Z$n <= nevels[i+1]
  txt <- paste0(nevels[i], "-", nevels[i+1])
  with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
  with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                     adj=c(1.15, 0.5), offset=0.3, cex=0.7, col="grey10", font=2, srt=90))
}

nevels <- seq(100,200, by=100)
for(i in seq_len(length(nevels)-1)) {
  wnt <- nevels[i] <= Z$n & Z$n <= nevels[i+1]
  txt <- paste0(nevels[i], "-", nevels[i+1])
  with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
  with(Z[wnt,], text(  weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), txt,
                     adj=c(1.15, 0.5), offset=0.3, cex=0.7, col="grey10", font=2, srt=90))
}

nevels <- seq(200,1000, by=100)
for(i in seq_len(length(nevels)-1)) {
  wnt <- nevels[i] <= Z$n & Z$n <= nevels[i+1]
  txt <- paste0(nevels[i], "-", nevels[i+1])
  with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
  #with(Z[wnt,], text(  mean(logittau3), mean(logittau4), txt,
  #                   pos=3, offset=0.3, cex=0.7, col="grey10"))
}
text(0.305, 0.1775, "Halfway between the\ntwo distributions", srt=44, cex=0.9)
dev.off()
