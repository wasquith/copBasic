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
plot(Z$n, Z$logittau3, log="x", cex=log10(Z$nsim)-2.9, col=1)
#plot(Z$n, Z$logittau4, log="x", cex=log10(Z$nsim)-2.9, col=1)
#stop()


plotlmrdia(lmrdia(), xlim=c(0.05,0.35), ylim=c(0.1,0.20), empty=TRUE,
           xaxs="i", yaxs="i", las=1)
points(Z$tau3,      Z$tau4,      pch=21, cex=0.8, col="turquoise4", bg="turquoise")
plotlmrdia(lmrdia(usrtrim=TRUE), add=TRUE, nopoints=TRUE, autolegend=TRUE, xleg="topleft",
           noaep4=TRUE, nogpa=TRUE, nogov=TRUE, noglo=TRUE, nopdq3=TRUE,
           nolimits=TRUE, lwd.cex=2, expand.names=TRUE)
nevels <- seq(3,14, by=1)
for(i in seq_len(length(nevels))) {
  wnt <- nevels[i] == Z$n
  with(Z[wnt,], points(weighted.mean(tau3, nsim), weighted.mean(tau4, nsim), pch=16, cex=0.5))
  with(Z[wnt,], text(  weighted.mean(tau3, nsim), weighted.mean(tau4, nsim), nevels[i],
                     pos=1, offset=0.3, cex=0.7, col="grey10", font=2))
}


plotlmrdia(lmrdia(), xlim=c(0.05,0.35), ylim=c(0.1,0.20), empty=TRUE,
           xaxs="i", yaxs="i", las=1)
points(Z$logittau3, Z$logittau4, pch=21, cex=0.8, col="salmon4",    bg="salmon1"  )

xy <- NULL
for(n in 3:100) {
  xy <- rbind(xy, data.frame(tau3=wolfCOPtest(0, n)$lmoms_logit_sigma[3],
                             tau4=wolfCOPtest(0, n)$lmoms_logit_sigma[4]))
}
row.names(xy) <- NULL
lines( xy[,1], xy[,2], col="wheat2", lwd=2)
points(xy[,1], xy[,2], col="wheat4", pch=16, cex=0.3)

plotlmrdia(lmrdia(usrtrim=TRUE), add=TRUE, nopoints=TRUE, autolegend=TRUE, xleg="topleft",
           noaep4=TRUE, nogev=TRUE, nogpa=TRUE, nogov=TRUE, noglo=TRUE, nopdq3=TRUE,
           nolimits=TRUE, lwd.cex=2, expand.names=TRUE)
gno <- lmrdia()$gno; pe3 <- lmrdia()$pe3
x <- pe3[,1]; y <- (gno[,2]+pe3[,2])/2
x[x < par()$usr[1] | x > par()$usr[2]] <- NA
y[y < par()$usr[1] | y > par()$usr[2]] <- NA
lines(x, y, lty=3)
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

nevels <- seq(200,2000, by=100)
for(i in seq_len(length(nevels)-1)) {
  wnt <- nevels[i] < Z$n & Z$n < nevels[i+1]
  txt <- paste0(nevels[i], "-", nevels[i+1])
  with(Z[wnt,], points(weighted.mean(logittau3, nsim), weighted.mean(logittau4, nsim), pch=16, cex=0.5))
  #with(Z[wnt,], text(  mean(logittau3), mean(logittau4), txt,
  #                   pos=3, offset=0.3, cex=0.7, col="grey10"))
}
