setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(copBasic)
library(copula)
library(vioplot)
library(randtoolbox)

hyphen <- "\u00ad"
options(scipen=16)

n <- 10; nsim <- 20000
J <- copula::indepTestSim(n, 2, N=100000)
WW <- II <- RR <- TT <- vector(mode="numeric", nsim)
for(i in seq_len(nsim)) {
  uv <- randtoolbox::halton(n, d=2, usetime=TRUE)
  I <- W <- R <- KT <- NULL
  try(I  <-   copula::indepTest(  uv,         J                    )$pvalues[1])
  try(W  <- copBasic::wolfCOPtest(uv[,1], uv[,2],       asuv=TRUE  )$p.value[2])
  try(R  <-              cor.test(uv[,1], uv[,2], method="spearman")$p.value[1])
  try(KT <-              cor.test(uv[,1], uv[,2], method="kendall" )$p.value[1])
  WW[i] <- W; II[i] <- I; RR[i] <- R; TT[i] <- KT
}
df <- data.frame(method=c(rep("copBasic::\nwolfCOPtest", nsim),
                          rep("copula::\nindepTest",     nsim),
                          rep("Spearman\nRho",           nsim),
                          rep("Kendall\nTau",            nsim)),
                 pvalue=c(WW, II, RR, TT))
df <- df[order(df$method),]
df$method <- as.factor(df$method)


pdffile <- "halton_pitest2_vio.pdf"
pdf(pdffile, useDingbats=FALSE)
  vioplot(df$pvalue~df$method, ylim=c(0,1), range=4, ljoin=2, las=1, tcl=0.5,
          mgp=c(3, 1.2, 0), bty="n", h=0.02, areaEqual=TRUE, border=NA, cex.axis=0.95,
          col    =c("lightgreen", "lightblue", "palevioletred", "salmon1"),
          rectCol=c("seagreen4",  "lightblue4", "palevioletred4", "salmon4"),
          lineCol=c("seagreen4",  "lightblue4", "palevioletred4", "salmon4"),
          yaxt="n", yaxs="i", xaxs="i", xlab="", ylab="")
  par(las=1, mgp=c(3, 0.5, 0))
  axis(2, at=axTicks(2), labels=TRUE,  lwd=0, lwd.ticks=1, tcl=0.5); par(las=1)
  axis(4, at=axTicks(2), labels=FALSE, lwd=0, lwd.ticks=1, tcl=0.5); par(las=1)
  axis(3, at=axTicks(1), labels=FALSE, lwd=0, lwd.ticks=1, tcl=0.5); par(las=1)
  tix <- seq(0.1, 0.9, by=0.1)
  axis(2, at=tix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=0.3)
  axis(4, at=tix, labels=FALSE, lwd=0, lwd.ticks=1, tcl=0.3)

  mtext("Hypothesis test method agaist bivariate Independence (Product copula)", side=1, line=2.6)
  txt <- paste0("Distribution of p", hyphen, "values by simulation at indicated sample size")
  par(las=0); mtext(txt, side=2, line=2)

  par(las=1, mgp=c(3, 0.15, 0))
  axis(3, at=1, labels=length(unique(WW)), lwd=0, lwd.ticks=0,
          col.axis="seagreen4",      font.axis=2)
  axis(3, at=2, labels=length(unique(II)), lwd=0, lwd.ticks=0,
          col.axis="lightblue4",     font.axis=2)
  axis(3, at=3, labels=length(unique(TT)), lwd=0, lwd.ticks=0,
          col.axis="palevioletred4", font.axis=2)
  axis(3, at=4, labels=length(unique(RR)), lwd=0, lwd.ticks=0,
          col.axis="salmon4",        font.axis=2)
  txt <- paste0("Number of unique p", hyphen, "values for sample n = 10 at ",
                               prettyNum(nsim, big.mark=","), " simulations")
  mtext(txt, line=1.2, font=1, col=grey(0.3))
dev.off()
system(paste0("pdfcrop --clip --margins 5 ", pdffile, " ", pdffile))



hyphen <- "\u00ad"
pvt <- paste0("p", hyphen, "value")
seed <- 2; set.seed(2); n <- 10
J <- copula::indepTestSim(n, 2, N=100000)
    pdffile <- "halton_pitest1_quad.pdf"
pdf(pdffile, useDingbats=FALSE, width=7, height=6)
  xlab <- "U, nonexceedance probability"; ylab <- "V, nonexceedance probability"
  tix <- seq(0, 1, by=0.1)
  layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE))
  par(xpd=NA, bg="white", las=1, lend=2, mgp=c(2.5, 0.8, 0), mai=c(0.24, 0.7, 0.24, 0.06))
  uv1 <- randtoolbox::halton(n, d=2, usetime=TRUE)
  uv2 <- randtoolbox::halton(n, d=2, usetime=TRUE)
  uv3 <- randtoolbox::halton(n, d=2, usetime=TRUE)
  uv4 <- randtoolbox::halton(n, d=2, usetime=TRUE)
  i1  <- sprintf("%0.2f", round(    copula::indepTest(uv1, J)$pvalues[1],                  digits=2))
  i2  <- sprintf("%0.2f", round(    copula::indepTest(uv2, J)$pvalues[1],                  digits=2))
  i3  <- sprintf("%0.2f", round(    copula::indepTest(uv3, J)$pvalues[1],                  digits=2))
  i4  <- sprintf("%0.2f", round(    copula::indepTest(uv4, J)$pvalues[1],                  digits=2))
  r1  <- sprintf("%0.2f", round(    cor.test(uv1[,1], uv1[,2], method="spearman")$p.value, digits=2))
  r2  <- sprintf("%0.2f", round(    cor.test(uv2[,1], uv2[,2], method="spearman")$p.value, digits=2))
  r3  <- sprintf("%0.2f", round(    cor.test(uv3[,1], uv3[,2], method="spearman")$p.value, digits=2))
  r4  <- sprintf("%0.2f", round(    cor.test(uv4[,1], uv4[,2], method="spearman")$p.value, digits=2))
  t1  <- sprintf("%0.2f", round(    cor.test(uv1[,1], uv1[,2], method="kendall" )$p.value, digits=2))
  t2  <- sprintf("%0.2f", round(    cor.test(uv2[,1], uv2[,2], method="kendall" )$p.value, digits=2))
  t3  <- sprintf("%0.2f", round(    cor.test(uv3[,1], uv3[,2], method="kendall" )$p.value, digits=2))
  t4  <- sprintf("%0.2f", round(    cor.test(uv4[,1], uv4[,2], method="kendall" )$p.value, digits=2))
  w1  <- sprintf("%0.2f", round( wolfCOPtest(uv1[,1], uv1[,2])$p.value[2],                 digits=2))
  w2  <- sprintf("%0.2f", round( wolfCOPtest(uv2[,1], uv2[,2])$p.value[2],                 digits=2))
  w3  <- sprintf("%0.2f", round( wolfCOPtest(uv3[,1], uv3[,2])$p.value[2],                 digits=2))
  w4  <- sprintf("%0.2f", round( wolfCOPtest(uv4[,1], uv4[,2])$p.value[2],                 digits=2))
  plot(uv1[,1], uv1[,2], type="n", las=1, bty="n", xlim=c(0,1), ylim=c(0,1),
       xlab="", ylab=ylab, xaxs="i", yaxs="i", xaxt="n", yaxt="n")
  axis(1, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(3, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(2, at=tix, labels=TRUE,  lwd=0, lwd.ticks=1)
  axis(4, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  polygon(c(par()$usr[1], par()$usr[2], par()$usr[2], par()$usr[1], par()$usr[1]),
          c(par()$usr[3], par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), lwd=1)
  points(uv1[,1], uv1[,2], pch=21, col="#7D0112", bg="white", lwd=1.5)
  legend("topleft", c(paste0('Halton sequence "random" variate machine time to start'),
                      paste0("  sample Kendall Tau ",    pvt, " = ",  t1),
                      paste0("  sample Spearman Rho ",   pvt, " = ",  r1),
                      paste0("  copBasic::wolfCOPtest ", pvt, " = ",  w1),
                      paste0("  copula::indepTest ",     pvt, " = ",  i1)),
         bty="o", box.lty=0, bg=grey(1, 0.7), cex=0.7, y.intersp=1, inset=0.01)
  text(par()$usr[2], par()$usr[3], "(A)", adj=c(1.2,-1.1), font=2)
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
  legend("topleft", c(paste0('Halton sequence "random" variate machine time to start'),
                      paste0("  sample Kendall Tau ",    pvt, " = ",  t2),
                      paste0("  sample Spearman Rho ",   pvt, " = ",  r2),
                      paste0("  copBasic::wolfCOPtest ", pvt, " = ",  w2),
                      paste0("  copula::indepTest ",     pvt, " = ",  i2)),
         bty="o", box.lty=0, bg=grey(1, 0.7), cex=0.7, y.intersp=1, inset=0.01)
  text(par()$usr[2], par()$usr[3], "(B)", adj=c(1.2,-1.1), font=2)
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
  legend("topleft", c(paste0('Halton sequence "random" variate machine time to start'),
                      paste0("  sample Kendall Tau ",    pvt, " = ",  t3),
                      paste0("  sample Spearman Rho ",   pvt, " = ",  r3),
                      paste0("  copBasic::wolfCOPtest ", pvt, " = ",  w3),
                      paste0("  copula::indepTest ",     pvt, " = ",  i3)),
         bty="o", box.lty=0, bg=grey(1, 0.7), cex=0.7, y.intersp=1, inset=0.01)
   text(par()$usr[2], par()$usr[3], "(C)", adj=c(1.2,-1.1), font=2)
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
  points(uv4[,1], uv4[,2], pch=21, col="darkorchid2", bg="white", lwd=1.5)
  legend("topleft", c(paste0('Halton sequence "random" variate machine time to start'),
                      paste0("  sample Kendall Tau ",    pvt, " = ",  t4),
                      paste0("  sample Spearman Rho ",   pvt, " = ",  r4),
                      paste0("  copBasic::wolfCOPtest ", pvt, " = ",  w4),
                      paste0("  copula::indepTest ",     pvt, " = ",  i4)),
         bty="o", box.lty=0, bg=grey(1, 0.7), cex=0.7, y.intersp=1, inset=0.01)
   text(par()$usr[2], par()$usr[3], "(D)", adj=c(1.2,-1.1), font=2)
dev.off()
system(paste0("pdfcrop --clip --margins 5 ", pdffile, " ", pdffile))

pdffile <- "halton_pitest3_haldemo.pdf"
pdf(pdffile, useDingbats=FALSE)
  par(xpd=NA)
  uv1 <- randtoolbox::halton(200, d=2, usetime=TRUE)
  plot(uv1[,1], uv1[,2], type="n", las=1, bty="n", xlim=c(0,1), ylim=c(0,1),
       xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n")
  axis(1, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(3, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  axis(2, at=tix, labels=FALSE,  lwd=0, lwd.ticks=1)
  axis(4, at=tix, labels=FALSE, lwd=0, lwd.ticks=1)
  polygon(c(par()$usr[1], par()$usr[2], par()$usr[2], par()$usr[1], par()$usr[1]),
          c(par()$usr[3], par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), lwd=1)
  points(uv1[,1], uv1[,2], pch=21, col="grey50", bg="white", lwd=1.5)
dev.off()
system(paste0("pdfcrop --clip --margins 5 ", pdffile, " ", pdffile))
