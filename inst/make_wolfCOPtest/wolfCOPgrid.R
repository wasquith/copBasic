"wolfCOPgrid" <- function(uv, ctype="1/n", deluv=1/1000, ...) {
  ef <- EMPIRgrid(para=uv, gridonly=TRUE, ctype=ctype, deluv=deluv, ...); pf <- ef
  gu <- as.numeric(rownames(pf)); gv <- as.numeric(colnames(pf))
      for(i in seq_len(nrow(pf))) pf[i,] <- gv * gu[i]
  return(12 * (1/deluv) * mean(abs(ef - pf)) / (cumprod(dim(pf))[1] - 1))
}
