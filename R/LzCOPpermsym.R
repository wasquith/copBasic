"LzCOPpermsym" <-
function(cop=NULL, para=NULL, n=5E4,
         empfast=FALSE, empfunc=function(x) x,
         type=c("halton", "sobol", "torus", "runif"),
         as.abs=TRUE, as.vec=FALSE, as.mat=FALSE, plot=FALSE, ...) {

  if(empfast) {
    if(ncol(para) != 2) {
      warning("empfast is true, so para needs to have 2 columns as cbind(u,v)")
      return(NULL)
    }
    if(is.null(empfunc))       empfunc <- function(x) x
    if(! is.function(empfunc)) empfunc <- function(x) x
    ifelse(as.abs,  absf <- abs,  absf <- function(x) x)
    if(as.mat) {
      return( 3 * absf( EMPIRgrid_fast(para[,1:2], para=para, gridonly=TRUE, ...) -
                        EMPIRgrid_fast(para[,2:1], para=para, gridonly=TRUE, ...) ) )
    }
        H <-  3 * absf( EMPIRgrid_fast(para[,1:2], para=para, gridonly=TRUE, ...) -
                        EMPIRgrid_fast(para[,2:1], para=para, gridonly=TRUE, ...) )
        H <- as.vector(H); names(H) <- NULL
      return( empfunc( H ) )
  }

  type <- match.arg(type)
  n <- as.integer( n ) # to silence certain upstream warnings()

  if(type != "runif") {
    if(! "randtoolbox" %in% installed.packages()) {
      warning("randtoolbox package not installed, switch type to 'runif'")
      type <- 'runif'
    }
  }

  if(type == "runif") {
    ruv <- matrix(runif(2*n),    ncol=2)
  } else if(type == "halton") {
    ruv <- randtoolbox::halton(n, dim=2)
  } else if(type == "sobol" ) {
    ruv <- randtoolbox::sobol(n,  dim=2)
  } else if(type == "torus" ) {
    ruv <- randtoolbox::torus(n,  dim=2)
  }

  if(plot) {
    plot(ruv, main=type,
         xlab="U, NONEXCEEDANCE PROBABILITY", ylab="V, NONEXCEEDANCE")
  }

  if(as.vec | as.mat) as.abs <- FALSE

  ifelse(as.abs, absf <- abs, absf <- function(x) x)

  sym <- absf( COP(ruv[,1], ruv[,2], cop=cop, para=para, ...) -
               COP(ruv[,2], ruv[,1], cop=cop, para=para, ...) )
  if(as.mat) {
    return( cbind(ruv, 3 *        sym ) )
  } else {
    if(as.vec) {
      return(          3 *        sym   )
    } else {
      if(as.abs) {
        return(        3 * max(   sym ) )
      } else {
        return(        3 * range( sym ) )
      }
    }
  }
}


#"ofunc" <- function(p, cop=NA, rho=NA, permsym=NA, uv=NA) {
#      p <- pnorm(p)
#   trho <- rhoCOP(     cop=cop, para=p)
#   rmse <- rmseCOP(uv, cop=cop, para=p)
#   if(is.na(permsym)) return((rho - trho)^2 + rmse^2)
#   tsym <- mean(LzCOPpermsym2(cop=cop, para=p, as.abs=TRUE, as.vec=TRUE))
#   return((rho - trho)^2 + (permsym - tsym)^2 + rmse^2)
#}

#uv <- rCOP(n=1000, cop=MOcop, para=c(0.9,0.2))
#M <- mean(LzCOPpermsym2(cop=EMPIRcop, para=uv, as.abs=TRUE, as.vec=TRUE))
#R <- rhoCOP(as.sample=TRUE, para=uv)
#E <- rmseCOP(uv, cop=EMPIRcop, para=uv)
#init.para <- c(0,0)
#rt <- optim(init.para, ofunc, cop=MOcop, rho=R, permsym=NA, uv=uv)
#fitt.para <- pnorm(rt$par)
#uv2 <- simCOP(n=100, cop=MOcop, para=fitt.para)
#R2  <- rhoCOP(cop=MOcop, para=fitt.para)
#M2  <- mean(LzCOPpermsym2(cop=MOcop, para=fitt.para, as.abs=TRUE, as.vec=TRUE))
