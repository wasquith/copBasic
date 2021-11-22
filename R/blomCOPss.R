blomCOPss <-
function(cop=NULL, para=NULL, as.sample=FALSE,
         uu=rep(0.5, 2), vv=rep(0.5, 2), trap.nan=TRUE, ...) {
   # Schmid and Schmidt doi:10.1007/s00184-006-0114-3 --- Abbreviated SS here
   if(as.sample) {
      if(is.null(para)) {
         warning("Sample Blomqvist's Beta desired but para is NULL, returning NULL")
         return(NULL)
      }
      if(length(names(para)) != 2) {
        warning("para argument must be data.frame having only two columns, returning NULL")
        return(NULL)
      }
      warning("sample estimation is not implemented")
      return(NULL)
   } else {
      if(is.null(cop)) {
         warning("must have copula argument specified, returning NULL")
         return(NULL)
      }
      if(uu <= vv && (uu > 0 | vv < 1)) {      # SS pp.4
        hdgd <- function(u, v) {               # SS description of equation 5
           n <- length(u)
          gd <- cumprod(u)[n] + cumprod(1-v)[n]
          return(list( hd=1 / ( min(u) + min(1-v) - gd), gd=gd))
        }
        hhgg <- hdgd(uu, vv)

        C    <- cop(uu[1], uu[2], para=para, ...)
        Cbar <- surfuncCOP(vv[1], vv[2], cop=cop, para=para, ...)
        if(trap.nan) { # Asquith handling as explained. We know the copBasic::PSP needs
          if(is.nan(C))    C    <- 0 # this as PSP does not trap its own NaN
          if(is.nan(Cbar)) Cbar <- 0 # interred with a PSP reflected
        }
        return(hhgg$hd * (C + Cbar - hhgg$gd))
      } else {
        warning('invalid uu and vv, must be uu <= vv & (uu > 0 & vv < 1)')
      }
   }
}
