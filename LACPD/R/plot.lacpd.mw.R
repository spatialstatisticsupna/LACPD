#' Classes and methods for class lacpd.mw
#'
#' Classes and methods for class lacpd.mw
#'
#' @param x an object of class lacpd.mw
#'
#'
#' @details
#' plot functions for class lacpd.mw
#'
#'
#'
#' @author Mehdi Moradi \email{m2.moradi@yahoo.com}
#' @seealso
#' \link[LACPD]{lacpd_mw}
#'
#' @export
plot.lacpd.mw <- function(x,...){

  stopifnot(any(class(x)=="lacpd.mw"))

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  Z <- attr(x,"zs")
  p <- attr(x,"ps")
  m <- attr(x,"mags")
  t <- x$s

  par(mfrow = c(1, 3), pty = "s")

  plot(t,Z,main="Z",...)
  plot(t,p,main="p.value",...)
  plot(t,m,main="magnitude",...)

}
