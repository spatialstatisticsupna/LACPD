#' Classes and methods for class lacpd.mw
#'
#' Classes and methods for class lacpd.mw
#'
#' @param x an object of class lacpd.mw
#'
#'
#' @details
#' print function for class lacpd.mw
#'
#'
#'
#' @author Mehdi Moradi \email{m2.moradi@yahoo.com}
#' @seealso
#' \link[LACPD]{lacpd_mw}
#' @export
print.lacpd.mw <- function(x,round=5,...){
  cat("LACPD Mann-Whitney \n");
  cat("adaptive window:", " ", paste0(x$window),"\n");
  cat("   change-point:", " ", paste0(x$cp),"\n");
  cat("     statistics:",paste0(round(x$z,round)),"\n");
  cat("      magnitude:",paste0(round(x$mag,round)),"\n");
  cat("        p.value:",paste0(round(x$p,round)),"\n");
}

