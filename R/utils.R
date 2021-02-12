# Util functions
#------------------------------------------------------------------------------------
#' Finds all roots of a function within an interval
#'
#' Searches the interval from lower to upper for all roots of a function f with respect
#'     to its first argument.  Based on rootSolve::uniroot.all.
#' @param f The function for which the root is sought..
#' @param interval A vector containing the end-points of the interval to be
#'     searched for the root.
#' @param lower Lower end point of interval.
#' @param upper Upper end point of interval.
#' @param tol Desired accuracy of root
#' @param maxiter Maximum number of iterations
#' @param trace Integer number passed to function uniroot. If positive, tracing information is produced.
#'     Higher values giving more details.
#' @param n Number of subintervals in which the root is sought.
#' @param ...	Additional named or unnamed arguments to be passed to f.
#' @export

all_roots <- function (f, interval, lower = min(interval), upper = max(interval),
                       tol = .Machine$double.eps^0.2, maxiter = 1000, trace = 0,
                       n = 1000, ...)
{
  if (!missing(interval) && length(interval) != 2)
    stop("'Interval' must be a vector of length 2")
  if (!is.numeric(lower) || !is.numeric(upper) || lower >= upper)
    stop("lower < upper  is not fulfilled")
  xseq <- seq(lower, upper, len = n + 1)
  mod <- f(xseq, ...) #Try and remove sapply
  Equi <- xseq[mod == 0]
  len_Equi <- length(Equi)
  ss <- mod[1:n] * mod[2:(n + 1)]
  ii <- which(ss < 0)
  Equi <- c(Equi, rep(NA, length(ii)))
  if (length(ii) > 0){
    for (i in 1:length(ii)){
      Equi[len_Equi + i] <- stats::uniroot(f, lower = xseq[ii[i]], upper = xseq[ii[i]+1],
                                           maxiter = maxiter,
                                           tol = tol, trace = trace, ...)$root
    }
  }
  return(Equi)
}
