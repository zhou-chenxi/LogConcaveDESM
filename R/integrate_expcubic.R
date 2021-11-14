#' Integration of an Exponential-Cubic Function
#'
#' Integrates a function of the form \eqn{\exp(a*x^3+b*x^2+c*x+d)} over an interval.
#'
#' @param a,b,c,d The coefficients of \eqn{x^3}, \eqn{x^2}, \eqn{x}, and the constant term, respectively.
#' @param lower,upper Numeric values to specify the limits of integration.
#' @param minus_const A numeric to be subtracted in the exponent to
#' ensure the finite-ness of the integration result. Default is \code{0}.
#'
#' @return A numeric of the integral of \eqn{\exp(a*x^3+b*x^2+c*x+d)} from \code{lower} to \code{upper}.
#' @export
#'
#' @examples
#' integrate_expcubic(a = 1, b = -2, c = 3, d = -2, lower = -2, upper = 2)
#'
integrate_expcubic <- function(a, b, c, d, lower, upper, minus_const = 0) {

    fun <- function(x) {
        exp(a * x ** 3 + b * x ** 2 + c * x + d - minus_const)
    }

    result <- integrate(fun, lower, upper)$value

    return(result)

}
