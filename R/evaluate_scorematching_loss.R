#' Evaluation of the Score Matching Loss Function
#'
#' Evaluates the score matching loss function.
#'
#' @param scorematching_logconcave An object of class "LogConcaveDESM",
#' usually the output of \code{\link{lcd_scorematching}} or \code{\link{cv_optimal_density_estimate}}.
#' @param new_data A numeric vector of real numbers at which the score matching loss function
#' should be evaluated.
#'
#' @return A numeric of the score matching loss function evaluated at the new_data.
#' @export
#'
#' @examples
#' set.seed(1119)
#' N <- 100
#' data <- rnorm(N)
#' domain <- c(-5, 5)
#' result <- lcd_scorematching(data, domain, penalty_param = 1e-10)
#'
#' evaluate_scorematching_loss(scorematching_logconcave = result, new_data = result$data)
#'
#' evaluate_scorematching_loss(scorematching_logconcave = result,
#' new_data = seq(-5, 5, by = 0.01))
#'
evaluate_scorematching_loss <- function(scorematching_logconcave, new_data) {

    # Evaluate the score matching loss function
    deriv2 <- evaluate_logden_deriv2(scorematching_logconcave, newx = new_data)$logderiv2_vals
    deriv1 <- evaluate_logden_deriv1(scorematching_logconcave, newx = new_data)$logderiv1_vals
    result <- mean(deriv1 ** 2 / 2 + deriv2)

    return(result)
}
