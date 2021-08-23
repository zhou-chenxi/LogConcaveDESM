#' Computation the Normalizing Constant of a Log-concave Score Matching Density Estimate
#'
#' Computes the normalizing constant of a log-concave score matching density estimate.
#'
#' @param scorematching_logconcave An object of class "LogConcaveDESM",
#' usually the output of \code{\link{lcd_scorematching}} or \code{\link{cv_optimal_density_estimate}}.
#' @param minus_const A numeric to be subtracted in the exponent to
#' ensure the finite-ness of the integration result. Default is \code{0}.
#'
#' @return The normalizing constant of a log-concave score matching density estimate.
#' @export
#'
#' @examples
#' set.seed(1119)
#' N <- 100
#' data <- rnorm(N)
#' domain <- c(-5, 5)
#' # no penalty term
#' result <- lcd_scorematching(data, domain, penalty_param = 0)
#' normalizing_const(result, minus_const = 500)
#'
normalizing_const <- function(scorematching_logconcave, minus_const = 0) {

    sorted_data <- scorematching_logconcave$sorted_unique_data
    domain <- scorematching_logconcave$domain
    all_sorted_data <- c(domain[1], sorted_data, domain[2])
    opt_theta <- scorematching_logconcave$opt_theta

    # compute (g'(a))^*
    weighted_col_AB <- sweep(
        scorematching_logconcave$matrix_A + scorematching_logconcave$matrix_B,
        MARGIN = 2,
        STATS = scorematching_logconcave$data_weights,
        FUN = '*'
    )
    g_a_opt <- -sum(scorematching_logconcave$opt_theta * rowSums(weighted_col_AB)) / 2

    result <- 0
    # integrate each sub-interval
    for (i in 1:(length(scorematching_logconcave$sorted_unique_data) + 1)) {

        a0 <- (opt_theta[i + 1] - opt_theta[i]) / (all_sorted_data[i + 1] - all_sorted_data[i])

        a <- a0 / 6

        b <- (opt_theta[i] - a0 * all_sorted_data[i]) / 2

        if (i == 1) {

            c <- a0 * all_sorted_data[i] ** 2 / 2 - opt_theta[i] * all_sorted_data[i] + g_a_opt

            d <- (opt_theta[i] * all_sorted_data[i] ** 2 / 2 - domain[1] * g_a_opt -
                      a0 * all_sorted_data[i] ** 3 / 6)

        } else {

            data_subset_i_diff <- diff(all_sorted_data[1:i])
            term2_1 <- (sum(scorematching_logconcave$opt_theta[1:(i - 1)] * data_subset_i_diff) +
                            sum(scorematching_logconcave$opt_theta[2:i] * data_subset_i_diff)) / 2

            c <- (a0 * all_sorted_data[i] ** 2 / 2 - opt_theta[i] * all_sorted_data[i] + g_a_opt +
                      term2_1)

            term2_2 <- (sum(scorematching_logconcave$opt_theta[1:(i - 1)] * data_subset_i_diff *
                                all_sorted_data[2:i]) +
                            sum(scorematching_logconcave$opt_theta[2:i] * data_subset_i_diff *
                                    all_sorted_data[2:i])) / 2
            term3 <- sum(diff(scorematching_logconcave$opt_theta[1:i]) * (data_subset_i_diff) ** 2) / 6
            term4 <- sum(scorematching_logconcave$opt_theta[1:(i - 1)] * (data_subset_i_diff) ** 2) / 2

            d <- (opt_theta[i] * all_sorted_data[i] ** 2 / 2 - domain[1] * g_a_opt -
                      a0 * all_sorted_data[i] ** 3 / 6 - term2_2 + term3 + term4)

        }

        result <- result + integrate_expcubic(
            -a, -b, -c, -d,
            lower = all_sorted_data[i], upper = all_sorted_data[i + 1],
            minus_const = minus_const
        )

    }

    return(result)

}
