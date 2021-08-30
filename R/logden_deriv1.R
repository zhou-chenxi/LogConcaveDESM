#' Computation and Plotting the First Derivative of the Logarithm of the Log-concave Score Matching Density Estimate
#'
#' Based on a "LogConcaveDESM" object, evaluates and plots the first derivative of the logarithm of
#' the log-concave score matching density estimate.
#'
#' @param scorematching_logconcave An object of class "LogConcaveDESM",
#' usually the output of \code{\link{lcd_scorematching}} or \code{\link{cv_optimal_density_estimate}}.
#' @param newx A numeric vector of real numbers at which the first derivative of the logarithm of
#' the log-concave score matching density estimate should be evaluated.
#'
#' @import ggplot2
#' @return A data frame with the first column being the sorted newx and
#' the second column being the corresponding first derivative values of the logarithm of
#' the log-concave score matching density estimate.
#'
#' @export
#'
#' @examples
#' set.seed(1119)
#' N <- 100
#' data <- rnorm(N)
#' domain <- c(-5, 5)
#' result <- lcd_scorematching(data, domain, penalty_param = 1e-10)
#' # evaluationn
#' evaluate_logden_deriv1(scorematching_logconcave = result,
#' newx = seq(result$domain[1], result$domain[2], 0.01))
#'
#' # plot
#' plot_logden_deriv1(scorematching_logconcave = result,
#' plot_domain = result$domain, plot_points_cnt = 500)
#'
evaluate_logden_deriv1 <- function(scorematching_logconcave, newx) {

    domain <- scorematching_logconcave$domain
    if (min(newx) < domain[1] || max(newx) > domain[2]) {
        stop('newx is outside of the domain.')
    }

    sorted_data <- scorematching_logconcave$sorted_unique_data
    all_sorted_data <- c(domain[1], sorted_data, domain[2])
    newx <- sort(newx)
    opt_theta <- scorematching_logconcave$opt_theta

    # compute (g'(a))^*
    weighted_col_AB <- sweep(
        scorematching_logconcave$matrix_A + scorematching_logconcave$matrix_B,
        MARGIN = 2,
        STATS = scorematching_logconcave$data_weights,
        FUN = '*'
    )
    g_a_opt <- -sum(opt_theta * rowSums(weighted_col_AB)) / 2

    # integral of second derivative part
    interval_member <- data.frame(
        newx = newx,
        interval = base::findInterval(
            x = newx,
            vec = all_sorted_data,
            rightmost.closed = TRUE,
            left.open = TRUE))

    logderiv1_vals <- vector()

    for (j in unique(interval_member$interval)) {

        # summation terms involving the previous intervals
        if (j == 1) {
            val1 <- 0
        } else {

            data_subset_j_diff <- diff(all_sorted_data[1:j])
            val1 <- (sum(scorematching_logconcave$opt_theta[1:(j - 1)] * data_subset_j_diff) +
                         sum(scorematching_logconcave$opt_theta[2:j] * data_subset_j_diff)) / 2

        }

        # terms involving newx - Xi
        theta_i <- opt_theta[j]
        theta_i1 <- opt_theta[j + 1]
        Xi <- all_sorted_data[j]
        Xi1 <- all_sorted_data[j + 1]

        newx_subset <- interval_member[interval_member$interval == j, 'newx']

        val2 <- (0.5 * (theta_i1 - theta_i) / (Xi1 - Xi) * (newx_subset - Xi) ** 2 +
                     theta_i * (newx_subset - Xi))

        logderiv1_vals <- append(logderiv1_vals, val1 + val2)

    }

    result <- data.frame(
        newx_sorted = newx,
        logderiv1_vals = -(g_a_opt + logderiv1_vals)
    )

    return(result)
}

#' @rdname evaluate_logden_deriv1
#' @param plot_domain A numeric vector to indicate the domain of the plot.
#' @param plot_points_cnt A numeric to indicate the number of points for evaluating and plotting.
#' Default is \code{100}.
#'
#' @return A ggplot2 plot of the first derivative of the logarithm of
#' the log-concave score matching density estimate over the specified plot domain.
#' @export
#'
plot_logden_deriv1 <- function(scorematching_logconcave, plot_domain, plot_points_cnt = 100) {

    stopifnot(length(plot_domain) == 2)
    if (plot_domain[1] > plot_domain[2]) {
        plot_domain <- sort(plot_domain)
    } else if (plot_domain[1] == plot_domain[2]) {
        stop('The two points in plot_domain are identical. Please provide a meaningful plot_domain.')
    }

    plot_df <- evaluate_logden_deriv1(
        scorematching_logconcave = scorematching_logconcave,
        newx = seq(plot_domain[1], plot_domain[2], length.out = plot_points_cnt)
    )

    plot <- ggplot2::ggplot() +
        ggplot2::geom_line(
            data = plot_df,
            ggplot2::aes(x = newx_sorted, y = logderiv1_vals),
            size = 0.8,
            color = 'red') +
        ggplot2::geom_rug(
            data = data.frame(original_data = scorematching_logconcave$data),
            ggplot2::aes(x = original_data),
            color = 'black',
            alpha = 0.5
        ) +
        ggplot2::labs(
            x = 'x',
            y = 'first derivative',
            title = 'First Derivative of Log-density Estimate') +
        ggplot2::theme_bw()

    plot <- plot + ggplot2::coord_cartesian(xlim = plot_domain)

    return(plot)

}
