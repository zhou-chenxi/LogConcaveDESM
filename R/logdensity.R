#' Computation and Plotting the Logarithm of the Un-normalized Log-concave Score Matching Density Estimate
#'
#' Based on a "LogConcaveDESM" object, evaluates and plots the logarithm of
#' the un-normalized log-concave score matching density estimate.
#'
#' @param scorematching_logconcave An object of class "LogConcaveDESM",
#' usually the output of \code{\link{lcd_scorematching}} or \code{\link{cv_optimal_density_estimate}}.
#' @param newx A numeric vector of real numbers at which the logarithm of
#' the lun-normalized og-concave score matching density estimate should be evaluated.
#'
#' @import ggplot2
#' @return A data frame with the first column being the sorted newx and
#' the second column being the corresponding values of the logarithm of
#' the un-normalized log-concave score matching density estimate.
#'
#' @export
#'
#' @examples
#' set.seed(1119)
#' N <- 100
#' data <- rnorm(N)
#' domain <- c(-5, 5)
#' result <- lcd_scorematching(data, domain, penalty_param = 1e-10)
#'
#' # evaluation
#' evaluate_logdensity(scorematching_logconcave = result,
#' newx = seq(result$domain[1], result$domain[2], 0.01))
#'
#' # plot
#' plot_logdensity(scorematching_logconcave = result,
#' plot_domain = result$domain, plot_points_cnt = 500)
#'
evaluate_logdensity <- function(scorematching_logconcave, newx) {

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
    g_a_opt <- -sum(scorematching_logconcave$opt_theta * rowSums(weighted_col_AB)) / 2
    term1 <- g_a_opt * (newx - domain[1])

    # the remaining terms
    interval_member <- data.frame(
        newx = newx,
        interval = base::findInterval(
            x = newx,
            vec = all_sorted_data,
            rightmost.closed = TRUE,
            left.open = TRUE))

    logdensity_vals <- vector()

    for (j in unique(interval_member$interval)) {

        if (j == 1) {
            val1 <- 0
        } else {

            # term 2
            newx_subset <- interval_member[interval_member$interval == j, 'newx']

            data_subset_j_diff <- diff(all_sorted_data[1:j])
            term2_coef <- (sum(scorematching_logconcave$opt_theta[1:(j - 1)] * data_subset_j_diff) +
                               sum(scorematching_logconcave$opt_theta[2:j] * data_subset_j_diff)) / 2
            term2_1 <- term2_coef * newx_subset
            term2_2 <- (sum(scorematching_logconcave$opt_theta[1:(j - 1)] * data_subset_j_diff *
                                all_sorted_data[2:j]) +
                            sum(scorematching_logconcave$opt_theta[2:j] * data_subset_j_diff *
                                    all_sorted_data[2:j])) / 2

            term2 <- term2_1 - term2_2

            # term 3
            term3 <- sum(diff(scorematching_logconcave$opt_theta[1:j]) * (data_subset_j_diff) ** 2) / 6

            # term 4
            term4 <- sum(scorematching_logconcave$opt_theta[1:(j - 1)] * (data_subset_j_diff) ** 2) / 2

            val1 <- term2 + term3 + term4

        }

        # terms involving newx - Xi
        theta_i <- opt_theta[j]
        theta_i1 <- opt_theta[j + 1]
        Xi <- all_sorted_data[j]
        Xi1 <- all_sorted_data[j + 1]

        newx_subset <- interval_member[interval_member$interval == j, 'newx']

        term5 <- ((theta_i1 - theta_i) / (Xi1 - Xi) * (newx_subset - Xi) ** 3 / 6 +
                      theta_i * (newx_subset - Xi) ** 2 / 2)

        logdensity_vals <- append(logdensity_vals, val1 + term5)

    }

    result <- data.frame(
        newx_sorted = newx,
        logdensity_vals = -(term1 + logdensity_vals)
    )

    return(result)
}

#' @rdname evaluate_logdensity
#' @param plot_domain A numeric vector to indicate the domain of the plot.
#' @param plot_points_cnt A numeric to indicate the number of points for evaluating and plotting.
#' Default is \code{100}.
#'
#' @return A ggplot2 plot of the logarithm of
#' the un-normalized log-concave score matching density estimate over the specified plot domain.
#' @export
#'
plot_logdensity <- function(scorematching_logconcave, plot_domain, plot_points_cnt = 100) {

    stopifnot(length(plot_domain) == 2)
    if (plot_domain[1] > plot_domain[2]) {
        plot_domain <- sort(plot_domain)
    } else if (plot_domain[1] == plot_domain[2]) {
        stop('The two points in plot_domain are identical. Please provide a meaningful plot_domain.')
    }

    plot_df <- evaluate_logdensity(
        scorematching_logconcave = scorematching_logconcave,
        newx = seq(plot_domain[1], plot_domain[2], length.out = plot_points_cnt)
    )

    plot <- ggplot2::ggplot() +
        ggplot2::geom_line(
            data = plot_df,
            ggplot2::aes(x = newx_sorted, y = logdensity_vals),
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
            y = 'log-density',
            title = 'Log-density Estimate') +
        ggplot2::theme_bw()

    plot <- plot + ggplot2::coord_cartesian(xlim = plot_domain)

    return(plot)

}
