#' Computation and Plotting the Log-concave Score Matching Density Estimate
#'
#' Based on a "LogConcaveDESM" object, evaluates and plots the log-concave score matching
#' density estimate.
#'
#' @param scorematching_logconcave An object of class "LogConcaveDESM",
#' usually the output of \code{\link{lcd_scorematching}} or \code{\link{cv_optimal_density_estimate}}.
#' @param newx A numeric vector of real numbers at which the log-concave score matching
#' density estimate should be evaluated.
#' @param minus_const A numeric to be subtracted in the exponent to
#' ensure the finite-ness of the integration result. Default is \code{0}.
#'
#' @import ggplot2
#' @import stats
#'
#' @return A data frame with the first column being the sorted newx and
#' the second column being the corresponding values of
#' the log-concave score matching density estimate.
#'
#' @export
#'
#' @examples
#' set.seed(1119)
#' N <- 100
#' data <- rnorm(N)
#' domain <- c(-5, 5)
#' result <- lcd_scorematching(data, domain, penalty_param = 1e-1)
#'
#' # evaluation
#' evaluate_density(
#' scorematching_logconcave = result,
#' newx = seq(result$domain[1], result$domain[2], by = 0.01))
#'
#' # plot
#' plot_density(result, minus_const = 0, plot_domain = result$domain, plot_points_cnt = 500,
#' plot_hist = TRUE)
#'
#'
evaluate_density <- function(scorematching_logconcave, newx, minus_const = 0) {

    newx <- sort(newx)

    domain <- scorematching_logconcave$domain

    valid_newx <- newx[newx >= domain[1] & newx <= domain[2]]
    invalid_newx <- newx[newx < domain[1] | newx > domain[2]]

    # compute density estimates for invalid_newx
    invalid_newx_den <- data.frame(
        newx_sorted = invalid_newx,
        density_vals = rep(0, length(invalid_newx))
    )

    # compute density estimates for valid_newx
    log_den_vals <- evaluate_logdensity(
        scorematching_logconcave = scorematching_logconcave,
        newx = valid_newx)

    den_vals <- exp(log_den_vals$logdensity_vals - minus_const)

    norm_const <- normalizing_const(scorematching_logconcave, minus_const)

    valid_newx_den <- data.frame(
        newx_sorted = log_den_vals$newx_sorted,
        density_vals = den_vals / norm_const
    )

    result <- rbind(invalid_newx_den, valid_newx_den)
    result <- result[order(result$newx_sorted), ]

    return(result)

}

#' @rdname evaluate_density
#' @param plot_domain A numeric vector to indicate the domain of the plot.
#' @param plot_points_cnt A numeric to indicate the number of points for evaluating and plotting.
#' Default is \code{100}.
#' @param plot_hist A logical value to indicate whether to the histogram. Default is \code{FALSE}.
#' @param hist_binwidth A numeric to indicate the bin width parameter in plotting the histogram;
#' only works when \code{plot_hist == TRUE}. Default is \code{NULL}, under which condition
#' the Freedman-Diaconis rule is used.
#' @param hist_alpha A numeric to indicate the opacity in plotting the histogram;
#' must range from 0 to 1, inclusively; only works when \code{plot_hist == TRUE}.
#' Default is \code{0.2}.
#'
#' @return A ggplot2 plot of the log-concave score matching density estimate
#' over the specified plot domain.
#' @export
#'
plot_density <- function(scorematching_logconcave, minus_const = 0, plot_domain, plot_points_cnt = 100,
                         plot_hist = FALSE, hist_binwidth = NULL, hist_alpha = 0.2) {

    stopifnot(length(plot_domain) == 2)
    if (plot_domain[1] > plot_domain[2]) {
        plot_domain <- sort(plot_domain)
    } else if (plot_domain[1] == plot_domain[2]) {
        stop('The two points in plot_domain are identical. Please provide a meaningful plot_domain.')
    }

    plot_df <- evaluate_density(
        scorematching_logconcave = scorematching_logconcave,
        newx = seq(plot_domain[1], plot_domain[2], length.out = plot_points_cnt),
        minus_const = minus_const
    )

    plot <- ggplot2::ggplot() +
        ggplot2::geom_line(
            data = plot_df,
            ggplot2::aes(x = newx_sorted, y = density_vals),
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
            y = 'density',
            title = 'Density Estimate') +
        ggplot2::theme_bw()

    if (plot_hist) {

        # use the Freedman-Diaconis rule for the binwidth
        if (is.null(hist_binwidth)) {
            NN <- length(scorematching_logconcave$data)
            hist_binwidth <- 2 * stats::IQR(scorematching_logconcave$data) / (NN ** (1/3))
        }
        plot <- plot + ggplot2::geom_histogram(
            data = data.frame(x = scorematching_logconcave$data),
            ggplot2::aes(x = x, y = ..density..),
            binwidth = hist_binwidth,
            color = "darkblue",
            fill = "lightblue",
            alpha = hist_alpha)

    }

    plot <- plot + ggplot2::coord_cartesian(xlim = plot_domain)

    return(plot)

}
