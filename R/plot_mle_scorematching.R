#' Plotting the Log-concave Maximum Likelihood and Score Matching Density Estimates
#'
#' Plots the log-concave maximum likelihood and (penalized) score matching density estimates together.
#'
#' @param scorematching_logconcave An object of class "LogConcaveDESM",
#' usually the output of \code{\link{lcd_scorematching}} or \code{\link{cv_optimal_density_estimate}}.
#' @param plot_domain A numeric vector to indicate the domain of the plot.
#' @param smoothed A logical value to indicate whether to compute the smoothed log-concave maximum
#' likelihood density estimate. Default is \code{FALSE}.
#' @param plot_log A logical value to indicate whether to plot the density estimate in the log scale.
#' Default is \code{FALSE}.
#' @param plot_points_cnt A numeric to indicate the number of points for evaluating and plotting.
#' Default is \code{100}.
#' @param minus_const A numeric to be subtracted in the exponent to
#' ensure the finite-ness of the integration result. Default is \code{0}.
#' @param plot_hist A logical value to indicate whether to the histogram. Default is \code{FALSE}.
#' @param hist_binwidth A numeric to indicate the bin width parameter in plotting the histogram;
#' only works when \code{plot_hist == TRUE}. Default is \code{NULL}, under which condition
#' the Freedman-Diaconis rule is used.
#' @param hist_alpha A numeric to indicate the opacity in plotting the histogram;
#' must range from 0 to 1, inclusively; only works when \code{plot_hist == TRUE}.
#' Default is \code{0.2}.
#'
#' @import ggplot2
#' @import logcondens
#' @import stats
#' @import tidyr
#'
#' @seealso \href{https://cran.r-project.org/web/packages/logcondens/index.html}{logcondens}
#'
#' @return A ggplot2 plot of the log-concave maximum likelihood and (penalized) score matching
#' density estimates over the specified \code{plot_domain}.
#' @export
#'
#' @examples
#' set.seed(1119)
#' N <- 100
#' data <- rnorm(N)
#' domain <- c(-5, 5)
#' result <- lcd_scorematching(data, domain, penalty_param = 1e-1)
#'
#' plot_mle_scorematching(result, result$domain, plot_points_cnt = 500,
#' plot_hist = TRUE)
#'
plot_mle_scorematching <- function(scorematching_logconcave, plot_domain, smoothed = FALSE, plot_log = FALSE,
                                   plot_points_cnt = 100, minus_const = 0,
                                   plot_hist = FALSE, hist_binwidth = NULL, hist_alpha = 0.2) {

    stopifnot(length(plot_domain) == 2)
    if (plot_domain[1] > plot_domain[2]) {
        plot_domain <- sort(plot_domain)
    } else if (plot_domain[1] == plot_domain[2]) {
        stop('The two points in plot_domain are identical. Please provide a meaningful plot_domain.')
    }

    # prepare for the MLE
    data <- scorematching_logconcave$data
    mle <- logcondens::logConDens(data, smoothed = smoothed, print = FALSE)
    newxs <- seq(plot_domain[1], plot_domain[2], length.out = plot_points_cnt)
    mle_result <- logcondens::evaluateLogConDens(newxs, res = mle, which = 1:3)[, c('xs', 'density')]
    colnames(mle_result) <- c('newx_sorted', 'density_vals')

    # prepare for SM
    sm_result <- evaluate_density(
        scorematching_logconcave = scorematching_logconcave,
        newx = newxs,
        minus_const = minus_const)

    merged_result <- merge(x = mle_result, y = sm_result, by = 'newx_sorted')
    colnames(merged_result) <- c('newx_sorted', 'mle_denvals', 'sm_denvals')
    # head(merged_result)

    gather_merged_result <- tidyr::gather(
        merged_result,
        key = 'key',
        value = 'value',
        -newx_sorted)

    if (plot_log == FALSE) {

        density_plot_combined <- ggplot2::ggplot() +
            ggplot2::geom_line(
                data = gather_merged_result,
                ggplot2::aes(x = newx_sorted, y = value, color = key),
                size = 0.8) +
            ggplot2::geom_rug(
                data = data.frame(original_data = scorematching_logconcave$data),
                ggplot2::aes(x = original_data),
                color = 'black',
                alpha = 0.5
            ) +
            ggplot2::scale_color_manual(
                labels = c("Maximum Likelihood", "Score Matching"),
                values = c("deeppink4", #"dodgerblue4",
                           "darkgreen")) +
            ggplot2::labs(
                x = 'x',
                y = 'density',
                title = 'Density Estimate') +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position = "bottom")

        if (plot_hist) {

            # use the Freedman-Diaconis rule for the binwidth
            if (is.null(hist_binwidth)) {
                NN <- length(scorematching_logconcave$data)
                hist_binwidth <- 2 * stats::IQR(scorematching_logconcave$data) / (NN ** (1/3))
            }
            density_plot_combined <- density_plot_combined + ggplot2::geom_histogram(
                data = data.frame(x = data),
                ggplot2::aes(x = x, y = ..density..),
                binwidth = hist_binwidth,
                color = "darkblue",
                fill = "lightblue",
                alpha = hist_alpha)

        }

    } else {

        density_plot_combined <- ggplot2::ggplot() +
            ggplot2::geom_line(
                data = gather_merged_result,
                ggplot2::aes(x = newx_sorted, y = log(value), color = key),
                size = 0.8) +
            ggplot2::geom_rug(
                data = data.frame(original_data = scorematching_logconcave$data),
                ggplot2::aes(x = original_data),
                color = 'black',
                alpha = 0.5
            ) +
            ggplot2::scale_color_manual(
                labels = c("Maximum Likelihood", "Score Matching"),
                values = c("deeppink4", #"dodgerblue4",
                           "darkgreen")) +
            ggplot2::labs(
                x = 'x',
                y = 'density',
                title = 'Density Estimate') +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position = "bottom")

        if (plot_hist) {

            message("You chose to plot the log-density. It does not really make to plot the histogram.")

        }

    }

    return(density_plot_combined)

}

