#' Computation and Plotting the First Derivative of the Logarithm of the Log-concave Score Matching Density Estimate
#'
#' Based on a "LogConcaveDESM" object, evaluates and plots the first derivative of the logarithm of
#' the (penalized) log-concave score matching density estimate.
#'
#' @param scorematching_logconcave An object of class "LogConcaveDESM",
#' usually the output of \code{\link{lcd_scorematching}} or \code{\link{cv_optimal_density_estimate}}.
#' @param newx A numeric vector of real numbers at which the first derivative of the logarithm of
#' the (penalized) log-concave score matching density estimate should be evaluated.
#'
#' #' @details The functions \code{evaluate_logdensity_deriv1_bounded}, \code{evaluate_logdensity_deriv1_R},
#' \code{evaluate_logdensity_deriv1_ninfb},
#' and \code{evaluate_logdensity_deriv1_ainf} evaluates the first derivative of the logarithm of
#' the (penalized) log-concave score matching density estimate
#' at \code{newx} when the underlying \code{domain} is a bounded interval, the entire real line, an interval of the form \eqn{(-\infty, b)}
#' for some \eqn{b < \infty},
#' and an interval of the form \eqn{(a, \infty)} for some \eqn{a > -\infty}, respectively.
#' The function \code{evaluate_logdensity_deriv2} encompasses all four cases.
#'
#' The function \code{plot_logdensity_deriv1} plots the first derivative of the logarithm of the (penalized)
#' log-concave score matching density estimate within the plot_domain.
#'
#' @import ggplot2
#' @return A data frame with the first column being the sorted \code{newx} and
#' the second column being the corresponding first derivative values of the logarithm of
#' the (penalized) log-concave score matching density estimate.
#'
#' @examples
#' set.seed(1119)
#' N <- 100
#' data <- rnorm(N)
#' domain <- c(-5, 5)
#' result <- lcd_scorematching(data, domain, penalty_param = 1e-10)
#' # evaluationn
#' evaluate_logdensity_deriv1(scorematching_logconcave = result,
#' newx = seq(result$domain[1], result$domain[2], 0.01))
#'
#' # plot
#' plot_logdensity_deriv1(scorematching_logconcave = result,
#' plot_domain = result$domain, plot_points_cnt = 500)
#'
#' @name logdensity_deriv1
NULL

#' @rdname logdensity_deriv1
#' @export
evaluate_logdensity_deriv1 <- function(scorematching_logconcave, newx) {

    # preprocess data
    newx <- as.numeric(newx)
    newx <- newx[!is.nan(newx)]

    domain <- scorematching_logconcave$domain

    domain1 <- domain[1]
    domain2 <- domain[2]

    if ((domain1 == -Inf) & (domain2 == Inf)) {

        # R case
        result <- evaluate_logdensity_deriv1_R(
            scorematching_logconcave = scorematching_logconcave,
            newx = newx)

    } else if (is.finite(domain1) & is.finite(domain2)) {

        # bounded interval case
        result <- evaluate_logdensity_deriv1_bounded(
            scorematching_logconcave = scorematching_logconcave,
            newx = newx)

    } else if (is.finite(domain1) & (domain2 == Inf)) {

        # [a, Inf) case
        result <- evaluate_logdensity_deriv1_ainf(
            scorematching_logconcave = scorematching_logconcave,
            newx = newx)

    }  else if ((domain1 == -Inf) & is.finite(domain2)) {

        # (-Inf, b] case
        result <- evaluate_logdensity_deriv1_ninfb(
            scorematching_logconcave = scorematching_logconcave,
            newx = newx)

    } else {

        stop(paste0("The domain entered, ", domain, ", is not valid."))

    }

    return(result)

}

#' @rdname logdensity_deriv1
#' @export
evaluate_logdensity_deriv1_bounded <- function(scorematching_logconcave, newx) {

    domain <- scorematching_logconcave$domain
    stopifnot(all(is.finite(domain)))
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
            val1 <- (sum(opt_theta[1:(j - 1)] * data_subset_j_diff) +
                         sum(opt_theta[2:j] * data_subset_j_diff)) / 2

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

#' @rdname logdensity_deriv1
#' @export
evaluate_logdensity_deriv1_R <- function(scorematching_logconcave, newx) {

    domain <- scorematching_logconcave$domain
    stopifnot(domain == c(-Inf, Inf))

    sorted_data <- scorematching_logconcave$sorted_unique_data
    all_sorted_data <- sorted_data
    newx <- sort(newx)
    opt_theta <- scorematching_logconcave$opt_theta

    # compute (g'(-Inf))^*
    weighted_col_AB <- sweep(
        scorematching_logconcave$matrix_A + scorematching_logconcave$matrix_B,
        MARGIN = 2,
        STATS = scorematching_logconcave$data_weights[2:length(scorematching_logconcave$data_weights)],
        FUN = '*'
    )
    g_a_opt <- -sum(opt_theta * rowSums(weighted_col_AB)) / 2

    # less than X_{(1)} part
    newx_part1 <- newx[newx < all_sorted_data[1]]
    if (length(newx_part1) == 0) {

        result1 <- data.frame()

    } else {

        result1 <- data.frame(
            newx_sorted = newx_part1,
            logderiv1_vals = g_a_opt
        )

    }

    # between X_{(1)} and X_{(m)}
    newx_part2 <- newx[newx >= all_sorted_data[1] & newx <= all_sorted_data[length(all_sorted_data)]]
    if (length(newx_part2) == 0) {

        result2 <- data.frame()

    } else {

        interval_member <- data.frame(
            newx = newx_part2,
            interval = base::findInterval(
                x = newx_part2,
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
                val1 <- (sum(opt_theta[1:(j - 1)] * data_subset_j_diff) +
                             sum(opt_theta[2:j] * data_subset_j_diff)) / 2

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

        result2 <- data.frame(
            newx_sorted = newx_part2,
            logderiv1_vals = g_a_opt + logderiv1_vals
        )

    }


    # larger than X_{(m)}
    newx_part3 <- newx[newx > all_sorted_data[length(all_sorted_data)]]
    if (length(newx_part3) == 0) {

        result3 <- data.frame()

    } else {

        data_subset_all_diff <- diff(all_sorted_data)
        term2 <- (sum(opt_theta[1:(length(opt_theta) - 1)] * data_subset_all_diff) +
                      sum(opt_theta[2:length(opt_theta)] * data_subset_all_diff)) / 2
        result3 <- data.frame(
            newx_sorted = newx_part3,
            logderiv1_vals = g_a_opt + term2
        )

    }

    result <- dplyr::bind_rows(result1, result2, result3)
    result <- dplyr::arrange(result, newx_sorted)
    result$logderiv1_vals <- result$logderiv1_vals * (-1)

    return(result)
}

#' @rdname logdensity_deriv1
#' @export
evaluate_logdensity_deriv1_ninfb <- function(scorematching_logconcave, newx) {

    domain <- scorematching_logconcave$domain
    stopifnot(domain[1] == -Inf, is.finite(domain[2]))

    if (max(newx) > domain[2]) {
        stop("newx is outside of the domain.")
    }

    sorted_data <- scorematching_logconcave$sorted_unique_data
    all_sorted_data <- c(sorted_data, domain[2])
    newx <- sort(newx)
    opt_theta <- scorematching_logconcave$opt_theta

    # compute (g'(-Inf))^*
    weighted_col_AB <- sweep(
        scorematching_logconcave$matrix_A + scorematching_logconcave$matrix_B,
        MARGIN = 2,
        STATS = scorematching_logconcave$data_weights[2:length(scorematching_logconcave$data_weights)],
        FUN = '*'
    )
    g_a_opt <- -sum(opt_theta * rowSums(weighted_col_AB)) / 2

    # less than X_{(1)} part
    newx_part1 <- newx[newx < all_sorted_data[1]]
    if (length(newx_part1) == 0) {

        result1 <- data.frame()

    } else {

        result1 <- data.frame(
            newx_sorted = newx_part1,
            logderiv1_vals = g_a_opt
        )

    }

    # between X_{(1)} and X_{(m+1)}
    newx_part2 <- newx[newx >= all_sorted_data[1]]
    if (length(newx_part2) == 0) {

        result2 <- data.frame()

    } else {

        interval_member <- data.frame(
            newx = newx_part2,
            interval = base::findInterval(
                x = newx_part2,
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
                val1 <- (sum(opt_theta[1:(j - 1)] * data_subset_j_diff) +
                             sum(opt_theta[2:j] * data_subset_j_diff)) / 2

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

        result2 <- data.frame(
            newx_sorted = newx_part2,
            logderiv1_vals = g_a_opt + logderiv1_vals
        )

    }

    result <- dplyr::bind_rows(result1, result2)
    result <- dplyr::arrange(result, newx_sorted)
    result$logderiv1_vals <- result$logderiv1_vals * (-1)

    return(result)
}

#' @rdname logdensity_deriv1
#' @export
evaluate_logdensity_deriv1_ainf <- function(scorematching_logconcave, newx) {

    domain <- scorematching_logconcave$domain
    stopifnot(is.finite(domain[1]), domain[2] == Inf)

    if (min(newx) < domain[1]) {
        stop("newx is outside of the domain.")
    }

    sorted_data <- scorematching_logconcave$sorted_unique_data
    all_sorted_data <- c(domain[1], sorted_data)
    newx <- sort(newx)
    opt_theta <- scorematching_logconcave$opt_theta

    # compute (g'(-Inf))^*
    weighted_col_AB <- sweep(
        scorematching_logconcave$matrix_A + scorematching_logconcave$matrix_B,
        MARGIN = 2,
        STATS = scorematching_logconcave$data_weights,
        FUN = '*'
    )
    g_a_opt <- -sum(opt_theta * rowSums(weighted_col_AB)) / 2

    # between a = X_{(0)} and X_{(m)}
    newx_part2 <- newx[newx <= all_sorted_data[length(all_sorted_data)]]
    if (length(newx_part2) == 0) {

        result2 <- data.frame()

    } else {

        interval_member <- data.frame(
            newx = newx_part2,
            interval = base::findInterval(
                x = newx_part2,
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
                val1 <- (sum(opt_theta[1:(j - 1)] * data_subset_j_diff) +
                             sum(opt_theta[2:j] * data_subset_j_diff)) / 2

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

        result2 <- data.frame(
            newx_sorted = newx_part2,
            logderiv1_vals = g_a_opt + logderiv1_vals
        )

    }

    # larger than X_{(m)}
    newx_part3 <- newx[newx > all_sorted_data[length(all_sorted_data)]]
    if (length(newx_part3) == 0) {

        result3 <- data.frame()

    } else {

        data_subset_all_diff <- diff(all_sorted_data)
        term2 <- (sum(opt_theta[1:(length(opt_theta) - 1)] * data_subset_all_diff) +
                      sum(opt_theta[2:length(opt_theta)] * data_subset_all_diff)) / 2
        result3 <- data.frame(
            newx_sorted = newx_part3,
            logderiv1_vals = g_a_opt + term2
        )

    }

    result <- dplyr::bind_rows(result2, result3)
    result <- dplyr::arrange(result, newx_sorted)
    result$logderiv1_vals <- result$logderiv1_vals * (-1)

    return(result)
}

#' @rdname logdensity_deriv1
#' @param plot_domain A numeric vector to indicate the domain of the plot.
#' @param plot_points_cnt A numeric to indicate the number of points for evaluating and plotting.
#' Default is \code{100}.
#'
#' @return A ggplot2 plot of the first derivative of the logarithm of
#' the (penalized) log-concave score matching density estimate over the specified plot domain.
#' @export
#'
plot_logdensity_deriv1 <- function(scorematching_logconcave, plot_domain, plot_points_cnt = 100) {

    stopifnot(length(plot_domain) == 2)
    if (plot_domain[1] > plot_domain[2]) {
        plot_domain <- sort(plot_domain)
    } else if (plot_domain[1] == plot_domain[2]) {
        stop('The two points in plot_domain are identical. Please provide a meaningful plot_domain.')
    }

    plot_df <- evaluate_logdensity_deriv1(
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

    return(plot)

}
