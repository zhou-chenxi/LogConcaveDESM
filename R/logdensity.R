#' Computation and Plotting the Logarithm of the Un-normalized Log-concave Score Matching Density Estimate
#'
#' Based on a "LogConcaveDESM" object, evaluates and plots the logarithm of
#' the (penalized) log-concave score matching density estimate up to a normalizing constant.
#'
#' @param scorematching_logconcave An object of class "LogConcaveDESM",
#' usually the output of \code{\link{lcd_scorematching}} or \code{\link{cv_optimal_density_estimate}}.
#' @param newx A numeric vector of real numbers at which the logarithm of
#' the (penalized) log-concave score matching density estimate should be evaluated.
#'
#' @details The functions \code{evaluate_logdensity_bounded}, \code{evaluate_logdensity_R}, \code{evaluate_logdensity_ninfb},
#' and \code{evaluate_logdensity_ainf} evaluates the logarithm of the (penalized) log-concave score matching density estimate
#' at \code{newx} when the underlying \code{domain} is a bounded interval, the entire real line, an interval of the form \eqn{(-\infty, b)}
#' for some \eqn{b < \infty}, and an interval of the form \eqn{(a, \infty)} for some \eqn{a > -\infty}, respectively.
#' The function \code{evaluate_logdensity} encompasses all four cases.
#'
#' The function \code{plot_logdensity} plots the logarithm of the (penalized) log-concave score matching density estimate
#' within the plot_domain.
#'
#' @import ggplot2
#' @return A data frame with the first column being the sorted \code{newx} and
#' the second column being the corresponding values of the logarithm of
#' the (penalized) log-concave score matching density estimate.
#'
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
#' @name logdensity
NULL

#' @rdname logdensity
#' @export
evaluate_logdensity <- function(scorematching_logconcave, newx) {

    # preprocess data
    newx <- as.numeric(newx)
    newx <- newx[!is.nan(newx)]

    domain <- scorematching_logconcave$domain

    domain1 <- domain[1]
    domain2 <- domain[2]

    if ((domain1 == -Inf) & (domain2 == Inf)) {

        # R case
        result <- evaluate_logdensity_R(
            scorematching_logconcave = scorematching_logconcave,
            newx = newx)

    } else if (is.finite(domain1) & is.finite(domain2)) {

        # bounded interval case
        result <- evaluate_logdensity_bounded(
            scorematching_logconcave = scorematching_logconcave,
            newx = newx)

    } else if (is.finite(domain1) & (domain2 == Inf)) {

        # [a, Inf) case
        result <- evaluate_logdensity_ainf(
            scorematching_logconcave = scorematching_logconcave,
            newx = newx)

    }  else if ((domain1 == -Inf) & is.finite(domain2)) {

        # (-Inf, b] case
        result <- evaluate_logdensity_ninfb(
            scorematching_logconcave = scorematching_logconcave,
            newx = newx)

    } else {

        stop(paste0("The domain entered, ", domain, ", is not valid."))

    }

    return(result)

}

#' @rdname logdensity
#' @export
evaluate_logdensity_bounded <- function(scorematching_logconcave, newx) {

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
            term2_coef <- (sum(opt_theta[1:(j - 1)] * data_subset_j_diff) + sum(opt_theta[2:j] * data_subset_j_diff)) / 2
            term2_1 <- term2_coef * newx_subset
            term2_2 <- (sum(opt_theta[1:(j - 1)] * data_subset_j_diff * all_sorted_data[2:j]) +
                            sum(opt_theta[2:j] * data_subset_j_diff * all_sorted_data[2:j])) / 2
            term2 <- term2_1 - term2_2

            # term 3
            term3 <- sum(diff(opt_theta[1:j]) * (data_subset_j_diff) ** 2) / 6

            # term 4
            term4 <- sum(opt_theta[1:(j - 1)] * (data_subset_j_diff) ** 2) / 2

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

#' @rdname logdensity
#' @export
evaluate_logdensity_R <- function(scorematching_logconcave, newx) {

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
            logdensity_vals = g_a_opt * (newx_part1 - all_sorted_data[1])
        )

    }

    # between X_{(1)} and X_{(m)}
    newx_part2 <- newx[newx >= all_sorted_data[1] & newx <= all_sorted_data[length(all_sorted_data)]]
    if (length(newx_part2) == 0) {

        result2 <- data.frame()

    } else {

        # term 2
        interval_member <- data.frame(
            newx = newx_part2,
            interval = base::findInterval(
                x = newx_part2,
                vec = all_sorted_data,
                rightmost.closed = TRUE,
                left.open = TRUE))

        logdensity_vals <- vector()

        for (j in unique(interval_member$interval)) {

            newx_subset <- interval_member[interval_member$interval == j, 'newx']

            term1 <- g_a_opt * (newx_subset - all_sorted_data[1])

            # summation terms involving the previous intervals
            if (j == 1) {

                val1 <- 0

            } else {

                data_subset_j_diff <- diff(all_sorted_data[1:j])
                val1 <- (sum(opt_theta[1:(j - 1)] * data_subset_j_diff) +
                             sum(opt_theta[2:j] * data_subset_j_diff)) / 2

            }

            term2 <- val1 * newx_subset

            # terms involving newx - Xi
            theta_i <- opt_theta[j]
            theta_i1 <- opt_theta[j + 1]
            Xi <- all_sorted_data[j]
            Xi1 <- all_sorted_data[j + 1]

            term3 <- ((theta_i1 - theta_i) / (Xi1 - Xi) * (newx_subset - Xi) ** 3 / 6 +
                         theta_i * (newx_subset - Xi) ** 2 / 2)

            # additional added terms
            if (j == 1) {

                term4 <- 0

            } else {

                data_subset_j_diff <- diff(all_sorted_data[1:j])
                theta_subset_2j <- opt_theta[2:j]
                theta_subset_1j <- opt_theta[1:(j - 1)]
                term4 <- sum((theta_subset_2j / 6 - theta_subset_1j / 6 + theta_subset_1j / 2) * data_subset_j_diff ** 2 -
                                (theta_subset_1j + theta_subset_2j) * data_subset_j_diff * all_sorted_data[2:j] / 2)

            }

            logdensity_vals <- append(
                logdensity_vals,
                term1 + term2 + term3 + term4)

        }

        result2 <- data.frame(
            newx_sorted = newx_part2,
            logdensity_vals = logdensity_vals
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

        theta_subset_2 <- opt_theta[2:length(opt_theta)]
        theta_subset_1 <- opt_theta[1:(length(opt_theta) - 1)]
        term3 <- sum((theta_subset_2 / 6 - theta_subset_1 / 6 + theta_subset_1 / 2) * data_subset_all_diff ** 2 -
                         (theta_subset_1 + theta_subset_2) * data_subset_all_diff * all_sorted_data[2:length(all_sorted_data)] / 2)

        result3 <- data.frame(
            newx_sorted = newx_part3,
            logdensity_vals = g_a_opt * (newx_part3 - all_sorted_data[1]) + term2 * newx_part3 + term3
        )

    }

    result <- dplyr::bind_rows(result1, result2, result3)
    result <- dplyr::arrange(result, newx_sorted)
    result$logdensity_vals <- result$logdensity_vals * (-1)

    return(result)
}

#' @rdname logdensity
#' @export
evaluate_logdensity_ninfb <- function(scorematching_logconcave, newx) {

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
            logdensity_vals = g_a_opt * (newx_part1 - all_sorted_data[1])
        )

    }

    # between X_{(1)} and X_{(m+1)}
    newx_part2 <- newx[newx >= all_sorted_data[1]]
    if (length(newx_part2) == 0) {

        result2 <- data.frame()

    } else {

        ## term 1
        # term1 <- g_a_opt * (newx_part2 - all_sorted_data[1])

        # term 2
        interval_member <- data.frame(
            newx = newx_part2,
            interval = base::findInterval(
                x = newx_part2,
                vec = all_sorted_data,
                rightmost.closed = TRUE,
                left.open = TRUE))

        logdensity_vals <- vector()

        for (j in unique(interval_member$interval)) {

            newx_subset <- interval_member[interval_member$interval == j, 'newx']
            term1 <- g_a_opt * (newx_subset - all_sorted_data[1])

            # summation terms involving the previous intervals
            if (j == 1) {

                val1 <- 0

            } else {

                data_subset_j_diff <- diff(all_sorted_data[1:j])
                val1 <- (sum(opt_theta[1:(j - 1)] * data_subset_j_diff) +
                             sum(opt_theta[2:j] * data_subset_j_diff)) / 2

            }
            term2 <- val1 * newx_subset

            # terms involving newx - Xi
            theta_i <- opt_theta[j]
            theta_i1 <- opt_theta[j + 1]
            Xi <- all_sorted_data[j]
            Xi1 <- all_sorted_data[j + 1]

            term3 <- ((theta_i1 - theta_i) / (Xi1 - Xi) * (newx_subset - Xi) ** 3 / 6 +
                         theta_i * (newx_subset - Xi) ** 2 / 2)

            # additional added terms
            if (j == 1) {

                term4 <- 0

            } else {

                data_subset_j_diff <- diff(all_sorted_data[1:j])
                theta_subset_2j <- opt_theta[2:j]
                theta_subset_1j <- opt_theta[1:(j - 1)]
                term4 <- sum((theta_subset_2j / 6 - theta_subset_1j / 6 + theta_subset_1j / 2) * data_subset_j_diff ** 2 -
                                (theta_subset_1j + theta_subset_2j) * data_subset_j_diff * all_sorted_data[2:j] / 2)

            }

            logdensity_vals <- append(logdensity_vals, term1 + term2 + term3 + term4)

        }

        result2 <- data.frame(
            newx_sorted = newx_part2,
            logdensity_vals = logdensity_vals
        )

    }

    result <- dplyr::bind_rows(result1, result2)
    result <- dplyr::arrange(result, newx_sorted)
    result$logdensity_vals <- result$logdensity_vals * (-1)

    return(result)
}

#' @rdname logdensity
#' @export
evaluate_logdensity_ainf <- function(scorematching_logconcave, newx) {

    domain <- scorematching_logconcave$domain
    stopifnot(is.finite(domain[1]), domain[2] == Inf)

    if (min(newx) < domain[1]) {
        stop("newx is outside of the domain.")
    }

    sorted_data <- scorematching_logconcave$sorted_unique_data
    all_sorted_data <- c(domain[1], sorted_data)
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

    # between X_{(1)} and X_{(m)}
    newx_part2 <- newx[newx <= all_sorted_data[length(all_sorted_data)]]
    if (length(newx_part2) == 0) {

        result2 <- data.frame()

    } else {

        ## term 1
        # term1 <- g_a_opt * (newx_part2 - all_sorted_data[1])

        # term 2
        interval_member <- data.frame(
            newx = newx_part2,
            interval = base::findInterval(
                x = newx_part2,
                vec = all_sorted_data,
                rightmost.closed = TRUE,
                left.open = TRUE))

        logdensity_vals <- vector()

        for (j in unique(interval_member$interval)) {

            newx_subset <- interval_member[interval_member$interval == j, 'newx']

            term1 <- g_a_opt * (newx_subset - all_sorted_data[1])

            # summation terms involving the previous intervals
            if (j == 1) {

                val1 <- 0

            } else {

                data_subset_j_diff <- diff(all_sorted_data[1:j])
                val1 <- (sum(opt_theta[1:(j - 1)] * data_subset_j_diff) +
                             sum(opt_theta[2:j] * data_subset_j_diff)) / 2

            }

            term2 <- val1 * newx_subset

            # terms involving newx - Xi
            theta_i <- opt_theta[j]
            theta_i1 <- opt_theta[j + 1]
            Xi <- all_sorted_data[j]
            Xi1 <- all_sorted_data[j + 1]

            term3 <- ((theta_i1 - theta_i) / (Xi1 - Xi) * (newx_subset - Xi) ** 3 / 6 +
                         theta_i * (newx_subset - Xi) ** 2 / 2)

            # additional added terms
            if (j == 1) {

                term4 <- 0

            } else {

                data_subset_j_diff <- diff(all_sorted_data[1:j])
                theta_subset_2j <- opt_theta[2:j]
                theta_subset_1j <- opt_theta[1:(j - 1)]
                term4 <- sum((theta_subset_2j / 6 - theta_subset_1j / 6 + theta_subset_1j / 2) * data_subset_j_diff ** 2 -
                                (theta_subset_1j + theta_subset_2j) * data_subset_j_diff * all_sorted_data[2:j] / 2)

            }

            logdensity_vals <- append(logdensity_vals, term1 + term2 + term3 + term4)

        }

        result2 <- data.frame(
            newx_sorted = newx_part2,
            logdensity_vals = logdensity_vals
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

        theta_subset_2 <- opt_theta[2:length(opt_theta)]
        theta_subset_1 <- opt_theta[1:(length(opt_theta) - 1)]
        term3 <- sum((theta_subset_2 / 6 - theta_subset_1 / 6 + theta_subset_1 / 2) * data_subset_all_diff ** 2 -
                         (theta_subset_1 + theta_subset_2) * data_subset_all_diff * all_sorted_data[2:length(all_sorted_data)] / 2)

        result3 <- data.frame(
            newx_sorted = newx_part3,
            logdensity_vals = g_a_opt * (newx_part3 - all_sorted_data[1]) + term2 * newx_part3 + term3
        )

    }

    result <- dplyr::bind_rows(result2, result3)
    result <- dplyr::arrange(result, newx_sorted)
    result$logdensity_vals <- result$logdensity_vals * (-1)

    return(result)

}

#' @rdname logdensity
#' @param plot_domain A numeric vector to indicate the domain of the plot.
#' @param plot_points_cnt A numeric to indicate the number of points for evaluating and plotting.
#' Default is \code{100}.
#'
#' @return A ggplot2 plot of the logarithm of
#' the (penalized) log-concave score matching density estimate over the specified plot domain.
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

    return(plot)

}
