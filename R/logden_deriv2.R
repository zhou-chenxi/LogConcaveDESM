#' Computation and Plotting the Second Derivative of the Logarithm of the Log-concave Score Matching Density Estimate
#'
#' Based on a "LogConcaveDESM" object, evaluates and plots the second derivative of the logarithm of
#' the (penalized) log-concave score matching density estimate.
#'
#' @param scorematching_logconcave An object of class "LogConcaveDESM",
#' usually the output of \code{\link{lcd_scorematching}} or \code{\link{cv_optimal_density_estimate}}.
#' @param newx A numeric vector of real numbers at which the second derivative of the logarithm of
#' the log-concave score matching density estimate should be evaluated.
#'
#' @details The functions \code{evaluate_logdensity_deriv2_bounded}, \code{evaluate_logdensity_deriv2_R}, \code{evaluate_logdensity_deriv2_ninfb},
#' and \code{evaluate_logdensity_deriv2_ainf} evaluates the second derivative of the logarithm of
#' the (penalized) log-concave score matching density estimate
#' at \code{newx} when the underlying \code{domain} is a bounded interval, the entire real line, an interval of the form \eqn{(-\infty, b)}
#' for some \eqn{b < \infty}, and an interval of the form \eqn{(a, \infty)} for some \eqn{a > -\infty}, respectively.
#' The function \code{evaluate_logdensity_deriv2} encompasses all four cases.
#'
#' The function \code{plot_logdensity_deriv2} plots the second derivative of the logarithm of the (penalized) log-concave score matching density estimate
#' within the plot_domain.
#'
#' @import ggplot2
#' @return A data frame with the first column being the sorted \code{newx} and
#' the second column being the corresponding second derivative values of the logarithm of
#' the (penalized) log-concave score matching density estimate.
#'
#' @examples
#' set.seed(1119)
#' N <- 100
#' data <- rnorm(N)
#' domain <- c(-5, 5)
#' result <- lcd_scorematching(data, domain, penalty_param = 1e-1)
#' # evaluation
#' evaluate_logdensity_deriv2(scorematching_logconcave = result,
#' newx = seq(result$domain[1], result$domain[2], 0.01))
#'
#' # plot
#' plot_logdensity_deriv2(scorematching_logconcave = result,
#' plot_domain = result$domain, plot_points_cnt = 500)
#'
#' @name logdensity_deriv2
NULL

#' @rdname logdensity_deriv2
#' @export
evaluate_logdensity_deriv2 <- function(scorematching_logconcave, newx) {

    # preprocess data
    newx <- as.numeric(newx)
    newx <- newx[!is.nan(newx)]

    domain <- scorematching_logconcave$domain

    domain1 <- domain[1]
    domain2 <- domain[2]

    if ((domain1 == -Inf) & (domain2 == Inf)) {

        # R case
        result <- evaluate_logdensity_deriv2_R(
            scorematching_logconcave = scorematching_logconcave,
            newx = newx)

    } else if (is.finite(domain1) & is.finite(domain2)) {

        # bounded interval case
        result <- evaluate_logdensity_deriv2_bounded(
            scorematching_logconcave = scorematching_logconcave,
            newx = newx)

    } else if (is.finite(domain1) & (domain2 == Inf)) {

        # [a, Inf) case
        result <- evaluate_logdensity_deriv2_ainf(
            scorematching_logconcave = scorematching_logconcave,
            newx = newx)

    }  else if ((domain1 == -Inf) & is.finite(domain2)) {

        # (-Inf, b] case
        result <- evaluate_logdensity_deriv2_ninfb(
            scorematching_logconcave = scorematching_logconcave,
            newx = newx)

    } else {

        stop(paste0("The domain entered, ", domain, ", is not valid."))

    }

    return(result)

}

#' @rdname logdensity_deriv2
#' @export
evaluate_logdensity_deriv2_bounded <- function(scorematching_logconcave, newx) {

    # case of bounded interval

    domain <- scorematching_logconcave$domain
    if (min(newx) < domain[1] || max(newx) > domain[2]) {
        stop('newx is outside of the domain.')
    }

    sorted_data <- scorematching_logconcave$sorted_unique_data
    all_sorted_data <- c(domain[1], sorted_data, domain[2])
    newx <- sort(newx)
    opt_theta <- scorematching_logconcave$opt_theta

    interval_member <- data.frame(
        newx = newx,
        interval = base::findInterval(
            x = newx,
            vec = all_sorted_data,
            rightmost.closed = TRUE,
            left.open = TRUE))

    logderiv2_vals <- vector()

    for (j in unique(interval_member$interval)) {

        theta_i <- opt_theta[j]
        theta_i1 <- opt_theta[j + 1]
        Xi <- all_sorted_data[j]
        Xi1 <- all_sorted_data[j + 1]

        newx_subset <- interval_member[interval_member$interval == j, 'newx']

        vals <- (theta_i1 - theta_i) / (Xi1 - Xi) * (newx_subset - Xi) + theta_i
        logderiv2_vals <- append(logderiv2_vals, vals)

    }

    result <- data.frame(
        newx_sorted = newx,
        logderiv2_vals = -logderiv2_vals
    )

    return(result)
}


#' @rdname logdensity_deriv2
#' @export
evaluate_logdensity_deriv2_R <- function(scorematching_logconcave, newx) {

    domain <- scorematching_logconcave$domain
    # check the validity of the domain
    stopifnot(domain == c(-Inf, Inf))

    sorted_data <- scorematching_logconcave$sorted_unique_data
    all_sorted_data <- sorted_data
    newx <- sort(newx)
    opt_theta <- scorematching_logconcave$opt_theta

    # less than the min and larger than the max
    newx_beyond <- newx[newx < sorted_data[1] | newx > sorted_data[length(sorted_data)]]
    result_beyond <- data.frame(
        newx_sorted = newx_beyond,
        logderiv2_vals = rep(0, length(newx_beyond))
    )

    # between the min and the max
    newx_mid <- newx[newx >= sorted_data[1] & newx <= sorted_data[length(sorted_data)]]
    interval_member <- data.frame(
        newx = newx_mid,
        interval = base::findInterval(
            x = newx_mid,
            vec = sorted_data,
            rightmost.closed = TRUE,
            left.open = TRUE))

    logderiv2_vals <- vector()

    for (j in unique(interval_member$interval)) {

        theta_i <- opt_theta[j]
        theta_i1 <- opt_theta[j + 1]
        Xi <- all_sorted_data[j]
        Xi1 <- all_sorted_data[j + 1]

        newx_subset <- interval_member[interval_member$interval == j, 'newx']

        vals <- (theta_i1 - theta_i) / (Xi1 - Xi) * (newx_subset - Xi) + theta_i
        logderiv2_vals <- append(logderiv2_vals, vals)

    }

    result_mid <- data.frame(
        newx_sorted = newx_mid,
        logderiv2_vals = -logderiv2_vals
    )

    result <- dplyr::bind_rows(result_beyond, result_mid)
    result <- dplyr::arrange(result, newx_sorted)

    return(result)
}


#' @rdname logdensity_deriv2
#' @export
evaluate_logdensity_deriv2_ninfb <- function(scorematching_logconcave, newx) {

    # case of (-infty, b)

    domain <- scorematching_logconcave$domain
    # check the validity of the domain
    stopifnot(domain[1] == -Inf, is.finite(domain[2]))

    if (max(newx) > domain[2]) {
        stop("newx is outside of the domain.")
    }

    sorted_data <- scorematching_logconcave$sorted_unique_data
    all_sorted_data <- c(sorted_data, domain[2])
    newx <- sort(newx)
    opt_theta <- scorematching_logconcave$opt_theta

    # less than the min data
    newx_beyond <- newx[newx < sorted_data[1]]
    result_beyond <- data.frame(
        newx_sorted = newx_beyond,
        logderiv2_vals = rep(0, length(newx_beyond))
    )

    # between the min and b
    newx_mid <- newx[newx >= sorted_data[1]]
    interval_member <- data.frame(
        newx = newx_mid,
        interval = base::findInterval(
            x = newx_mid,
            vec = all_sorted_data,
            rightmost.closed = TRUE,
            left.open = TRUE))

    logderiv2_vals <- vector()

    for (j in unique(interval_member$interval)) {

        theta_i <- opt_theta[j]
        theta_i1 <- opt_theta[j + 1]
        Xi <- all_sorted_data[j]
        Xi1 <- all_sorted_data[j + 1]

        newx_subset <- interval_member[interval_member$interval == j, 'newx']

        vals <- (theta_i1 - theta_i) / (Xi1 - Xi) * (newx_subset - Xi) + theta_i
        logderiv2_vals <- append(logderiv2_vals, vals)

    }

    result_mid <- data.frame(
        newx_sorted = newx_mid,
        logderiv2_vals = -logderiv2_vals
    )

    result <- dplyr::bind_rows(result_beyond, result_mid)
    result <- dplyr::arrange(result, newx_sorted)

    return(result)
}

#' @rdname logdensity_deriv2
#' @export
evaluate_logdensity_deriv2_ainf <- function(scorematching_logconcave, newx) {

    # case of (a, infty)

    domain <- scorematching_logconcave$domain
    # check the validity of the domain
    stopifnot(is.finite(domain[1]), domain[2] == Inf)

    if (min(newx) < domain[1]) {
        stop("newx is outside of the domain.")
    }

    sorted_data <- scorematching_logconcave$sorted_unique_data
    all_sorted_data <- c(domain[1], sorted_data)
    newx <- sort(newx)
    opt_theta <- scorematching_logconcave$opt_theta

    # larger than the max data
    newx_beyond <- newx[newx > sorted_data[length(sorted_data)]]
    result_beyond <- data.frame(
        newx_sorted = newx_beyond,
        logderiv2_vals = rep(0, length(newx_beyond))
    )

    # between the min and b
    newx_mid <- newx[newx <= sorted_data[length(sorted_data)]]
    interval_member <- data.frame(
        newx = newx_mid,
        interval = base::findInterval(
            x = newx_mid,
            vec = all_sorted_data,
            rightmost.closed = TRUE,
            left.open = TRUE))

    logderiv2_vals <- vector()

    for (j in unique(interval_member$interval)) {

        theta_i <- opt_theta[j]
        theta_i1 <- opt_theta[j + 1]
        Xi <- all_sorted_data[j]
        Xi1 <- all_sorted_data[j + 1]

        newx_subset <- interval_member[interval_member$interval == j, 'newx']

        vals <- (theta_i1 - theta_i) / (Xi1 - Xi) * (newx_subset - Xi) + theta_i
        logderiv2_vals <- append(logderiv2_vals, vals)

    }

    result_mid <- data.frame(
        newx_sorted = newx_mid,
        logderiv2_vals = -logderiv2_vals
    )

    result <- dplyr::bind_rows(result_beyond, result_mid)
    result <- dplyr::arrange(result, newx_sorted)

    return(result)
}


#' @rdname logdensity_deriv2
#' @param plot_domain A numeric vector to indicate the domain of the plot.
#' @param plot_points_cnt A numeric to indicate the number of points for evaluating and plotting.
#' Default is \code{100}.
#'
#' @return A ggplot2 plot of the second derivative of the logarithm of
#' the (penalized) log-concave score matching density estimate over the specified plot domain.
#' @export
#'
plot_logdensity_deriv2 <- function(scorematching_logconcave, plot_domain, plot_points_cnt = 100) {

    stopifnot(length(plot_domain) == 2)
    if (plot_domain[1] > plot_domain[2]) {
        plot_domain <- sort(plot_domain)
    } else if (plot_domain[1] == plot_domain[2]) {
        stop('The two points in plot_domain are identical. Please provide a meaningful plot_domain.')
    }

    plot_df <- evaluate_logdensity_deriv2(
        scorematching_logconcave = scorematching_logconcave,
        newx = seq(plot_domain[1], plot_domain[2], length.out = plot_points_cnt)
    )

    plot <- ggplot2::ggplot() +
        ggplot2::geom_line(
            data = plot_df,
            ggplot2::aes(x = newx_sorted, y = logderiv2_vals),
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
            y = 'second derivative',
            title = 'Second Derivative of Log-density Estimate') +
        ggplot2::theme_bw()

    return(plot)

}
