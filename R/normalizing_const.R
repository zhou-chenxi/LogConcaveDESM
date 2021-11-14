#' Computation the Normalizing Constant
#'
#' Computes the normalizing constant of a (penalized) un-normalized log-concave score matching density estimate.
#'
#' @param scorematching_logconcave An object of class "LogConcaveDESM",
#' usually the output of \code{\link{lcd_scorematching}} or \code{\link{cv_optimal_density_estimate}}.
#' @param minus_const A numeric to be subtracted in the exponent to
#' ensure the finite-ness of the integration result. Default is \code{0}.
#'
#' @details The functions \code{normalizing_const_bounded}, \code{normalizing_const_R}, \code{normalizing_const_ninfb},
#' and \code{normalizing_const_ainf} computes the normalizing constant of a (penalized) un-normalized
#' log-concave score matching density estimate when the underlying \code{domain} is a bounded interval,
#' the entire real line, an interval of the form \eqn{(-\infty, b)}
#' for some \eqn{b < \infty}, and an interval of the form \eqn{(a, \infty)} for some \eqn{a > -\infty}, respectively.
#' The function \code{normalizing_const} encompasses all four cases.
#'
#' @return The normalizing constant of a (penalized) un-normalized log-concave score matching density estimate.
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
#' @name normalizing_const
NULL

#' @rdname normalizing_const
#' @export
normalizing_const <- function(scorematching_logconcave, minus_const = 0) {

    domain <- scorematching_logconcave$domain

    domain1 <- domain[1]
    domain2 <- domain[2]

    if ((domain1 == -Inf) & (domain2 == Inf)) {

        # R case
        result <- normalizing_const_R(
            scorematching_logconcave = scorematching_logconcave,
            minus_const = minus_const)

    } else if (is.finite(domain1) & is.finite(domain2)) {

        # bounded interval case
        result <- normalizing_const_bounded(
            scorematching_logconcave = scorematching_logconcave,
            minus_const = minus_const)

    } else if (is.finite(domain1) & (domain2 == Inf)) {

        # [a, Inf) case
        result <- normalizing_const_ainf(
            scorematching_logconcave = scorematching_logconcave,
            minus_const = minus_const)

    }  else if ((domain1 == -Inf) & is.finite(domain2)) {

        # (-Inf, b] case
        result <- normalizing_const_ninfb(
            scorematching_logconcave = scorematching_logconcave,
            minus_const = minus_const)

    } else {

        stop(paste0("The domain entered, ", domain, ", is not valid."))

    }

    return(result)

}

#' @rdname normalizing_const
#' @export
normalizing_const_bounded <- function(scorematching_logconcave, minus_const = 0) {

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
    for (i in 1:(length(sorted_data) + 1)) {

        a0 <- (opt_theta[i + 1] - opt_theta[i]) / (all_sorted_data[i + 1] - all_sorted_data[i])

        a <- a0 / 6

        b <- (opt_theta[i] - a0 * all_sorted_data[i]) / 2

        if (i == 1) {

            c <- a0 * all_sorted_data[i] ** 2 / 2 - opt_theta[i] * all_sorted_data[i] + g_a_opt

            d <- (opt_theta[i] * all_sorted_data[i] ** 2 / 2 - domain[1] * g_a_opt -
                      a0 * all_sorted_data[i] ** 3 / 6)

        } else {

            data_subset_i_diff <- diff(all_sorted_data[1:i])
            term2_1 <- (sum(opt_theta[1:(i - 1)] * data_subset_i_diff) +
                            sum(opt_theta[2:i] * data_subset_i_diff)) / 2

            c <- (a0 * all_sorted_data[i] ** 2 / 2 - opt_theta[i] * all_sorted_data[i] + g_a_opt +
                      term2_1)

            term2_2 <- (sum(opt_theta[1:(i - 1)] * data_subset_i_diff * all_sorted_data[2:i]) +
                            sum(opt_theta[2:i] * data_subset_i_diff * all_sorted_data[2:i])) / 2
            term3 <- sum(diff(opt_theta[1:i]) * (data_subset_i_diff) ** 2) / 6
            term4 <- sum(opt_theta[1:(i - 1)] * (data_subset_i_diff) ** 2) / 2

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

#' @rdname normalizing_const
#' @export
normalizing_const_R <- function(scorematching_logconcave, minus_const = 0) {

    sorted_data <- scorematching_logconcave$sorted_unique_data
    domain <- scorematching_logconcave$domain
    stopifnot(domain == c(-Inf, Inf))
    all_sorted_data <- sorted_data
    opt_theta <- scorematching_logconcave$opt_theta

    # compute (g'(a))^*
    weighted_col_AB <- sweep(
        scorematching_logconcave$matrix_A + scorematching_logconcave$matrix_B,
        MARGIN = 2,
        STATS = scorematching_logconcave$data_weights[2:length(scorematching_logconcave$data_weights)],
        FUN = '*'
    )
    g_a_opt <- -sum(opt_theta * rowSums(weighted_col_AB)) / 2

    # compute g_b_opt = \hat{c} + \frac{1}{2} \sum_{j=1}^{i-1} \parens{\theta_{j} + \theta_{j+1}} \parens{X_{j+1} - X_{j}}
    data_subset_all_diff <- diff(all_sorted_data)
    term2 <- (sum(opt_theta[1:(length(opt_theta) - 1)] * data_subset_all_diff) +
                  sum(opt_theta[2:length(opt_theta)] * data_subset_all_diff)) / 2
    g_b_opt <- g_a_opt + term2

    # additive constant on (X_{(m)}, \infty)
    theta_subset_2 <- opt_theta[2:length(opt_theta)]
    theta_subset_1 <- opt_theta[1:(length(opt_theta) - 1)]
    tm <- sum((theta_subset_2 / 6 - theta_subset_1 / 6 + theta_subset_1 / 2) * data_subset_all_diff ** 2 -
                  (theta_subset_1 + theta_subset_2) * data_subset_all_diff * all_sorted_data[2:length(all_sorted_data)] / 2)

    result <- 0
    # integrate the region < X_{(1)}
    result <- result + integrate_expcubic(
        a = 0, b = 0, c = -g_a_opt, d = g_a_opt * all_sorted_data[1],
        lower = -Inf,
        upper = all_sorted_data[1],
        minus_const = minus_const)

    # integrate the region > X_{(m)}
    result <- result + integrate_expcubic(
        a = 0, b = 0, c = -g_b_opt, d = -(tm - g_a_opt * all_sorted_data[1]),
        lower = all_sorted_data[length(all_sorted_data)],
        upper = Inf,
        minus_const = minus_const)

    # integrate each sub-interval
    for (i in 1:(length(all_sorted_data) - 1)) {

        a0 <- (opt_theta[i + 1] - opt_theta[i]) / (all_sorted_data[i + 1] - all_sorted_data[i])

        a <- a0 / 6

        b <- (opt_theta[i] - a0 * all_sorted_data[i]) / 2

        if (i == 1) {

            c <- a0 * all_sorted_data[i] ** 2 / 2 - opt_theta[i] * all_sorted_data[i] + g_a_opt

            d <- (opt_theta[i] * all_sorted_data[i] ** 2 / 2 - g_a_opt * all_sorted_data[1] -
                      a0 * all_sorted_data[i] ** 3 / 6)

        } else {

            data_subset_i_diff <- diff(all_sorted_data[1:i])
            term2_1 <- (sum(opt_theta[1:(i - 1)] * data_subset_i_diff) +
                            sum(opt_theta[2:i] * data_subset_i_diff)) / 2

            c <- (a0 * all_sorted_data[i] ** 2 / 2 - opt_theta[i] * all_sorted_data[i] + g_a_opt +
                      term2_1)

            opt_theta_sub1 <- opt_theta[1:(i - 1)]
            opt_theta_sub2 <- opt_theta[2:i]

            term2_2 <- (sum(opt_theta_sub1 * data_subset_i_diff * all_sorted_data[2:i]) +
                            sum(opt_theta_sub2 * data_subset_i_diff * all_sorted_data[2:i])) / 2

            term3 <- sum(((opt_theta_sub2 - opt_theta_sub1) / 6 + opt_theta_sub1 / 2) * data_subset_i_diff ** 2)

            d <- (opt_theta[i] * all_sorted_data[i] ** 2 / 2 - g_a_opt * all_sorted_data[1] -
                      a0 * all_sorted_data[i] ** 3 / 6 - term2_2 + term3)

        }

        result <- result + integrate_expcubic(
            -a, -b, -c, -d,
            lower = all_sorted_data[i],
            upper = all_sorted_data[i + 1],
            minus_const = minus_const
        )

    }

    return(result)

}

#' @rdname normalizing_const
#' @export
normalizing_const_ninfb <- function(scorematching_logconcave, minus_const = 0) {

    sorted_data <- scorematching_logconcave$sorted_unique_data
    domain <- scorematching_logconcave$domain
    stopifnot(domain[1] == -Inf, is.finite(domain[2]))
    all_sorted_data <- c(sorted_data, domain[2])
    opt_theta <- scorematching_logconcave$opt_theta

    # compute (g'(a))^*
    weighted_col_AB <- sweep(
        scorematching_logconcave$matrix_A + scorematching_logconcave$matrix_B,
        MARGIN = 2,
        STATS = scorematching_logconcave$data_weights[2:length(scorematching_logconcave$data_weights)],
        FUN = '*'
    )
    g_a_opt <- -sum(opt_theta * rowSums(weighted_col_AB)) / 2

    result <- 0
    # integrate from -Inf to X_{(1)}
    result <- result + integrate_expcubic(
        a = 0, b = 0, c = -g_a_opt, d = g_a_opt * all_sorted_data[1],
        lower = -Inf,
        upper = all_sorted_data[1],
        minus_const = minus_const
    )

    # integrate each sub-interval
    for (i in 1:length(sorted_data)) {

        a0 <- (opt_theta[i + 1] - opt_theta[i]) / (all_sorted_data[i + 1] - all_sorted_data[i])

        a <- a0 / 6

        b <- (opt_theta[i] - a0 * all_sorted_data[i]) / 2

        if (i == 1) {

            c <- a0 * all_sorted_data[i] ** 2 / 2 - opt_theta[i] * all_sorted_data[i] + g_a_opt

            d <- (opt_theta[i] * all_sorted_data[i] ** 2 / 2 - all_sorted_data[1] * g_a_opt -
                      a0 * all_sorted_data[i] ** 3 / 6)

        } else {

            data_subset_i_diff <- diff(all_sorted_data[1:i])
            term2_1 <- (sum(opt_theta[1:(i - 1)] * data_subset_i_diff) +
                            sum(opt_theta[2:i] * data_subset_i_diff)) / 2

            c <- (a0 * all_sorted_data[i] ** 2 / 2 - opt_theta[i] * all_sorted_data[i] + g_a_opt +
                      term2_1)

            term2_2 <- sum((opt_theta[1:(i - 1)] + opt_theta[2:i]) * data_subset_i_diff * all_sorted_data[2:i]) / 2
            term3 <- sum((opt_theta[2:i] - opt_theta[1:(i - 1)]) * (data_subset_i_diff ** 2)) / 6
            term4 <- sum(opt_theta[1:(i - 1)] * (data_subset_i_diff ** 2)) / 2

            d <- (opt_theta[i] * all_sorted_data[i] ** 2 / 2 - all_sorted_data[1] * g_a_opt -
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

#' @rdname normalizing_const
#' @export
normalizing_const_ainf <- function(scorematching_logconcave, minus_const = 0) {

    sorted_data <- scorematching_logconcave$sorted_unique_data
    domain <- scorematching_logconcave$domain
    stopifnot(is.finite(domain[1]), domain[2] == Inf)
    all_sorted_data <- c(domain[1], sorted_data)
    opt_theta <- scorematching_logconcave$opt_theta

    # compute (g'(a))^*
    weighted_col_AB <- sweep(
        scorematching_logconcave$matrix_A + scorematching_logconcave$matrix_B,
        MARGIN = 2,
        STATS = scorematching_logconcave$data_weights,
        FUN = '*'
    )
    g_a_opt <- -sum(opt_theta * rowSums(weighted_col_AB)) / 2

    result <- 0
    # integrate each sub-interval
    for (i in 1:length(sorted_data)) {

        a0 <- (opt_theta[i + 1] - opt_theta[i]) / (all_sorted_data[i + 1] - all_sorted_data[i])

        a <- a0 / 6

        b <- (opt_theta[i] - a0 * all_sorted_data[i]) / 2

        if (i == 1) {

            c <- a0 * all_sorted_data[i] ** 2 / 2 - opt_theta[i] * all_sorted_data[i] + g_a_opt

            d <- (opt_theta[i] * all_sorted_data[i] ** 2 / 2 - domain[1] * g_a_opt -
                      a0 * all_sorted_data[i] ** 3 / 6)

        } else {

            data_subset_i_diff <- diff(all_sorted_data[1:i])
            term2_1 <- sum(opt_theta[1:(i - 1)] * data_subset_i_diff + opt_theta[2:i] * data_subset_i_diff) / 2
            c <- (a0 * all_sorted_data[i] ** 2 / 2 - opt_theta[i] * all_sorted_data[i] + g_a_opt + term2_1)

            term2_2 <- sum(opt_theta[1:(i - 1)] * data_subset_i_diff * all_sorted_data[2:i] +
                               opt_theta[2:i] * data_subset_i_diff * all_sorted_data[2:i]) / 2

            opt_theta_sub1 <- opt_theta[1:(i - 1)]
            opt_theta_sub2 <- opt_theta[2:i]
            term3 <- sum((opt_theta_sub2 / 6 - opt_theta_sub1 / 6 + opt_theta_sub1 / 2) * (data_subset_i_diff ** 2))

            d <- (opt_theta[i] * all_sorted_data[i] ** 2 / 2 - domain[1] * g_a_opt -
                      a0 * all_sorted_data[i] ** 3 / 6 - term2_2 + term3)

        }

        result <- result + integrate_expcubic(
            -a, -b, -c, -d,
            lower = all_sorted_data[i],
            upper = all_sorted_data[i + 1],
            minus_const = minus_const
        )

    }

    # integrate over [X_{(m)}, Inf)
    data_all_diff <- diff(all_sorted_data)
    g_b_opt <- sum(opt_theta[1:(length(opt_theta) - 1)] * data_all_diff + opt_theta[2:length(opt_theta)] * data_all_diff) / 2
    c_term <- g_a_opt + g_b_opt

    d_term1 <- sum((opt_theta[1:(length(opt_theta) - 1)] + opt_theta[2:length(opt_theta)]) * data_all_diff *
                       all_sorted_data[2:length(all_sorted_data)]) / 2
    d_term2 <- sum((opt_theta[2:length(opt_theta)] / 6 -
                        opt_theta[1:(length(opt_theta) - 1)] / 6 +
                        opt_theta[1:(length(opt_theta) - 1)] / 2) * data_all_diff ** 2)

    d_term <- -g_a_opt * domain[1] - d_term1 + d_term2

    result <- result + integrate_expcubic(
        0, 0, -c_term, -d_term,
        lower = all_sorted_data[length(all_sorted_data)],
        upper = Inf,
        minus_const = minus_const
        )

    return(result)

}
