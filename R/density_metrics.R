#' Computation of Various Distances between Two Probability Distributions
#'
#' A collection of functions to compute various distances between two probability distributions,
#' which are used to assess the quality of the density estimate.
#'
#' @param true_density An \code{R} object of class "truncated_normal", "truncated_gamma",
#' "truncated_lognormal" or "beta_dist", returned from the \code{truncated_normal},
#' \code{truncated_gamma}, \code{truncated_lognormal}, or \code{beta_dist}, respectively.
#' @param density_estimate An object of class "LogConcaveDESM",
#' usually the output of \code{\link{lcd_scorematching}} or
#' \code{\link{cv_optimal_density_estimate}}.
#' @param minus_const A numeric to be subtracted in the exponent to
#' ensure the finite-ness of the integration result. Default is \code{0}.
#' @param mc_batch_size A numeric to specify the batch size of each Monte Carlo draw of random samples.
#' Default is \code{100}.
#' @param mc_rel_tol_param A numeric to specify the criterion of terminating the batch Monte Carlo algorithm.
#' Default is \code{1e-2}.
#' @param print_error A logical value to specify whether to print the error after each Monte Carlo draw.
#' Default is \code{FALSE}.
#'
#' @return The approximate value of the probability distribution distance.
#'
#' @name density_metrics
NULL

#' @rdname density_metrics
#' @export
#'
kl_div <- function(true_density, density_estimate, minus_const = 0,
                   mc_batch_size = 100, mc_rel_tol_param = 1e-2, print_error = FALSE) {

    # check the domain
    true_domain <- true_density$domain
    esti_domain <- density_estimate$domain
    if ((true_domain[1] != esti_domain[1]) || (true_domain[2] != esti_domain[2])) {
        stop("The domains in true_density and density_estimate should be the same.")
    }

    # ----------------------------------------------------------------------------------------------
    batch <- true_density$sampling(mc_batch_size)
    true_vals <- log(true_density$evaluate_density(batch)$denvals)
    esti_vals <- log(evaluate_density(density_estimate, batch, minus_const)$density_vals)
    result1 <- mean(true_vals - esti_vals)
    rel_error <- .Machine$double.xmax
    batch_cnt <- 1
    while (rel_error > mc_rel_tol_param) {

        result2 <- result1

        new_batch <- true_density$sampling(mc_batch_size)
        true_vals <- log(true_density$evaluate_density(new_batch)$denvals)
        esti_vals <- log(evaluate_density(density_estimate, new_batch, minus_const)$density_vals)

        result1 <- mean(true_vals - esti_vals)

        result1 <- (result1 + result2) / 2
        rel_error <- abs((result1 - result2) / result2)

        batch_cnt <- batch_cnt + 1

        if (print_error) {
            message(paste0("Error of Batch ", batch_cnt, ": ", rel_error, "."))
        }

    }

    return(result1)

}

#' @rdname density_metrics
#' @export
#'
hyvarinen_div <- function(true_density, density_estimate,
                          mc_batch_size = 100, mc_rel_tol_param = 1e-2, print_error = FALSE) {

    # check the domain
    true_domain <- true_density$domain
    esti_domain <- density_estimate$domain
    if ((true_domain[1] != esti_domain[1]) || (true_domain[2] != esti_domain[2])) {
        stop("The domains in true_density and density_estimate should be the same.")
    }

    # ----------------------------------------------------------------------------------------------
    batch <- true_density$sampling(mc_batch_size)
    true_vals <- log(true_density$evaluate_logderiv1(batch)$logdervals)
    esti_vals <- log(evaluate_logden_deriv1(density_estimate, batch)$logderiv1_vals)

    result1 <- mean((true_vals - esti_vals) ** 2)
    rel_error <- .Machine$double.xmax
    batch_cnt <- 1
    while (rel_error > mc_rel_tol_param) {

        result2 <- result1

        new_batch <- true_density$sampling(mc_batch_size)
        true_vals <- log(true_density$evaluate_logderiv1(new_batch)$logdervals)
        esti_vals <- log(evaluate_logden_deriv1(density_estimate, new_batch)$logderiv1_vals)

        result1 <- mean((true_vals - esti_vals) ** 2)
        result1 <- (result1 + result2) / 2
        rel_error <- abs((result1 - result2) / result2)

        batch_cnt <- batch_cnt + 1

        if (print_error) {
            message(paste0("Error of Batch ", batch_cnt, ": ", rel_error, "."))
        }

    }

    return(result1)

}

#' @rdname density_metrics
#' @export
#'
tv_dist <- function(true_density, density_estimate, minus_const = 0,
                    mc_batch_size = 100, mc_rel_tol_param = 1e-2, print_error = FALSE) {

    # check the domain
    true_domain <- true_density$domain
    esti_domain <- density_estimate$domain
    if ((true_domain[1] != esti_domain[1]) || (true_domain[2] != esti_domain[2])) {
        stop("The domains in true_density and density_estimate should be the same.")
    }
    a <- true_domain[1]
    b <- true_domain[2]

    # ----------------------------------------------------------------------------------------------
    batch <- runif(mc_batch_size, a, b)
    true_vals <- true_density$evaluate_density(batch)$denvals
    esti_vals <- evaluate_density(density_estimate, batch, minus_const)$density_vals
    result1 <- mean(abs(true_vals - esti_vals) * (b - a))
    rel_error <- .Machine$double.xmax
    batch_cnt <- 1
    while (rel_error > mc_rel_tol_param) {

        result2 <- result1

        new_batch <- runif(mc_batch_size, a, b)
        true_vals <- true_density$evaluate_density(new_batch)$denvals
        esti_vals <- evaluate_density(density_estimate, new_batch, minus_const)$density_vals
        result1 <- mean(abs(true_vals - esti_vals) * (b - a))

        result1 <- (result1 + result2) / 2
        rel_error <- abs((result1 - result2) / result2)

        batch_cnt <- batch_cnt + 1

        if (print_error) {
            message(paste0("Error of Batch ", batch_cnt, ": ", rel_error, "."))
        }

    }

    return(result1)

}

#' @rdname density_metrics
#' @export
#'
hellinger_dist <- function(true_density, density_estimate, minus_const = 0,
                           mc_batch_size = 100, mc_rel_tol_param = 1e-2, print_error = FALSE) {

    # check the domain
    true_domain <- true_density$domain
    esti_domain <- density_estimate$domain
    if ((true_domain[1] != esti_domain[1]) || (true_domain[2] != esti_domain[2])) {
        stop("The domains in true_density and density_estimate should be the same.")
    }
    a <- true_domain[1]
    b <- true_domain[2]

    # ----------------------------------------------------------------------------------------------
    batch <- runif(mc_batch_size, a, b)
    true_vals <- true_density$evaluate_density(batch)$denvals
    esti_vals <- evaluate_density(density_estimate, batch, minus_const)$density_vals
    result1 <- mean((sqrt(true_vals) - sqrt(esti_vals)) ** 2 * (b - a)) / 2
    rel_error <- .Machine$double.xmax
    batch_cnt <- 1
    while (rel_error > mc_rel_tol_param) {

        result2 <- result1

        new_batch <- runif(mc_batch_size, a, b)
        true_vals <- true_density$evaluate_density(new_batch)$denvals
        esti_vals <- evaluate_density(density_estimate, new_batch, minus_const)$density_vals
        result1 <- mean((sqrt(true_vals) - sqrt(esti_vals)) ** 2 * (b - a)) / 2

        result1 <- (result1 + result2) / 2
        rel_error <- abs((result1 - result2) / result2)

        batch_cnt <- batch_cnt + 1

        if (print_error) {
            message(paste0("Error of Batch ", batch_cnt, ": ", rel_error, "."))
        }

    }

    return(result1)

}
