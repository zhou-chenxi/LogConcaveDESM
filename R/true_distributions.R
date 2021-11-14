#' True Distributions
#'
#' Functions for drawing random samples from the specified distribution and
#' evaluating the underlying density function and the first derivative of its logarithm.
#'
#' @import ggplot2
#' @import dplyr
#'
#' @return An object of class "truncated_normal", "truncated_gamma", "truncated_lognormal",
#' or "beta_dist", whose underlying structure is a list containing the following elements
#' \describe{
#'     \item{sampling}{function: generates random samples from the underlying distribution.
#'     The input is \code{sample_size}, the number of desired samples, and the output is
#'     a numeric vector of random samples. }
#'     \item{evaluate_density}{function: returns the underlying density
#'     values. The input is \code{newx}, the points at which the density function
#'     is to be evaluated, and the output is a data frame with the first column being
#'     the sorted \code{newx} and the second column being the corresponding density values.}
#'     \item{evaluate_logderiv1}{function: returns the first derivative values
#'     of the logarithm of the underlying density function.
#'     The input is \code{newx}, the points at which the first derivative
#'     is to be evaluated, and the output is a data frame with the first column being
#'     the sorted \code{newx} and the second column being the corresponding first derivative
#'     values.}
#'     \item{domain}{numeric: a numeric vector over which the density function is defined.}
#' }
#'
#' @name true_dist
NULL

#' @rdname true_dist
#' @export
#'
#' @param parent_mean A numeric of the mean parameter in the parent normal distribution.
#' @param parent_sd A numeric of the standard deviation parameter in the parent normal distribution;
#' must be strictly positive.
#' @param domain A numeric vector of the domain over which the distribution is defined.
#'
#' @examples
#' library(ggplot2)
#' ######################################
#' # truncated normal
#' parent_mean <- 4
#' parent_sd <- 1
#' domain <- c(2, 5)
#' t_normal <- truncated_normal(parent_mean, parent_sd, domain)
#'
#' # ------------------------------------
#' # sampling function
#' samples <- t_normal$sampling(10000)
#' df <- data.frame(x = samples)
#' x <- seq(domain[1], domain[2], by = 0.01)
#' y <- (dnorm(x, parent_mean, parent_sd) /
#' (pnorm(domain[2], parent_mean, parent_sd) -
#' pnorm(domain[1], parent_mean, parent_sd)))
#' density_df <- data.frame(x = x, y = y)
#' # use the Freedman-Diaconis rule for the binwidth
#' binwidth <- 2 * IQR(samples) / (length(samples) ** (1/3))
#' plot <- ggplot() +
#' geom_histogram(data = df, aes(x = x, y = ..density..),
#' binwidth = binwidth, color = "darkblue", fill = "lightblue") +
#' geom_line(data = density_df, aes(x = x, y = y), color = "black", size = 0.75) +
#' theme_bw()
#' plot
#'
#' # ------------------------------------
#' # evaluate density
#' t_normal$evaluate_density(seq(-5, 10, by = 0.1))
#'
#' # ------------------------------------
#' # evaluate derivative of log-density
#' t_normal$evaluate_logderiv1(seq(-5, 10, by = 0.1))
#'
truncated_normal <- function(parent_mean, parent_sd, domain) {

    # checking
    if (length(domain) != 2) {
        stop("The length of the domain should be 2.")
    }

    if (domain[1] > domain[2]) {
        stop("Numbers in the domain should be in the ascending order.")
    }

    if (parent_sd <= 0) {
        stop("parent_sd must be strictly positive.")
    }

    sampling <- function(sample_size) {

        unif_samples <- runif(sample_size)
        intermediate <- (
            pnorm(domain[1],
                  mean = parent_mean,
                  sd = parent_sd) +
                unif_samples * (
                    pnorm(domain[2],
                          mean = parent_mean,
                          sd = parent_sd) -
                        pnorm(domain[1],
                              mean = parent_mean,
                              sd = parent_sd)))

        result <- qnorm(intermediate, mean = parent_mean, sd = parent_sd)

        return(result)
    }

    evaluate_density <- function(newx) {

        valid_newx <- newx[(newx >= domain[1]) & (newx <= domain[2])]
        invalid_newx <- newx[(newx < domain[1]) | (newx > domain[2])]

        valid_denvals <- (
            dnorm(valid_newx,
                  mean = parent_mean,
                  sd = parent_sd) / (
                      pnorm(domain[2],
                            mean = parent_mean,
                            sd = parent_sd) -
                          pnorm(domain[1],
                                mean = parent_mean,
                                sd = parent_sd))
        )
        valid_result <- data.frame(newx = valid_newx, denvals = valid_denvals)

        invalid_denvals <- rep(0, length(invalid_newx))
        invalid_result <- data.frame(newx = invalid_newx, denvals = invalid_denvals)

        all_result <- dplyr::bind_rows(valid_result, invalid_result)
        all_result <- dplyr::arrange(all_result, newx)

        return(all_result)
    }

    evaluate_logderiv1 <- function(newx) {

        valid_newx <- newx[(newx >= domain[1]) & (newx <= domain[2])]
        invalid_newx <- newx[(newx < domain[1]) | (newx > domain[2])]

        valid_denvals <- -(valid_newx - parent_mean) / (parent_sd ** 2)
        valid_result <- data.frame(newx = valid_newx, logdervals = valid_denvals)

        invalid_denvals <- rep(NA, length(invalid_newx))
        invalid_result <- data.frame(newx = invalid_newx, logdervals = invalid_denvals)

        all_result <- dplyr::bind_rows(valid_result, invalid_result)
        all_result <- dplyr::arrange(all_result, newx)

        return(all_result)
    }

    obj <- list(
        sampling = sampling,
        evaluate_density = evaluate_density,
        evaluate_logderiv1 = evaluate_logderiv1,
        domain = domain
    )

    class(obj) <- "truncated_normal"

    return(obj)

}

#' @rdname true_dist
#' @export
#'
#' @param parent_shape A numeric of the shape parameter in the parent gamma distribution;
#' must be strictly positive.
#' @param parent_rate A numeric of the rate parameter in the parent gamma distribution;
#' must be strictly positive.
#'
#' @examples
#' ######################################
#' # truncated gamma
#' domain <- c(0, 10)
#' parent_shape <- 5
#' parent_rate <- 1
#' t_gamma <- truncated_gamma(parent_shape, parent_rate, domain)
#'
#' # ------------------------------------
#' # sampling function
#' samples <- t_gamma$sampling(10000)
#' df <- data.frame(x = samples)
#' x <- seq(domain[1], domain[2], by = 0.01)
#' y <- (dgamma(x, parent_shape, parent_rate) /
#' (pgamma(domain[2], parent_shape, parent_rate) -
#' pgamma(domain[1], parent_shape, parent_rate)))
#' density_df <- data.frame(x = x, y = y)
#' # use the Freedman-Diaconis rule for the binwidth
#' binwidth <- 2 * IQR(samples) / (length(samples) ** (1/3))
#' plot <- ggplot() +
#' geom_histogram(data = df, aes(x = x, y = ..density..),
#' binwidth = binwidth, color = "darkblue", fill = "lightblue") +
#' geom_line(data = density_df, aes(x = x, y = y), color = "black", size = 0.75) +
#' theme_bw()
#' plot
#'
#' # ------------------------------------
#' # evaluate density
#' t_gamma$evaluate_density(seq(-5, 10, by = 0.1))
#'
#' # ------------------------------------
#' # evaluate derivative of log-density
#' t_gamma$evaluate_logderiv1(seq(-5, 10, by = 0.1))
#'
truncated_gamma <- function(parent_shape, parent_rate, domain) {

    # checking
    if (length(domain) != 2) {
        stop("The length of the domain should be 2.")
    }

    domain <- sort(domain)

    if (any(domain < 0)) {
        stop("The domain must be a subset of the positive real line.")
    }

    if (parent_rate <= 0) {
        stop("parent_rate must be strictly positive.")
    }

    if (parent_shape <= 0) {
        stop("parent_shape must be strictly positive.")
    }

    sampling <- function(sample_size) {
        unif_samples <- runif(sample_size)
        intermediate <- (
            pgamma(domain[1],
                   shape = parent_shape,
                   rate = parent_rate) +
                unif_samples * (
                    pgamma(domain[2],
                           shape = parent_shape,
                           rate = parent_rate) -
                        pgamma(domain[1],
                               shape = parent_shape,
                               rate = parent_rate)))

        result <- qgamma(intermediate, shape = parent_shape, rate = parent_rate)

        return(result)
    }

    evaluate_density  <- function(newx) {

        valid_newx <- newx[(newx >= domain[1]) & (newx <= domain[2])]
        invalid_newx <- newx[(newx < domain[1]) | (newx > domain[2])]

        valid_denvals <- (
            dgamma(valid_newx,
                   shape = parent_shape,
                   rate = parent_rate) / (
                       pgamma(domain[2],
                              shape = parent_shape,
                              rate = parent_rate) -
                           pgamma(domain[1],
                                  shape = parent_shape,
                                  rate = parent_rate))
        )
        valid_result <- data.frame(newx = valid_newx, denvals = valid_denvals)

        invalid_denvals <- rep(0, length(invalid_newx))
        invalid_result <- data.frame(newx = invalid_newx, denvals = invalid_denvals)

        all_result <- dplyr::bind_rows(valid_result, invalid_result)
        all_result <- dplyr::arrange(all_result, newx)

        return(all_result)
    }

    evaluate_logderiv1 <- function(newx) {

        valid_newx <- newx[(newx >= domain[1]) & (newx <= domain[2])]
        invalid_newx <- newx[(newx < domain[1]) | (newx > domain[2])]

        valid_denvals <- (parent_shape - 1) / valid_newx - parent_rate
        valid_result <- data.frame(newx = valid_newx, logdervals = valid_denvals)

        invalid_denvals <- rep(NA, length(invalid_newx))
        invalid_result <- data.frame(newx = invalid_newx, logdervals = invalid_denvals)

        all_result <- dplyr::bind_rows(valid_result, invalid_result)
        all_result <- dplyr::arrange(all_result, newx)

        return(all_result)
    }

    obj <- list(
        sampling = sampling,
        evaluate_density = evaluate_density,
        evaluate_logderiv1 = evaluate_logderiv1,
        domain = domain
    )

    class(obj) <- "truncated_gamma"

    return(obj)

}


#' @rdname true_dist
#' @export
#'
#' @param parent_meanlog A numeric of the mean parameter in the parent log-normal distribution
#' in the log scale.
#' @param parent_sdlog A numeric of the standard deviation parameter in
#' the parent log-normal distribution in the log scale.
#'
#' @examples
#' ######################################
#' # truncated lognormal
#' parent_meanlog <- 0
#' parent_sdlog <- 1/4
#' domain <- c(0, 5)
#' t_lognormal <- truncated_lognormal(parent_meanlog, parent_sdlog, domain)
#'
#' # ------------------------------------
#' # sampling function
#' samples <- t_lognormal$sampling(10000)
#' df <- data.frame(x = samples)
#' x <- seq(domain[1], domain[2], by = 0.01)
#' y <- (dlnorm(x, parent_meanlog, parent_sdlog) /
#' (plnorm(domain[2], parent_meanlog, parent_sdlog) -
#' plnorm(domain[1], parent_meanlog, parent_sdlog)))
#' density_df <- data.frame(x = x, y = y)
#' # use the Freedman-Diaconis rule for the binwidth
#' binwidth <- 2 * IQR(samples) / (length(samples) ** (1/3))
#' plot <- ggplot() +
#' geom_histogram(data = df, aes(x = x, y = ..density..),
#' binwidth = binwidth, color = "darkblue", fill = "lightblue") +
#' geom_line(data = density_df, aes(x = x, y = y), color = "black", size = 0.75) +
#' theme_bw()
#' plot
#'
#' # ------------------------------------
#' # evaluate density
#' t_lognormal$evaluate_density(seq(-5, 10, by = 0.1))
#'
#' # ------------------------------------
#' # evaluate derivative of log-density
#' t_lognormal$evaluate_logderiv1(seq(-5, 10, by = 0.1))
#'
truncated_lognormal <- function(parent_meanlog, parent_sdlog, domain) {

    # checking
    if ((Inf %in% domain) || (-Inf %in% domain)) {
        stop("The domain must be bounded.")
    }

    if (length(domain) != 2) {
        stop("The length of the domain should be 2.")
    }

    if (domain[1] > domain[2]) {
        stop("Numbers in the domain should be in the ascending order.")
    }

    if (parent_sdlog <= 0) {
        stop("parent_sd must be strictly positive.")
    }

    sampling <- function(sample_size) {
        unif_samples <- runif(sample_size)
        intermediate <- (
            plnorm(domain[1],
                   meanlog = parent_meanlog,
                   sdlog = parent_sdlog) +
                unif_samples * (
                    plnorm(domain[2],
                           meanlog = parent_meanlog,
                           sdlog = parent_sdlog) -
                        plnorm(domain[1],
                               meanlog = parent_meanlog,
                               sdlog = parent_sdlog)))

        result <- qlnorm(intermediate,
                         meanlog = parent_meanlog,
                         sdlog = parent_sdlog)

        return(result)
    }

    evaluate_density <- function(newx) {

        valid_newx <- newx[(newx >= domain[1]) & (newx <= domain[2])]
        invalid_newx <- newx[(newx < domain[1]) | (newx > domain[2])]

        valid_denvals <- (
            dlnorm(valid_newx,
                   meanlog = parent_meanlog,
                   sdlog = parent_sdlog) / (
                       plnorm(domain[2],
                              meanlog = parent_meanlog,
                              sdlog = parent_sdlog) -
                           plnorm(domain[1],
                                  meanlog = parent_meanlog,
                                  sdlog = parent_sdlog))
        )
        valid_result <- data.frame(newx = valid_newx, denvals = valid_denvals)

        invalid_denvals <- rep(0, length(invalid_newx))
        invalid_result <- data.frame(newx = invalid_newx, denvals = invalid_denvals)

        all_result <- dplyr::bind_rows(valid_result, invalid_result)
        all_result <- dplyr::arrange(all_result, newx)

        return(all_result)

    }

    evaluate_logderiv1 <- function(newx) {

        valid_newx <- newx[(newx >= domain[1]) & (newx <= domain[2])]
        invalid_newx <- newx[(newx < domain[1]) | (newx > domain[2])]

        valid_denvals <- (-1 / valid_newx -
                              (log(valid_newx) - parent_meanlog) / parent_sdlog ** 2 / valid_newx)
        valid_result <- data.frame(newx = valid_newx, logdervals = valid_denvals)

        invalid_denvals <- rep(NA, length(invalid_newx))
        invalid_result <- data.frame(newx = invalid_newx, logdervals = invalid_denvals)

        all_result <- dplyr::bind_rows(valid_result, invalid_result)
        all_result <- dplyr::arrange(all_result, newx)

        return(all_result)
    }

    obj <- list(
        sampling = sampling,
        evaluate_density = evaluate_density,
        evaluate_logderiv1 = evaluate_logderiv1,
        domain = domain
    )

    class(obj) <- "truncated_lognormal"

    return(obj)

}


#' @rdname true_dist
#' @export
#'
#' @param shape1,shape2 Numeric values of the parameters in the parent beta distribution;
#' must be strictly positive.
#'
#' @examples
#' ######################################
#' # beta distribution
#' shape1 <- 4
#' shape2 <- 2
#' domain <- c(0, 1)
#' b <- beta_dist(shape1, shape2)
#'
#' # ------------------------------------
#' # sampling function
#' samples <- b$sampling(10000)
#' df <- data.frame(x = samples)
#' x <- seq(domain[1], domain[2], by = 0.01)
#' y <- dbeta(x, shape1, shape2)
#' density_df <- data.frame(x = x, y = y)
#' # use the Freedman-Diaconis rule for the binwidth
#' binwidth <- 2 * IQR(samples) / (length(samples) ** (1/3))
#' plot <- ggplot() +
#' geom_histogram(data = df, aes(x = x, y = ..density..),
#' binwidth = binwidth, color = "darkblue", fill = "lightblue") +
#' geom_line(data = density_df, aes(x = x, y = y), color = "black", size = 0.75) +
#' theme_bw()
#' plot
#'
#' # ------------------------------------
#' # evaluate density
#' b$evaluate_density(seq(0, 1, by = 0.1))
#'
#' # ------------------------------------
#' # evaluate derivative of log-density
#' b$evaluate_logderiv1(seq(0, 1, by = 0.1))
#'
beta_dist <- function(shape1, shape2) {

    domain <- c(0, 1)

    if ((shape1 <= 0) | (shape2 <= 0)) {
        stop("Both shape parameters shape1 and shape2 must be strictly positive.")
    }

    sampling <- function(sample_size) {
        result <- rbeta(sample_size,
                        shape1 = shape1,
                        shape2 = shape2)

        return(result)
    }

    evaluate_density <- function(newx) {

        valid_newx <- newx[(newx >= domain[1]) & (newx <= domain[2])]
        invalid_newx <- newx[(newx < domain[1]) | (newx > domain[2])]

        valid_denvals <- dbeta(valid_newx,
                               shape1 = shape1,
                               shape2 = shape2)
        valid_result <- data.frame(newx = valid_newx, denvals = valid_denvals)

        invalid_denvals <- rep(0, length(invalid_newx))
        invalid_result <- data.frame(newx = invalid_newx, denvals = invalid_denvals)

        all_result <- dplyr::bind_rows(valid_result, invalid_result)
        all_result <- dplyr::arrange(all_result, newx)

        return(all_result)
    }

    evaluate_logderiv1 <- function(newx) {

        valid_newx <- newx[(newx >= domain[1]) & (newx <= domain[2])]
        invalid_newx <- newx[(newx < domain[1]) | (newx > domain[2])]

        valid_denvals <- shape1 / valid_newx - shape2 / (1 - valid_newx)
        valid_result <- data.frame(newx = valid_newx, logdervals = valid_denvals)

        invalid_denvals <- rep(NA, length(invalid_newx))
        invalid_result <- data.frame(newx = invalid_newx, logdervals = invalid_denvals)

        all_result <- dplyr::bind_rows(valid_result, invalid_result)
        all_result <- dplyr::arrange(all_result, newx)

        return(all_result)
    }

    obj <- list(
        sampling = sampling,
        evaluate_density = evaluate_density,
        evaluate_logderiv1 = evaluate_logderiv1,
        domain = domain
    )

    class(obj) <- "beta"

    return(obj)

}
