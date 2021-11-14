#' Computation of Second Derivative of Score Matching Log-concave Density Estimate
#'
#' Computes the second derivative of the logarithm of the (penalized) log-concave score matching
#' density estimate based on i.i.d samples, assuming it is a continuous piecewise linear function
#' with knots at the samples.
#' The nature of the underlying density estimation problem is a constrained quadratic optimization problem,
#' which is solved by \code{CVXR} package.
#' The output is an object of class "\code{LogConcaveDESM}".
#'
#' @details The functions \code{lcd_scorematching_bounded}, \code{lcd_scorematching_R}, \code{lcd_scorematching_ninfb} and
#' \code{lcd_scorematching_ainf} compute the second derivative of the logarithm of the (penalized) log-concave score matching density estimate
#' at \code{data} when the underlying \code{domain} is a bounded interval, the entire real line,
#' an interval of the form \eqn{-\infty, b} for some \eqn{b < \infty}, and
#' an interval of the form \eqn{a, \infty} for some \eqn{a > -\infty}, respectively.
#' The function \code{lcd_scorematching} encompasses all four cases.
#'
#' @param data A numeric vector whose log-concave density function is to be estimated;
#' missing values are automatically removed.
#' @param domain A numeric vector of length 2 specifying the left and right
#' endpoints of the domain. The \code{domain} can be a bounded interval, the entire real line,
#' a left-bounded and right-unbounded interval, or a left-unbounded and right-bounded interval.
#' Its components cannot be \code{NaN}.
#' @param penalty_param Penalty parameter for computing the density estimate; must be non-negative.
#' Default is \code{1e-1}.
#' @param verbose An optional logical value to indicate whether to print detailed CVXR
#' solver output. Default is \code{FALSE}.
#'
#' @import CVXR
#' @import plyr
#'
#' @return An object of class "LogConcaveDESM" whose underlying structure is
#' a list containing the following elements
#' \describe{
#'     \item{call}{call: the call which produced the result.}
#'     \item{data}{numeric: the original data whose log-concave density function is to be estimated.}
#'     \item{sorted_unique_data}{numeric: the sorted unique data.}
#'     \item{data_weights}{numeric: weight vector of the sorted_unique_data.}
#'     \item{domain}{numeric: the bounded domain over the density function is assumed and estimated.}
#'     \item{opt_theta}{numeric: the optimal second derivative of the logarithm of the (penalized)
#'     log-concave score matching density estimate at the sorted unique data and two boundary points.}
#'     \item{matrix_A}{numeric: matrix A used to compute the log-concave score matching density estimate.}
#'     \item{matrix_B}{numeric: matrix B used to compute the log-concave score matching density estimate.}
#'     \item{penalty_param}{numeric: the penalty parameter used to compute the log-concave
#'     score matching density estimate.}
#' }
#'
#'
#' @examples
#' set.seed(1119)
#' N <- 100
#' data <- rnorm(N)
#' domain <- c(-5, 5)
#' # no penalty term
#' result <- lcd_scorematching(data, domain, penalty_param = 0)
#'
#' # with penalty term
#' result <- lcd_scorematching(data, domain, penalty_param = 1e-10)
#'
#' @name lcd
NULL

#' @rdname lcd
#' @export
lcd_scorematching <- function(data, domain, penalty_param = 1e-1, verbose = FALSE) {

    # preprocess data
    data <- as.numeric(data)
    data <- data[!is.na(data)]

    # check the validity of the penalty parameter
    stopifnot(penalty_param >= 0)

    # check the length of domain is 2
    stopifnot(length(domain) == 2)

    # check the validity of the domain
    if (any(is.nan(domain))) {
        stop("domain cannot contain 'NaN'.")
    }

    domain <- sort(domain)

    domain1 <- domain[1]
    domain2 <- domain[2]

    if ((domain1 == -Inf) & (domain2 == Inf)) {

        # R case
        result <- lcd_scorematching_R(
            data = data,
            domain = domain,
            penalty_param = penalty_param,
            verbose = verbose)

    } else if (is.finite(domain1) & is.finite(domain2)) {

        # bounded interval case
        result <- lcd_scorematching_bounded(
            data = data,
            domain = domain,
            penalty_param = penalty_param,
            verbose = verbose)

    } else if (is.finite(domain1) & (domain2 == Inf)) {

        # [a, Inf) case
        result <- lcd_scorematching_ainf(
            data = data,
            domain = domain,
            penalty_param = penalty_param,
            verbose = verbose)

    }  else if ((domain1 == -Inf) & is.finite(domain2)) {

        # (-Inf, b] case
        result <- lcd_scorematching_ninfb(
            data = data,
            domain = domain,
            penalty_param = penalty_param,
            verbose = verbose)

    } else {

        stop(paste0("The domain entered, ", domain, ", is not valid."))

    }

    class(result) <- 'LogConcaveDESM'

    return(result)

}


#' @rdname lcd
#' @export
lcd_scorematching_bounded <- function(data, domain, penalty_param = 1e-1, verbose = FALSE) {

    if ((Inf %in% domain) || (-Inf %in% domain)) {
        stop("domain cannot contain '-Inf' or 'Inf'.")
    }

    # -----------------------------------------------------------------------------
    # prepare data
    N <- length(data)
    sorted_data <- sort(data)
    freq_data <- plyr::count(sorted_data)
    unique_data <- freq_data$x
    data_weight <- freq_data$freq
    no_unique <- length(unique_data)

    # -----------------------------------------------------------------------------
    # build vectors A_i and B_i
    lower_bound <- domain[1]
    upper_bound <- domain[2]
    if (lower_bound > min(data) || upper_bound < max(data)) {
        stop('There exist data outside of the domain specified.')
    }

    diff_data <- diff(c(lower_bound, unique_data))
    A_i <- c(diff_data, 0, 0)

    matrix_A <- matrix(rep(A_i, no_unique), ncol = no_unique, byrow = FALSE)
    matrix_A[row(matrix_A) - col(matrix_A) >= 1] <- 0

    matrix_B <- rbind(
        matrix(0, nrow = 1, ncol = ncol(matrix_A)),
        matrix_A[1:(nrow(matrix_A) - 1), ]
    )

    # -----------------------------------------------------------------------------
    # form the weight matrix and vector
    weight_matrix <- base::diag(data_weight) / N
    inner_weight_matrix <- (weight_matrix -
                                matrix(data_weight, ncol = 1) %*% matrix(data_weight, nrow = 1) / N ** 2)
    weight_vector <- c(0, data_weight, 0) / N

    # -----------------------------------------------------------------------------
    # build the quadratic matrix for optimization
    quad_matrix <- 1/8 * (matrix_A + matrix_B) %*% inner_weight_matrix %*% t(matrix_A + matrix_B)

    # -----------------------------------------------------------------------------
    # construct the penalization matrix
    diff_data_all <- diff(c(lower_bound, unique_data, upper_bound))
    pen_matrix <- matrix(0, no_unique + 2, no_unique + 2)

    # diagonal elements
    diag_elements <- c(
        diff_data_all[1],
        diff(c(lower_bound, unique_data, upper_bound), lag = 2),
        diff_data_all[length(diff_data_all)]) / 3
    diag(pen_matrix) <- diag_elements

    # off-diagonal elements
    D_i <- c(diff_data_all) / 6
    pen_matrix[row(pen_matrix) - col(pen_matrix) == 1] <- D_i
    pen_matrix[row(pen_matrix) - col(pen_matrix) == -1] <- D_i

    # -----------------------------------------------------------------------------
    # solve the optimization problem
    theta_vec <- CVXR::Variable(no_unique + 2)
    objective <- CVXR::Minimize(
        CVXR::quad_form(theta_vec, quad_matrix) -
            t(theta_vec) %*% weight_vector +
            penalty_param * CVXR::quad_form(theta_vec, pen_matrix) / 2
    )
    if (penalty_param == 0) {
        constraint <- list(theta_vec >= 0, theta_vec[length(theta_vec)] == 0)
    } else {
        constraint <- list(theta_vec >= 0)
    }

    problem <- CVXR::Problem(objective, constraints = constraint)
    result <- CVXR::psolve(problem, verbose = verbose)

    message(paste0("The status of solving the constrained quadratic optimization problem is: ", result$status, ". "))

    # -----------------------------------------------------------------------------
    # build a class of results
    result <- list(
        call = match.call(),
        data = data,
        sorted_unique_data = sorted_data,
        data_weights = data_weight / length(data),
        domain = domain,
        opt_theta = result$getValue(theta_vec),
        matrix_A = matrix_A,
        matrix_B = matrix_B,
        penalty_param = penalty_param
    )

    return(result)
}

#' @rdname lcd
#' @export
lcd_scorematching_R <- function(data, domain, penalty_param = 1e-1, verbose = FALSE) {

    # double check the domain
    stopifnot(domain == c(-Inf, Inf))

    sorted_data <- sort(data)
    freq_data <- plyr::count(sorted_data)
    unique_data <- freq_data$x
    data_weight <- freq_data$freq
    no_unique <- length(unique_data)
    N <- length(data)

    # build vectors A_i and B_i
    diff_data <- diff(unique_data)
    A_i <- c(diff_data, 0)
    matrix_A <- matrix(rep(A_i, no_unique - 1), ncol = no_unique - 1, byrow = FALSE)
    matrix_A[row(matrix_A) - col(matrix_A) >= 1] <- 0

    matrix_B <- rbind(
        matrix(0, nrow = 1, ncol = ncol(matrix_A)),
        matrix_A[1:(nrow(matrix_A) - 1), ]
    )

    # form the weight matrix W = diag(w2, w3, ..., wm) / n
    # an (m-1) * (m-1) matrix
    weight_matrix <- base::diag(data_weight[2:length(data_weight)]) / N
    # inner_weight_matrix1 = W 1_m 1_m^\top W
    inner_weight_matrix1 <- weight_matrix %*% matrix(1, nrow = no_unique - 1, ncol = 1) %*%
        matrix(1, nrow = 1, ncol = no_unique - 1) %*% weight_matrix
    inner_weight_matrix <- weight_matrix - inner_weight_matrix1 + data_weight[1] / N * inner_weight_matrix1

    # form the weight vector w = (w1, w2, w3, ..., wm) / n
    # an R^m vector
    weight_vector <- data_weight / N

    # --------------------------------------------------------------------------------------------------
    # create the matrix for the penalty term
    diff_data_all <- diff(unique_data)
    pen_matrix <- matrix(0, no_unique, no_unique)

    # diagonal elements
    diag_elements <- c(
        diff_data_all[1],
        diff(unique_data, lag = 2),
        diff_data_all[length(diff_data_all)]) / 3
    diag(pen_matrix) <- diag_elements

    # off-diagonal elements
    D_i <- c(diff_data_all) / 6
    pen_matrix[row(pen_matrix) - col(pen_matrix) == 1] <- D_i
    pen_matrix[row(pen_matrix) - col(pen_matrix) == -1] <- D_i

    # ----------------------------------------------------------------------------------------------
    # solve the optimization problem
    theta_vec <- CVXR::Variable(no_unique)
    quad_matrix <- 1/8 * (matrix_A + matrix_B) %*% inner_weight_matrix %*% t(matrix_A + matrix_B)
    objective <- CVXR::Minimize(
        CVXR::quad_form(theta_vec, quad_matrix) -
            t(theta_vec) %*% weight_vector +
            penalty_param * CVXR::quad_form(theta_vec, pen_matrix) / 2
    )
    constraint <- list(theta_vec[1] == 0, theta_vec[length(theta_vec)] == 0, theta_vec >= 0)

    problem <- CVXR::Problem(objective, constraints = constraint)
    result <- CVXR::psolve(problem, verbose = verbose)

    message(paste0("The status of solving the constrained quadratic optimization problem is: ", result$status, ". "))

    # ----------------------------------------------------------------------------------------------
    # build a class of results
    result <- list(
        call = match.call(),
        data = data,
        sorted_unique_data = sorted_data,
        data_weights = data_weight / length(data),
        domain = domain,
        opt_theta = result$getValue(theta_vec),
        matrix_A = matrix_A,
        matrix_B = matrix_B,
        penalty_param = penalty_param
    )

    return(result)

}


#' @rdname lcd
#' @export
lcd_scorematching_ninfb <- function(data, domain, penalty_param = 1e-1, verbose = FALSE) {

    # double check the domain
    stopifnot(domain[1] == -Inf, is.finite(domain[2]))

    # check the data validity
    if (max(data) > domain[2]) {
        stop("There exist data outside the domain.")
    }

    sorted_data <- sort(data)
    freq_data <- plyr::count(sorted_data)
    unique_data <- freq_data$x
    data_weight <- freq_data$freq
    no_unique <- length(unique_data)
    N <- length(data)

    # build vectors A_i and B_i
    diff_data <- diff(unique_data)
    A_i <- c(diff_data, 0, 0)
    matrix_A <- matrix(rep(A_i, no_unique - 1), ncol = no_unique - 1, byrow = FALSE)
    matrix_A[row(matrix_A) - col(matrix_A) >= 1] <- 0

    matrix_B <- rbind(
        matrix(0, nrow = 1, ncol = ncol(matrix_A)),
        matrix_A[1:(nrow(matrix_A) - 1), ]
    )

    # form the weight matrix W = diag(w2, w3, ..., wm) / n
    # an (m-1) * (m-1) matrix
    weight_matrix <- base::diag(data_weight[2:length(data_weight)]) / N
    # inner_weight_matrix1 = W 1_m 1_m^\top W
    inner_weight_matrix1 <- weight_matrix %*% matrix(1, nrow = no_unique - 1, ncol = 1) %*%
        matrix(1, nrow = 1, ncol = no_unique - 1) %*% weight_matrix
    inner_weight_matrix <- weight_matrix - inner_weight_matrix1 + data_weight[1] / N * inner_weight_matrix1

    # form the weight vector w = (w1, w2, w3, ..., wm, 0) / n
    # an R^(m+1) vector
    weight_vector <- c(data_weight, 0) / N

    # --------------------------------------------------------------------------------------------------
    # create the matrix for the penalty term
    diff_data_all <- diff(c(unique_data, domain[2]))
    pen_matrix <- matrix(0, no_unique + 1, no_unique + 1)

    # diagonal elements
    diag_elements <- c(
        diff_data_all[1],
        diff(c(unique_data, domain[2]), lag = 2),
        diff_data_all[length(diff_data_all)]) / 3
    diag(pen_matrix) <- diag_elements

    # off-diagonal elements
    D_i <- c(diff_data_all) / 6
    pen_matrix[row(pen_matrix) - col(pen_matrix) == 1] <- D_i
    pen_matrix[row(pen_matrix) - col(pen_matrix) == -1] <- D_i

    # ----------------------------------------------------------------------------------------------
    # solve the optimization problem
    theta_vec <- CVXR::Variable(no_unique + 1)
    quad_matrix <- 1/8 * (matrix_A + matrix_B) %*% inner_weight_matrix %*% t(matrix_A + matrix_B)
    objective <- CVXR::Minimize(
        CVXR::quad_form(theta_vec, quad_matrix) -
            t(theta_vec) %*% weight_vector +
            penalty_param * CVXR::quad_form(theta_vec, pen_matrix) / 2
    )
    constraint <- list(theta_vec[1] == 0, theta_vec >= 0)

    problem <- CVXR::Problem(objective, constraints = constraint)
    result <- CVXR::psolve(problem, verbose = verbose)

    message(paste0("The status of solving the constrained quadratic optimization problem is: ", result$status, ". "))

    # ----------------------------------------------------------------------------------------------
    # build a class of results
    result <- list(
        call = match.call(),
        data = data,
        sorted_unique_data = sorted_data,
        data_weights = data_weight / length(data),
        domain = domain,
        opt_theta = result$getValue(theta_vec),
        matrix_A = matrix_A,
        matrix_B = matrix_B,
        penalty_param = penalty_param
    )

    return(result)

}


#' @rdname lcd
#' @export
lcd_scorematching_ainf <- function(data, domain, penalty_param = 1e-1, verbose = FALSE) {

    # double check the domain
    stopifnot(is.finite(domain[1]), domain[2] == Inf)

    # check the data validity
    if (min(data) < domain[1]) {
        stop("There exist data outside the domain.")
    }

    sorted_data <- sort(data)
    freq_data <- plyr::count(sorted_data)
    unique_data <- freq_data$x
    data_weight <- freq_data$freq
    no_unique <- length(unique_data)
    N <- length(data)

    # build vectors A_i and B_i
    diff_data <- diff(c(domain[1], unique_data))
    A_i <- c(diff_data, 0)
    matrix_A <- matrix(rep(A_i, no_unique), ncol = no_unique, byrow = FALSE)
    matrix_A[row(matrix_A) - col(matrix_A) >= 1] <- 0

    matrix_B <- rbind(
        matrix(0, nrow = 1, ncol = ncol(matrix_A)),
        matrix_A[1:(nrow(matrix_A) - 1), ]
    )

    # form the weight matrix W = diag(w1, w2, w3, ..., wm) / n
    # an m * m matrix
    weight_matrix <- base::diag(data_weight) / N
    # inner_weight_matrix1 = W 1_m 1_m^\top W
    inner_weight_matrix1 <- weight_matrix %*% matrix(1, nrow = no_unique, ncol = 1) %*%
        matrix(1, nrow = 1, ncol = no_unique) %*% weight_matrix
    inner_weight_matrix <- weight_matrix - inner_weight_matrix1

    # form the weight vector w = (0, w1, w2, w3, ..., wm) / n
    # an R^(m+1) vector
    weight_vector <- c(0, data_weight) / N

    # --------------------------------------------------------------------------------------------------
    # create the matrix for the penalty term
    diff_data_all <- diff(c(domain[1], unique_data))
    pen_matrix <- matrix(0, no_unique + 1, no_unique + 1)

    # diagonal elements
    diag_elements <- c(
        diff_data_all[1],
        diff(c(domain[1], unique_data), lag = 2),
        diff_data_all[length(diff_data_all)]) / 3
    diag(pen_matrix) <- diag_elements

    # off-diagonal elements
    D_i <- c(diff_data_all) / 6
    pen_matrix[row(pen_matrix) - col(pen_matrix) == 1] <- D_i
    pen_matrix[row(pen_matrix) - col(pen_matrix) == -1] <- D_i

    # ----------------------------------------------------------------------------------------------
    # solve the optimization problem
    theta_vec <- CVXR::Variable(no_unique + 1)
    quad_matrix <- 1/8 * (matrix_A + matrix_B) %*% inner_weight_matrix %*% t(matrix_A + matrix_B)
    objective <- CVXR::Minimize(
        CVXR::quad_form(theta_vec, quad_matrix) -
            t(theta_vec) %*% weight_vector +
            penalty_param * CVXR::quad_form(theta_vec, pen_matrix) / 2
    )
    constraint <- list(theta_vec[length(theta_vec)] == 0, theta_vec >= 0)

    problem <- CVXR::Problem(objective, constraints = constraint)
    result <- CVXR::psolve(problem, verbose = verbose)

    message(paste0("The status of solving the constrained quadratic optimization problem is: ", result$status, ". "))

    # ----------------------------------------------------------------------------------------------
    # build a class of results
    result <- list(
        call = match.call(),
        data = data,
        sorted_unique_data = sorted_data,
        data_weights = data_weight / length(data),
        domain = domain,
        opt_theta = result$getValue(theta_vec),
        matrix_A = matrix_A,
        matrix_B = matrix_B,
        penalty_param = penalty_param
    )

    return(result)

}


