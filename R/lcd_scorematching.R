#' Computation of Second Derivative of Score Matching Log-concave Density Estimate at Data
#'
#' Computes the second derivative of the (penalized) log-concave score matching
#' density estimate based on i.i.d samples, assuming it is a piecewise linear function
#' with knots at the samples.
#' The output is an object of class "\code{LogConcaveDESM}".
#'
#' @param data A numeric vector whose log-concave density function is to be estimated;
#' missing values are automatically removed.
#' @param domain A numeric vector of length 2 specifying the left and right
#' endpoints of the bounded domain;
#' its components cannot be \code{NA}, \code{NULL}, \code{-Inf}, or \code{Inf}.
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
#'     \item{data_weights}{numeric: weight vector of the sorted unique data.}
#'     \item{domain}{numeric: the bounded domain over the density function is assumed and estimated.}
#'     \item{opt_theta}{numeric: the optimal second derivative of the logarithm of the (penalized)
#'     log-concave score matching density estimate at the sorted unique data and two boundary points.}
#'     \item{matrix_A}{numeric: matrix A used to compute the log-concave score matching density estimate.}
#'     \item{matrix_B}{numeric: matrix B used to compute the log-concave score matching density estimate.}
#'     \item{penalty_param}{numeric: the penalty parameter used to compute the log-concave
#'     score matching density estimate.}
#' }
#'
#' @export
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
lcd_scorematching <- function(data, domain, penalty_param = 1e-1, verbose = FALSE) {

    # preprocess data
    data <- as.numeric(data)
    data <- data[!is.na(data)]

    # some checking
    stopifnot(penalty_param >= 0)
    data <- as.numeric(data)
    # check the length of domain is 2
    stopifnot(length(domain) == 2)
    if ((Inf %in% domain) || (-Inf %in% domain)) {
        stop("domain cannot contain '-Inf' or 'Inf'.")
    }
    if (any(is.na(domain))) {
        stop("domain cannot contain 'NA'.")
    }
    if (any(is.null(domain))) {
        stop("domain cannot contain 'NULL'.")
    }
    domain <- sort(domain)

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
    stopifnot(lower_bound <= min(data), upper_bound >= max(data))

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
    pen_matrix <- matrix(0, N + 2, N + 2)

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
    theta_vec <- CVXR::Variable(N + 2)
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

    class(result) <- 'LogConcaveDESM'

    return(result)
}
