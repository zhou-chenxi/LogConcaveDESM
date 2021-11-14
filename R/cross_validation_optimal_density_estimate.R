#' Cross Validation to Choose the Best Penalty Parameter
#'
#' Given a numeric vector of penalty parameter candidates, use the cross validation to
#' choose the best one.
#'
#' @param data A numeric vector whose log-concave density function is to be estimated;
#' missing values are automatically removed.
#' @param domain A numeric vector of length 2 specifying the left and right
#' endpoints of the bounded domain; its components cannot be \code{NaN}.
#' @param penalty_param_candidates A numeric vector of the penalty parameter candidates;
#' each element must be non-negative.
#' @param fold_number An integer to indicate the number of folds for cross validation.
#' Default is \code{5}.
#'
#' @return A object with class "LogConcaveDESM" with the penalty parameter being the optimal choice.
#' @export
#'
#' @seealso \code{\link{lcd_scorematching}}
#'
#' @examples
#' set.seed(1119)
#' N <- 100
#' data <- rnorm(N)
#' domain <- c(-5, 5)
#' lambda_cand <- exp(c(-Inf, seq(-3, 1, by = 0.5)))
#'
#' opt_result <- cv_optimal_density_estimate(data, domain, lambda_cand)
#'
cv_optimal_density_estimate <- function(data, domain, penalty_param_candidates, fold_number = 5) {

    data <- as.numeric(data)
    if (any(penalty_param_candidates < 0)) {
        stop("penalty_param_candidates should all be positive but contains negative elements.")
    }

    stopifnot(length(domain) == 2)
    if (any(is.nan(domain))) {
        stop("domain cannot contain 'NaN'.")
    }
    domain <- sort(domain)

    # split data for cross validation
    split_id <- sample(1:fold_number, length(data), replace = TRUE)
    loss_record <- rep(NA, fold_number)

    for (i in 1:length(penalty_param_candidates)) {

        lambda_val <- penalty_param_candidates[i]
        message(paste0("Penalty parameter value: ", lambda_val))

        loss_value <- 0

        for (j in 1:fold_number) {

            train_data <- data[split_id != j]
            test_data <- data[split_id == j]

            estimator <- lcd_scorematching(
                data = train_data,
                domain = domain,
                penalty_param = lambda_val)

            loss_value <- (loss_value +
                               evaluate_scorematching_loss(
                                   scorematching_logconcave = estimator,
                                   new_data = test_data)
            )
        }

        loss_record[i] <- loss_value

    }

    ddff <- data.frame(
        penalty_param = penalty_param_candidates,
        loss_vals = loss_record)

    message(ddff)

    opt_penalty_param <- penalty_param_candidates[which.min(loss_record)]
    message(paste0("Optimal penality parameter is ", opt_penalty_param, "."))

    opt_estimator <- lcd_scorematching(
        data = data,
        domain = domain,
        penalty_param = opt_penalty_param)

    return(opt_estimator)

}
