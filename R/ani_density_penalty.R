#' Building a gif to view the effects of the penalty parameter on the density estimate
#'
#' Creates a gif to visualize how different penalty parameters can affect the density estimates.
#' The resulting gif shows the density estimates with increasing values of the penalty parameters.
#' For the small penalty parameter, the resulting density estimate should be concentrating at
#' the place where the data are abundant. As the penalty parameter increases, the resulting density
#' estimate is getting closer to the density function of the uniform distribution over the domain.
#'
#' @param data A numeric vector whose log-concave density function is to be estimated;
#' missing values are automatically removed.
#' @param domain A numeric vector of length 2 specifying the left and right
#' endpoints of the bounded domain;
#' its components cannot be \code{NA}, \code{NULL}, \code{-Inf}, or \code{Inf}.
#' @param penalty_params_seq Penalty parameter for computing the density estimate; must be non-negative.
#' Default is \code{1e-1}.
#' @param plot_domain A numeric vector to indicate the domain of the plot.
#' @param plot_points_cnt A numeric to indicate the number of points for evaluating and plotting.
#' Default is \code{100}.
#'
#' @import ggplot2
#' @import gganimate
#' @import tidyr
#' @import forcats
#' @import transformr
#'
#' @return A list of plots with multiple layers for creating the gif.
#' To view the resulting gif, please use
#' \code{gganimate::animate(anim, renderer = gganimate::gifski_renderer())}.
#'
#' @export
#'
#' @examples
#' library(transformr)
#' data <- rnorm(200)
#' domain <- c(-5, 5)
#' penalty_params_seq <- c(0, exp(seq(-10, 1, length.out = 20)))
#' anim <- ani_density_penalty(
#' data = data,
#' domain = domain,
#' penalty_params_seq = penalty_params_seq,
#' plot_domain = domain,
#' plot_points_cnt = 500)
#'
#' gganimate::animate(anim, renderer = gganimate::gifski_renderer())
#'
#'
ani_density_penalty <- function(data, domain, penalty_params_seq, plot_domain, plot_points_cnt = 100) {

    if(any(penalty_params_seq < 0)) {
        stop("All values in penalty_params_seq must be non-negative.")
    }

    penalty_params_seq <- sort(penalty_params_seq)

    # prepare the newx to evaluate and plot
    newx <- seq(plot_domain[1], plot_domain[2], length.out = plot_points_cnt)
    all_denvals_list <- list()

    for (i in 1:length(penalty_params_seq)) {

        lambda_val <- penalty_params_seq[i]
        message(paste0("Penalty parameter value: ", lambda_val, "."))

        result <- lcd_scorematching(
            data = data,
            domain = domain,
            penalty_param = lambda_val,
            verbose = FALSE)

        den_vals <- evaluate_density(
            scorematching_logconcave = result,
            newx = newx)

        colnames(den_vals) <- c("sorted_newx", paste0("pen", i))
        all_denvals_list[[paste0("pen", i)]] <- den_vals

    }

    # concatenate all data frames into a large one
    df <- all_denvals_list[[paste0("pen", 1)]]
    for (i in 2:length(penalty_params_seq)) {

        df <- merge(
            df,
            all_denvals_list[[paste0("pen", i)]],
            by = "sorted_newx")

    }

    new_df <- tidyr::gather(
        df,
        key = "key",
        value = "value",
        -sorted_newx)
    new_df$key <- forcats::as_factor(new_df$key)
    old_code <- paste0("pen", 1:length(penalty_params_seq))
    new_code <- forcats::as_factor(as.character(round(penalty_params_seq, 10)))
    new_df$key <- new_code[match(new_df$key, old_code)]

    plot <- ggplot2::ggplot() +
        ggplot2::geom_line(
            data = new_df,
            ggplot2::aes(x = sorted_newx, y = value),
            color = "red",
            size = 0.8) +
        ggplot2::geom_rug(
            data = data.frame(original_data = data),
            ggplot2::aes(x = original_data),
            color = 'black',
            alpha = 0.5
        ) +
        ggplot2::labs(
            x = 'x',
            y = 'density') +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom")

    # if (plot_hist) {
    #
    #     # use the Freedman-Diaconis rule for the binwidth
    #     if (is.null(hist_binwidth)) {
    #         NN <- length(data)
    #         hist_binwidth <- 2 * stats::IQR(data) / (NN ** (1/3))
    #     }
    #
    #     plot <- plot + ggplot2::geom_histogram(
    #         data = data.frame(x = data),
    #         ggplot2::aes(x = x, y = ..density..),
    #         binwidth = hist_binwidth,
    #         color = "darkblue",
    #         fill = "lightblue",
    #         alpha = hist_alpha)
    #
    # }

    plot <- plot +
        gganimate::transition_states(
            key,
            transition_length = 1,
            state_length = 2) +
        gganimate::ease_aes('cubic-in-out') +
        ggplot2::ggtitle("Density Estimates with Increasing Penalty Parameters"
                         ,subtitle = "Penalty Parameter = {closest_state}")

    return(plot)

}
