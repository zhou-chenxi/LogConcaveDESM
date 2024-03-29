% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ani_density_penalty.R
\name{ani_density_penalty}
\alias{ani_density_penalty}
\title{Building an animation to view the effects of the penalty parameter on density estimates}
\usage{
ani_density_penalty(
  data,
  domain,
  penalty_params_seq,
  plot_domain,
  plot_points_cnt = 100
)
}
\arguments{
\item{data}{A numeric vector whose log-concave density function is to be estimated;
missing values are automatically removed.}

\item{domain}{A numeric vector of length 2 specifying the left and right
endpoints of the bounded domain;
its components cannot be \code{NA}, \code{NULL}, or \code{NaN}.}

\item{penalty_params_seq}{Penalty parameter for computing the density estimate; must be non-negative.}

\item{plot_domain}{A numeric vector to indicate the domain of the plot.}

\item{plot_points_cnt}{A numeric to indicate the number of points for evaluating and plotting.
Default is \code{100}.}
}
\value{
A list of plots with multiple layers for creating the animation.
To view the resulting animation, please use
\code{gganimate::animate(anim, renderer = gganimate::gifski_renderer())}.
}
\description{
Creates an animation to visualize how different penalty parameters can affect the density estimates.
The resulting animation shows the density estimates with increasing values of the penalty parameters.
}
\examples{
library(transformr)
data <- rnorm(200)
domain <- c(-5, 5)
penalty_params_seq <- c(0, exp(seq(-10, 1, length.out = 20)))
anim <- ani_density_penalty(
data = data,
domain = domain,
penalty_params_seq = penalty_params_seq,
plot_domain = domain,
plot_points_cnt = 500)

gganimate::animate(anim, renderer = gganimate::gifski_renderer())


}
