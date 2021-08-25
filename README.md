
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Log-Concave Density Estimation Using the Score Matching Loss Function

`LogConcaveDESM` is an package to compute and visualize the log-concave
score matching derocnsity estimate over a bounded domain. It also
contains functions to plot the first and second derivatives of the
log-density estimate. Furthermore, functions to compute various
distances between two probability distributions are provided to assess
the quality of the density estimate.

## Installation

The development version of the `LogConcaveDESM` can be installed from
the GitHub using `devtools`:

``` r
# install.packages("devtools")
devtools::install_github("zhoucx1119/LogConcaveDESM")
```

## Load the Package

``` r
library(LogConcaveDESM, warn.conflicts = FALSE)
library(patchwork)
library(ggplot2)
```

## Computing and Visualizing Unpenalized Log-concave Score Matching Density Estimate

We show how to use functions in `LogConcaveDESM` to compute and
visualize the log-concave score matching density estimate over a bounded
domain.

We first simulate data from a truncated normal distribution over the
interval \[ − 5, 5\].

``` r
set.seed(2021)
N <- 200
parent_mean <- 0
parent_sd <- 1
domain <- c(-5, 5)
t_normal <- truncated_normal(parent_mean, parent_sd, domain)
data <- t_normal$sampling(N)
```

We plot the histogram of the data below.

``` r
binwidth <- 2 * stats::IQR(data) / N ** (1/3)
ggplot() + geom_histogram(
  data = data.frame(x = data), 
  aes(x = x, y = ..density..), 
  binwidth = binwidth,
            color = "darkblue",
            fill = "lightblue") + 
  coord_cartesian(xlim = domain, ylim = c(-0.01, 0.5)) + 
  theme_bw()
```

![](README-hist_data-1.png)<!-- -->

We assume that the second derivative of the logarithm of the log-concave
score matching density estimate is a piecewise linear function with
knots at the data points. With no penalty, the resulting optimal second
derivative values at data points can be computed as below:

``` r
result <- lcd_scorematching(
  data = data, 
  domain = domain, 
  penalty_param = 0, 
  verbose = TRUE
)
#> -----------------------------------------------------------------
#>            OSQP v0.6.0  -  Operator Splitting QP Solver
#>               (c) Bartolomeo Stellato,  Goran Banjac
#>         University of Oxford  -  Stanford University 2019
#> -----------------------------------------------------------------
#> problem:  variables n = 202, constraints m = 203
#>           nnz(P) + nnz(A) = 20504
#> settings: linear system solver = qdldl,
#>           eps_abs = 1.0e-05, eps_rel = 1.0e-05,
#>           eps_prim_inf = 1.0e-04, eps_dual_inf = 1.0e-04,
#>           rho = 1.00e-01 (adaptive),
#>           sigma = 1.00e-06, alpha = 1.60, max_iter = 10000
#>           check_termination: on (interval 25),
#>           scaling: on, scaled_termination: off
#>           warm start: on, polish: on, time_limit: off
#> 
#> iter   objective    pri res    dua res    rho        time
#>    1  -4.3429e+00   1.35e+01   4.63e-03   1.00e-01   5.74e-03s
#>  200  -6.7500e+01   2.79e+00   1.88e-03   1.98e-02   2.28e-02s
#>  400  -7.1064e+01   4.32e+00   4.59e-04   3.81e-03   3.54e-02s
#>  600  -6.5517e+01   1.60e+00   1.22e-04   3.81e-03   4.64e-02s
#>  800  -6.7135e+01   4.42e-01   3.10e-05   3.81e-03   5.73e-02s
#> 1000  -6.6824e+01   1.01e-01   6.51e-06   3.81e-03   6.86e-02s
#> plsh  -6.6874e+01   4.27e-10   3.14e-09   --------   7.21e-02s
#> 
#> status:               solved
#> solution polish:      successful
#> number of iterations: 1000
#> optimal objective:    -66.8737
#> run time:             7.21e-02s
#> optimal rho estimate: 1.93e-03
#> The status of solving the constrained quadratic optimization problem is: optimal.
```

Then, the second derivative of the logarithm of the log-concave score
matching density estimate at arbitrary points can be computed using the
`evaluate_logden_deriv2` function, and can be visualized using the
`plot_logden_deriv2` function.

``` r
logden_deriv2_vals <- evaluate_logden_deriv2(
  scorematching_logconcave = result, 
  newx = seq(domain[1], domain[2], length.out = 11)
)
logden_deriv2_vals
#>    newx_sorted logderiv2_vals
#> 1           -5  -1.488096e-23
#> 2           -4  -7.214839e+01
#> 3           -3  -1.442968e+02
#> 4           -2   2.368584e-13
#> 5           -1   2.606746e-12
#> 6            0  -1.251196e-10
#> 7            1   1.092312e-10
#> 8            2  -5.711015e-11
#> 9            3  -7.400681e-12
#> 10           4  -3.700341e-12
#> 11           5   0.000000e+00

plot_deriv2 <- plot_logden_deriv2(
  scorematching_logconcave = result, 
  plot_domain = domain, 
  plot_points_cnt = 500)
plot_deriv2
```

![](README-deriv2-1.png)<!-- -->

Subsequently, the first derivative of the logarithm of the log-concave
score matching density estimate at arbitrary points can be computed
using the `evaluate_logden_deriv1` function, and can be visualized using
the `plot_logden_deriv1` function.

``` r
logden_deriv1_vals <- evaluate_logden_deriv1(
  scorematching_logconcave = result, 
  newx = seq(domain[1], domain[2], length.out = 11)
)
logden_deriv1_vals
#>    newx_sorted logderiv1_vals
#> 1           -5     294.125904
#> 2           -4     258.051708
#> 3           -3     149.829119
#> 4           -2      13.667716
#> 5           -1       9.167263
#> 6            0       5.598630
#> 7            1     -13.979790
#> 8            2     -33.035110
#> 9            3     -33.035110
#> 10           4     -33.035110
#> 11           5     -33.035110

plot_deriv1 <- plot_logden_deriv1(
  scorematching_logconcave = result, 
  plot_domain = domain, 
  plot_points_cnt = 500)
plot_deriv1
```

![](README-deriv1-1.png)<!-- -->

Then, the (un-normalized) log-density estimate at arbitrary point can be
computed using the `evaluate_logdensity` function, and can be visualized
using the `plot_logdensity` function.

``` r
logden_vals <- evaluate_logdensity(
  scorematching_logconcave = result, 
  newx = seq(domain[1], domain[2], length.out = 11)
)
logden_vals
#>    newx_sorted logdensity_vals
#> 1           -5          0.0000
#> 2           -4        282.1012
#> 3           -3        492.0540
#> 4           -2        562.3994
#> 5           -1        574.0298
#> 6            0        580.8876
#> 7            1        575.0199
#> 8            2        553.7940
#> 9            3        520.7589
#> 10           4        487.7238
#> 11           5        454.6887

plot_logden <- plot_logdensity(
  scorematching_logconcave = result, 
  plot_domain = domain, 
  plot_points_cnt = 500)
plot_logden
```

![](README-logden-1.png)<!-- -->

Finally, to evaluate and visualize the density estimate itself, we use
`evaluate_density` and `plot_density` functions, respectively.

``` r
den_vals <- evaluate_density(
  scorematching_logconcave = result, 
  newx = seq(domain[1], domain[2], length.out = 11)
)
den_vals
#>    newx_sorted  density_vals
#> 1           -5 5.718881e-253
#> 2           -4 1.871947e-130
#> 3           -3  2.842021e-39
#> 4           -2  1.009840e-08
#> 5           -1  1.135730e-03
#> 6            0  1.080366e+00
#> 7            1  3.056907e-03
#> 8            2  1.849261e-12
#> 9            3  8.318257e-27
#> 10           4  3.741679e-41
#> 11           5  1.683064e-55

plot_den <- plot_density(
  scorematching_logconcave = result,
  plot_domain = domain,
  plot_points_cnt = 500, 
  plot_hist = TRUE
  )
plot_den
```

![](README-den-1.png)<!-- -->

## Computing and Visualizing Penalized Log-concave Score Matching Density Estimate

As we can see above, the un-penalized log-concave score matching density
estimate is too concentrated at the place where data are abundant. To
remedy this, we consider penalized log-concave score matching density
estimate, where the optimal penalty parameter is chosen using the
`cv_optimal_density_estimate` function.

``` r
lambda_cand <- exp(seq(-5, 1, by = 0.5))
opt_den <- cv_optimal_density_estimate(
  data = data, 
  domain = domain, 
  penalty_param_candidates = lambda_cand, 
  fold_number = 5
)
#> Penalty parameter value: 0.00673794699908547
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal.
#> Penalty parameter value: 0.0111089965382423
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal.
#> Penalty parameter value: 0.0183156388887342
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal.
#> Penalty parameter value: 0.0301973834223185
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal.
#> Penalty parameter value: 0.0497870683678639
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal.
#> Penalty parameter value: 0.0820849986238988
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal.
#> Penalty parameter value: 0.135335283236613
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal.
#> Penalty parameter value: 0.22313016014843
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal.
#> Penalty parameter value: 0.367879441171442
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal.
#> Penalty parameter value: 0.606530659712633
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal.
#> Penalty parameter value: 1
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal.
#> Penalty parameter value: 1.64872127070013
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal.
#> Penalty parameter value: 2.71828182845905
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal. 
#> The status of solving the constrained quadratic optimization problem is: optimal.
#>    penalty_param  loss_vals
#> 1    0.006737947 10.4393950
#> 2    0.011108997  6.5639880
#> 3    0.018315639  3.6106832
#> 4    0.030197383  1.3993944
#> 5    0.049787068  0.1192177
#> 6    0.082084999 -0.6599348
#> 7    0.135335283 -0.9945474
#> 8    0.223130160 -1.0987366
#> 9    0.367879441 -1.0469678
#> 10   0.606530660 -0.9045135
#> 11   1.000000000 -0.7083408
#> 12   1.648721271 -0.5121017
#> 13   2.718281828 -0.3486011
#> [1] "Optimal penality parameter is 0.22313016014843."
#> The status of solving the constrained quadratic optimization problem is: optimal.
```

The resulting density estimate is the following.

``` r
plot_opt_den <- plot_density(
  scorematching_logconcave = opt_den,
  plot_domain = domain,
  plot_points_cnt = 500, 
  plot_hist = TRUE
  )
plot_opt_den
```

![](README-plot-opt-1.png)<!-- -->

We could assess the quality of this density estimate using various
metrics between two probability distributions, for example, the
Kullback-Leibler divergence, the Hyvarinen divergence, the *L*1
distance, and the Hellinger distance.

``` r
kl <- kl_div(
  true_density = t_normal, 
  density_estimate = opt_den) 

hyva <- hyvarinen_div(
  true_density = t_normal, 
  density_estimate = opt_den) 

l1 <- L1_dist(
  true_density = t_normal, 
  density_estimate = opt_den) 

he <- hellinger_dist(
  true_density = t_normal, 
  density_estimate = opt_den) 

metric_table <- cbind(
  Metrics = c("Kullback-Leibler Divergence", "Hyvarinen divergence", 
              "L1 Distance", "Hellinger Distance"), 
  Values = round(c(kl, hyva, l1, he), 5)
)

knitr::kable(metric_table)
```

| Metrics                     | Values  |
|:----------------------------|:--------|
| Kullback-Leibler Divergence | 0.10629 |
| Hyvarinen divergence        | 0.10367 |
| L1 Distance                 | 0.30989 |
| Hellinger Distance          | 0.1748  |