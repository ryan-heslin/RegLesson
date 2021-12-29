
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RegLesson

<!-- badges: start -->

<!-- badges: end -->

RegLesson originated from my preparations the summer before taking a
linear regression course. I decided to implement linear regression from
the ground up, using only base R. Instead of relying on `lm`, I built a
`LinRegProblem` class to store a linear model and associated data. I
extended this with functions designed to make it easier to solve
textbook problems. I have no illusions that my work is better than R’s
native linear modeling functions, but implementing it was a useful
exercise.

You can create `LinRegProblem` objects like so:

``` r
library(RegLesson)

X <- LinRegProblem(mtcars[, -1], mtcars$mpg)
```

This represents a regression of `mpg` on all other `mtcars` variables.

the `do_problems` function evaluates expressions in an environment and
prints the results. Required, omitted arguments are replaced with their
values in the environment. Now we can check whether the model passes the
overall \(F\) test and get Working-Hoteling confidence bands for the
observations.

``` r
do_problems(
  X,
  regression_relation(),
  Working_Hoteling(X_h = cbind(1, as.matrix(mtcars[, -1])))
)
#> F test of regression relation
#> alpha = 0.05
#> Test statistic: 13.9324635520384
#> Threshold: 0.361846621741306
#> Decision: Reject null
#>  13.93246 
#>  
#>            [,1]     [,2]     [,3]
#>  [1,] 19.621795 22.59951 25.57722
#>  [2,] 19.195269 22.11189 25.02850
#>  [3,] 23.604901 26.25064 28.89639
#>  [4,] 18.653752 21.23740 23.82106
#>  [5,] 15.275194 17.69343 20.11167
#>  [6,] 17.506579 20.38304 23.25950
#>  [7,] 11.295468 14.38626 17.47705
#>  [8,] 19.384840 22.49601 25.60719
#>  [9,] 19.754631 24.41909 29.08354
#> [10,] 15.151641 18.69903 22.24642
#> [11,] 15.876924 19.19165 22.50638
#> [12,] 11.190642 14.17216 17.15368
#> [13,] 13.226585 15.59957 17.97256
#> [14,] 13.181823 15.74222 18.30263
#> [15,]  8.720951 12.03401 15.34708
#> [16,]  7.926723 10.93644 13.94615
#> [17,]  7.495369 10.49363 13.49189
#> [18,] 25.482659 27.77291 30.06315
#> [19,] 26.023080 29.89674 33.77040
#> [20,] 26.899764 29.51237 32.12497
#> [21,] 20.078869 23.64310 27.20734
#> [22,] 14.415191 16.94305 19.47092
#> [23,] 15.470953 17.73218 19.99341
#> [24,]  9.847552 13.30602 16.76449
#> [25,] 14.238578 16.69168 19.14478
#> [26,] 26.252151 28.29347 30.33479
#> [27,] 21.878922 26.15295 30.42699
#> [28,] 24.081571 27.63627 31.19098
#> [29,] 14.460895 18.87004 23.27918
#> [30,] 16.308397 19.69383 23.07926
#> [31,]  9.600629 13.94112 18.28161
#> [32,] 21.450210 24.36827 27.28633
```

This package contains functions for several diagnostic tests as well.
Have fun\!

## Installation

You can install the development version of RegLesson from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ryan-heslin/RegLesson")
```
