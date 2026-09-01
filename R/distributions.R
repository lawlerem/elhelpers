#' Dirichlet distribution
#' 
#' @param x 
#'     A probability vector
#' @param pi 
#'     The expected probability vector
#' @param concentration 
#'     The expected closeness of x to y
#' @param log 
#'     Should the log density be returned?
#' 
#' @return The (log) density of x
#' 
#' @export
ddirichlet<- function(x, pi, concentration, log = TRUE) {
    if( !requireNamespace("RTMB", quietly = TRUE) ) {
        stop("Must have RTMB installed.")
    }
    alpha<- concentration * pi
    alg<- 0 * alpha
    for( i in alg |> seq_along() ) alg[i]<- lgamma(alpha[i])
    logconst<- alg |>
        sum() |>
        (\(logdenom) lgamma(concentration) - logdenom)() 
    ll<- ((alpha - 1) * log(x)) |>
        sum() |>
        (\(y) y + logconst)()
    if( !log ) ll<- ll |> exp()
    return(ll)
}

#' Dirichlet distribution
#' 
#' @param n 
#'     The number of Dirichlet samples
#' @param pi 
#'     The expected probability vector
#' @param concentration 
#'     The expected closeness of x to y
#' 
#' @return A n * length(pi) matrix
#' 
#' @export
rdirichlet<- function(n, pi, concentration) {
    alpha<- concentration * pi
    Y<- alpha |> sapply(\(a) rgamma(n, shape = a, scale = 1))
    dim(Y)<- c(n, length(pi))
    Y<- sweep(Y, 1, rowSums(Y), `/`)
    return(Y)
}


#' Gamma distribution parameterized with mean and variance
#'
#' @param x
#'     A positive numeric variable
#' @param mean
#'     Mean of the Gamma distribution
#' @param var
#'     Variance of the Gamma distribution
#' @param log
#'     Should the log density be returned?
#'
#' @export
dgammamv<- function(x, mean, var, log = TRUE) {
    # shape = a; scale = s
    # mu = a * s; var = a * s^2
    # mu = a * s; var = mu * s
    # s = mu / a; var = mu * mu / a
    # s = mu / a; a = mu^2 / var
    # s = mu / (mu^2 / var) = var / mu; a = mu^2 / var
    stats::dgamma(x, shape = mean^2 / var, scale = var / mean, log = log)
}

#' Gamma distribution parameterized with mean and variance
#'
#' @param n
#'     The number of Gamma samples
#' @param mean
#'     Mean of the Gamma distribution
#' @param var
#'     Variance of the Gamma distribution
#'
#' @export
rgammamv<- function(n, mean, var) {
    # shape = a; scale = s
    # mu = a * s; var = a * s^2
    # mu = a * s; var = mu * s
    # s = mu / a; var = mu * mu / a
    # s = mu / a; a = mu^2 / var
    # s = mu / (mu^2 / var) = var / mu; a = mu^2 / var
    stats::rgamma(n, shape = mean^2 / var, scale = var / mean)
}