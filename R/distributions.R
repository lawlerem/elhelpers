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
    ll<- alpha |>
        sapply(lgamma) |>
        sum() |>
        (\(x) lgamma(concentration) - x)() |>
        (\(x) x + (alpha - 1) * log(x) |> sum())()
    if( !log ) ll<- ll |> exp()
    return(ll)
}