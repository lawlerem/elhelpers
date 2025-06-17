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