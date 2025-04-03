#' Convert parameters <-> probability vector using continuation ratio logits
#'
#' When using crl, the output vector will have one additional element added to 
#'     the end compared to the input vector to ensure probabilities add to 1.
#' When using icrl, the last element of the input is discarded.
#'
#' @param qpi 
#'     An unbounded parameter vector to be transformed into a probability
#'     vector.
#' @param p 
#'     A probability vector to be transformed to an unbounded probability
#'     vector.
#' 
#' @return A numeric vector
#' 
#' @export
crl<- function(qpi) {
    if( requireNamespace("RTMB", quietly = TRUE) ) {
        plogis<- RTMB::plogis
    } else {
        plogis<- stats::plogis
    }
    pi<- qpi |> plogis() |> c(1)
    qi<- 1 - pi
    cqi<- qi |> head(-1) |> (\(x) c(1, x))() |> cumprod()
    p<- pi * cqi
    return( p )
}
#' @rdname crl
#' @export
icrl<- function(p) {
    pi<- 0 * p |> utils::head(-1)
    for( i in pi |> seq_along() ) {
        pi[i] <- p[i]
        for( j in seq_len(i - 1) ) {
            pi[i]<- pi[i] / (1 - pi[j])
        }
    }

    if( requireNamespace("RTMB", quietly = TRUE) ) {
        qlogis<- RTMB::qlogis
    } else {
        qlogis<- stats::qlogis
    }
    qpi<- pi |> stats::qlogis()
    return(qpi)
}