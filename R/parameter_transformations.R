#' Convert parameters <-> probability vector using continuation ratio logits
#'
#' When using crl, the output vector will have one additional element added to 
#'     the end compared to the input vector to ensure probabilities add to 1.
#' When using icrl, the last element of the input is discarded.
#' Can be vectorized, in which case elements with the same (n-1) indices
#'     will be treated as parameters for the same probability vector.
#' For example if x is a 3 dimensional array, then the operation is applied to
#'     each vector x[i, j, ].
#'
#' @param qpi 
#'     An unbounded parameter vector to be transformed into a probability
#'         vector.
#'     
#' @param p 
#'     A probability vector to be transformed to an unbounded probability
#'     vector.
#' 
#' @return A numeric vector
#' 
#' @export
crl<- function(qpi) {
    do_ad<- (
            requireNamespace("RTMB", quietly = TRUE) && 
            "advector" %in% class(qpi)
        )
    if( qpi |> dim() |> is.null() ) dim(qpi)<- c(1, qpi |> length())
    if( do_ad ) {
        plogis<- RTMB::plogis
    } else {
        plogis<- stats::plogis
    }
    pi<- qpi |> plogis()
    qi<- 1 - pi

    p_dim<- qpi |> dim()
    p_dim[length(p_dim)]<- p_dim[length(p_dim)] + 1
    p<- 0 |> array(dim = p_dim)

    cqi<- 1 |> array(dim = p_dim |> head(-1))
    if( do_ad ) {
        p<- p |> RTMB::AD()
        cqi<- cqi |> RTMB::AD()
    } else {}
    selector<- p_dim |>
        head(-1) |>
        lapply(seq) |>
        expand.grid() |>
        as.matrix()

    for( i in p |> dim() |> tail(1) |> seq() ) {
        if( i > 1 ) {
            cqi[selector]<- cqi[selector] *  qi[selector |> cbind(i - 1)]
        }
        if( i == (p |> dim() |> tail(1)) ) {
            p[selector |> cbind(i)]<- cqi[selector]
            next
        }
        p[selector |> cbind(i)]<- pi[selector |> cbind(i)] * cqi[selector]
    }
    return( p )
}
#' @rdname crl
#' @export
icrl<- function(p) {
    do_ad<- (
            requireNamespace("RTMB", quietly = TRUE) && 
            "advector" %in% class(p)
        )
    if( p |> dim() |> is.null() ) dim(p)<- c(1, p |> length())
    if( do_ad ) {
        qlogis<- RTMB::qlogis
    } else {
        qlogis<- stats::qlogis
    }
    qpi_dim<- p |> dim()
    qpi_dim[length(qpi_dim)]<- qpi_dim[length(qpi_dim)] - 1
    qpi<- 0 |> array(dim = qpi_dim)

    denom<- 1 |> array(dim = qpi_dim |> head(-1))
    if( do_ad ) {
        qpi<- qpi |> RTMB::AD()
        denom<- cnp |> RTMB::AD()
    } else {}
    selector<- qpi_dim |>
        head(-1) |>
        lapply(seq) |>
        expand.grid() |>
        as.matrix()

    for( i in p |> dim() |> tail(1) |> (\(x) x - 1)() |> seq() ) {
        if( i > 1 ) {
            denom[selector]<- denom[selector] - p[selector |> cbind(i - 1)]
        }
        qpi[selector |> cbind(i)]<- (
                p[selector |> cbind(i)] / denom[selector]
            ) |>
            qlogis()
    }
    return(qpi)
}