#' An RTMB-compatible version of base::lgamma
#' 
#' On package-load, makes an AD-taped version of lgamma if RTMB is installed.
#' 
#' @param x 
#'     A numeric vector
lgamma<- function(x) NULL

.onLoad<- function(
        libname, 
        pkgname
    ) {
    if( requireNamespace("RTMB", quietly = TRUE) ) {
        lgamma<- base::lgamma |> RTMB::MakeTape(1)
        assign(
            "lgamma",
            lgamma,
            environment() |> parent.env()
        )
    }
    invisible()
}