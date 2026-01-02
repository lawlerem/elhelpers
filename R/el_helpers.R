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
        .orig_TapeConfig<- RTMB::TapeConfig()
        RTMB::TapeConfig(comparison = "tape")
        sp<- RTMB::MakeTape(
            \(x) (x < 1) * exp(x - 1) + (x >= 1) * x,
            0
        )
        assign("sp", sp, .GlobalEnv)
        softplus<- \(x) {
            y<- x |> RTMB::sapply(sp)
            if( !is.null(dim(x)) ) y<- array(y, dim = dim(x))
            y
        }
        isp<- RTMB::MakeTape(
            \(y) (y < 1) * (log(y) + 1) + (y >= 1) * y,
            1
        )
        isoftplus<- \(x) {
            y<- x |> RTMB::sapply(isp)
            if( !is.null(dim(x)) ) y<- array(y, dim = dim(x))
            y
        }
        assign(
            "softplus",
            softplus,
            environment() |> parent.env()
        )
        assign(
            "isoftplus",
            isoftplus,
            environment() |> parent.env()
        )
        RTMB::TapeConfig |> do.call(as.list(.orig_TapeConfig))
    }
    invisible()
}