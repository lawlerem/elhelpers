#' Make a tessellation of a polygon.
#' 
#' @param shape 
#'     A polygon of class sf. If multiple polygons are given, they are combined
#'     into a single polygon before creating the tesselation.
#' @param n 
#'     The number of polygons in the tesselation
#' @param type 
#'     The type of tesselation. See sf::st_sample.
#' 
#' @return An sf data.frame containing the tesselation polygons and an index
#'     column.
#' 
#' @export
tesselate<- function(shape, n, type = "hexagonal") {
    shape<- shape |> sf::st_union()
    tesselation<- shape |>
        sf::st_sample(
            n,
            type = type
        ) |>
        sf::st_combine() |>
        sf::st_voronoi() |>
        sf::st_cast() |>
        sf::st_as_sf() |>
        sf::st_intersection(shape)
    tesselation$tesselation<- tesselation |> nrow() |> seq()
    return(tesselation)
}

#' Create a directed adjacency matrix from a tesselation
#' 
#' The tesselation is converted into a directed acyclic graph
#' 
#' @param tesselation
#'     An sf data.frame containing polygons
#' @param sf_predicate
#'     The function used to determine node adjacency, see ?sf::geos_binary_pred.
#' @param ...
#'     Additional arguments to pass to sf_predicate.
#' 
#' @return A direct acyclic graph from the igraph package
#' 
#' @export
dagify<- function(tesselation, sf_predicate = sf::st_touches, ...) {
    graph<- tesselation |>
        sf_predicate(sparse = FALSE, ...) |>
        (\(x) {x[x |> lower.tri(diag = TRUE)]<- FALSE; return(x)})() |>
        igraph::graph_from_adjacency_matrix()
    return(graph)
}

#' Create a key to translate indices between rasters and tesselations
#' 
#' @param target_raster A stars raster
#' @param prediction_raster (Optional) A stars raster
#' @param tesselation (Optional) An sf data.frame containing polygons
#' 
#' @return A stars raster with the same geometry as the target raster with one 
#'     or more of the following layers:
#' - target_cell: The cell index for the target raster
#' - prediction_cell: The index for the cell of the prediction raster that 
#'     coincides with the target cell
#' - tesselation_cell: The index for the polygon of the tesselation that
#'     coincides with the target cell
#' 
#' @export
make_cell_key<- function(
        target_raster,
        prediction_raster,
        tesselation
    ) {
    cell_key<- target_raster |>
        _[1] |>
        (\(x) {x[[1]][]<- seq_along(x[[1]]); return(x)})() |>
        stats::setNames("target_cell")

    if( !missing(prediction_raster) ) {
        prediction_raster<- prediction_raster |>
            _[1] |>
            (\(x) {x[[1]][]<- seq_along(x[[1]]); return(x)})() |>
            stats::setNames("prediction_cell") |>
            stars::st_warp(target_raster)
        cell_key<- c(cell_key, prediction_raster)
    }

    if( !missing(tesselation) ) {
        tesselation_raster<- tesselation |>
            (\(x) 
                sf::st_sf(
                    tesselation_cell = x |> nrow() |> seq(),
                    geometry = x |> sf::st_geometry()
                )
            )() |>
            stars::st_rasterize(cell_key, align = TRUE) |>
            stars::st_warp(cell_key)
        cell_key<- c(cell_key, tesselation_raster)
    }
    return(cell_key)
}

