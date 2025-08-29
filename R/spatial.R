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
        sf::st_intersection(shape) |>
        sf::st_cast() |>
        sf::st_cast("POLYGON")
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
#' @param adjacency_power
#'     The (undirected) adjacency matrix will be raised to this power to
#'     increase the neighbour size. E.g. if p1 -> p2 -> ... -> pk is a chain of
#'     parents according to sf_predicate, then p1 -> pk if the adjacency_power
#'     is k.
#' @param ...
#'     Additional arguments to pass to sf_predicate.
#' 
#' @return A direct acyclic graph from the igraph package
#' 
#' @export
dagify<- function(
        tesselation, 
        sf_predicate = sf::st_touches,
        adjacency_power = 1,
        ...
    ) {
    graph<- tesselation |>
        sf_predicate(sparse = FALSE, ...) |>
        (\(x) {
            adjacency_power<- adjacency_power |> round()
            if( adjacency_power == 1 ) return(x)
            diag(x)<- 1
            x<- x |> Matrix::Matrix(sparse = TRUE)
            y<- x
            i<- 2
            while( i <= adjacency_power ) {
                y<- y %*% x
                i<- i + 1
            }
            x<- y |> as.matrix() |> (\(x) (x != 0))()
            return(x)
        })() |>
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



#' Transform one polygon geometry to another using area-weighted averages
#' 
#' @param x An 'sf' or 'stars' object with POLYGON geometries
#' @param y An 'sf' object with the target POLYGON geometries
#' 
#' @return An 'sf' or 'stars' object with the same attributes and dimensions as
#'     x and the same geometry as y.
#' 
#' @export
st_repolygonize<- function(
        source,
        target
    ) {
    if( source |> inherits("sf") ) return(st_repolygonize_sf(source, target))
    if( source |> inherits("stars") ) return(st_repolygonize_stars(source, target))
    stop("source must be an 'sf' or 'stars' object.")
}
st_repolygonize_sf<- function(
        source,
        target
    ) {
    source<- source |>
        (\(x) {
            is_complete<- x |>
                sf::st_drop_geometry() |>
                complete.cases()
            x[is_complete, ]
        })()
    sourcegeo<- source |>
        sf::st_geometry() |> 
        (\(x) sf::st_sf(source = x |> seq_along(), geometry = x))()
    targetgeo<- target |>
        sf::st_geometry() |>
        (\(y) sf::st_sf(target = y |> seq_along(), geometry = y))()
    refined<- sourcegeo |> sf::st_intersection(targetgeo)
    refined$area<- refined |> sf::st_area()
    refined$weight<- refined |>
        with(tapply(area, target, \(x) x / sum(x))) |>
        do.call(c, args = _) |>
        units::drop_units()

    A<- Matrix::sparseMatrix(
            x = refined |> _$weight,
            i = refined |> _$target,
            j = refined |> _$source,
            dims = c(targetgeo |> nrow(), sourcegeo |> nrow())
        )
    ans<- sf::st_as_sf(
            source |> 
                sf::st_drop_geometry() |> 
                as.matrix() |> 
                (\(x) A %*% x)() |>
                as.matrix() |>
                as.data.frame(),
            geometry = targetgeo |> sf::st_geometry()
        )
    return(ans)
}
st_repolygonize_stars<- function(
        x,
        y
    ) {
    xgeo<- x |>
        sf::st_geometry() |>
        (\(x) sf::st_sf(x = x |> seq_along(), geometry = x))()
    ygeo<- y |>
        sf::st_geometry() |>
        sf::st_intersection(xgeo |> sf::st_geometry() |> sf::st_union()) |>
        (\(y) sf::st_sf(y = y |> seq_along(), geometry = y))()
    ydim<- x |> stars::st_dimensions()

    geodim<- ydim |>
        seq_along() |>
        sapply(
            function(v) {
                ydim |> 
                    stars::st_get_dimension_values(v) |>
                    inherits("sfc")
            }
        ) |>
        which() |>
        min()
    ydim[[geodim]]<- ygeo |>
        stars::st_as_stars() |>
        stars::st_dimensions() |>
        _[[1]]
    y<- x |>
        names() |>
        lapply(function(v) array(0, dim = ydim |> dim())) |>
        setNames(x |> names()) |>
        stars::st_as_stars(dimensions = ydim)

    refined<- xgeo |> sf::st_intersection(ygeo)
    refined$area<- refined |> sf::st_area()
    refined$weight<- refined$area / sf::st_area(ygeo)[refined$y]

    A<- Matrix::sparseMatrix(
            x = refined |> _$weight |> units::drop_units(),
            i = refined |> _$y,
            j = refined |> _$x,
            dims = c(ygeo |> nrow(), xgeo |> nrow())
        )
    
    put_geodim_first<- c(geodim, ydim |> seq_along() |> _[-geodim])
    put_first_to_geodim<- ydim |> seq_along() |> match(put_geodim_first)

    for( v in y |> names() |> seq_along() ) {
        y[[v]]<- x |>
            _[[v]] |>
            aperm(put_geodim_first) |>
            matrix(nrow = xgeo |> nrow()) |>
            (\(x) A %*% x)() |>
            array(dim = ydim |> dim() |> _[put_geodim_first]) |>
            aperm(put_first_to_geodim)
    }
    return(y)
}
