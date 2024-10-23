#' A minimal class to represent Zarr files
#'
#' @docType class
#'
#' @return A `TENxZarr` class object
#'
#' @exportClass TENxZarr
.TENxZarr <- setClass(
    Class = "TENxZarr",
    contains = "TENxFile"
)

#' Import 10X data from the Zarr format
#'
#' @param resource `character(1)` The path to the Zarr file. Can be zipped.
#'
#' @return A `TENxZarr` class object
#'
#' @examples
#' zarrzip <- paste0(
#'     "~/data/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs/",
#'     "cell_feature_matrix.zarr.zip"
#' )
#' TENxZarr(zarrzip) |>
#'     import()
#'
#' @export
TENxZarr <- function(resource) {
    if (!is(resource, "TENxFile"))
        resource <- TENxFile(resource)
    compressed <- grepl("\\.zarr\\.zip$", path(resource), ignore.case = TRUE)
    .TENxZarr(resource, compressed = compressed)
}

.TENxUnzip <- function(con) {
    dir.create(tempdir <- tempfile(fileext = ".rarr"))
    unzip(path(con), exdir = tempdir)
    tempdir
}

#' @rdname TENxZarr
#' @exportMethod decompress
setMethod("decompress", "TENxZarr", function(manager, con, ...) {
    if (!identical(tools::file_ext(path(con)), "zip"))
        stop("The file '", basename(path(con)), "' is not a zip file.")
    .TENxUnzip(con)
})

#' @rdname TENxZarr
#' @importFrom Rarr read_zarr_array
#' @exportMethod import
setMethod("import", "TENxZarr", function(con, format, text, ...) {
    checkInstalled("Rarr")
    if (con@compressed)
        rarr <- decompress(con = con)
    else
        rarr <- path(con)
    zfold <- list.files(rarr, full.names = TRUE)
    if (!startsWith(basename(zfold), "cell_features"))
        stop(
            "'cell_features' directory not found in '", basename(path(con)), "'"
        )
    dfolds <- list.files(zfold, full.names = TRUE)
    names(dfolds) <- basename(dfolds)
    lapply(dfolds, read_zarr_array)
})
