#' @include XeniumFile.R
#' @importClassesFrom TENxIO TENxFileList TENxH5
setClassUnion("TENxFileList_OR_TENxH5", members = c("TENxFileList", "TENxH5"))

#' @docType class
#'
#' @title A class to represent Xenium output data
#'
#' @description This class is a composed class of [TENxFileList] which can
#'   contain a list of [TENxFile] objects for the cell-feature matrix. It is
#'   meant to handle a single Xenium sample from 10X Genomics.
#'
#' @slot resources A [TENxFileList] or [TENxH5] object containing the cell
#'   feature matrix.
#'
#' @slot coordNames `character()` A vector specifying the names
#'   of the columns in the spatial data containing the spatial coordinates.
#'
#' @slot sampleId `character(1)` A scalar specifying the sample identifier.
#'
#' @return A [SpatialExperiment] object
#'
#' @seealso <https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs>
#'
#' @exportClass TENxXenium
.TENxXenium <- setClass(
    Class = "TENxXenium",
    slots = c(
        resources = "TENxFileList_OR_TENxH5",
        coordNames = "character",
        sampleId = "character",
        metadata = "XeniumFile"
    )
)

.FEATURE_BC_MATRIX_FILES <- c("barcodes.tsv", "features.tsv", "matrix.mtx")

.FEATURE_BC_MATRIX_FILES_PRINT <- paste0(
    sQuote(.FEATURE_BC_MATRIX_FILES), collapse = ", "
)

.FEATURE_MATRIX_FILE_STUB <- "cell_feature_matrix"

#' @rdname TENxXenium-class
#'
#' @inheritParams VisiumIO::TENxVisium
#'
#' @param xeniumOut `character(1)` The path to the Xenium output directory.
#'
#' @importFrom methods is new
#' @importFrom BiocBaseUtils isScalarCharacter
#'
#' @export
TENxXenium <- function(
    resources,
    xeniumOut,
    sample_id = "sample01",
    format = c("mtx", "h5"),
    spatialCoordsNames = c("x_centroid", "y_centroid"),
    ...
) {
    format <- match.arg(format)

    if (!missing(xeniumOut)) {
        if (isScalarCharacter(xeniumOut) && !dir.exists(xeniumOut))
            stop(
                "The '", xeniumOut, "' directory was not found.",
                "\n Verify 'xeniumOut' input directory.",
                call. = FALSE
            )
        resources <- .file_for_format(xeniumOut, format)
    } else {
        stopifnot(
            (isScalarCharacter(resources) && file.exists(resources)) ||
                is(resources, "TENxFileList_OR_TENxH5")
        )
        if (
            !is(resources, "TENxFileList_OR_TENxH5") &&
            identical(tools::file_ext(resources), "h5")
        )
            resources <- TENxH5(resources, ...)
        else if (isScalarCharacter(resources))
            resources <- TENxFileList(resources, ...)
        xeniumOut <- dirname(resources)
    }

    xeniumfile <- .filter_xenium_file(xeniumOut)

    .TENxXenium(
        resources = resources,
        coordNames = spatialCoordsNames,
        sampleId = sample_id,
        metadata = xeniumfile
    )
}

#' @importFrom S4Vectors metadata
setMethod("import", "TENxXenium", function(con, format, text, ...) {
    sce <- import(con@resources)
    metadata <- import(con@metadata)
    metadata(sce) <- metadata
    sce
})
