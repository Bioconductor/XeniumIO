#' @include XeniumFile.R
#' @importClassesFrom TENxIO TENxFileList TENxH5
setClassUnion("TENxFileList_OR_TENxH5", members = c("TENxFileList", "TENxH5"))

setClassUnion(
    "TENxSpatialParquet_OR_TENxSpatialCSV",
    members = c("TENxSpatialParquet", "TENxSpatialCSV")
)

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
#' @slot colData `TENxSpatialParquet` A [TENxSpatialParquet] object containing
#'  the spatial coordinates data.
#'
#' @slot metadata `XeniumFile` A [XeniumFile] object containing the metadata
#'  information.
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
        boundaries = "TENxSpatialParquet_OR_TENxSpatialCSV",
        coordNames = "character",
        sampleId = "character",
        colData = "TENxSpatialParquet",
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
#' @param boundaries_format `character(1)` Either "parquet" or "csv.gz" to
#'   specify the file extension of the boundaries file. Default is "parquet".
#'
#' @importFrom methods is new
#' @importFrom BiocBaseUtils isScalarCharacter
#'
#' @examples
#' if (interactive()) {
#'     download.file(
#'         url = paste0(
#'             "https://cf.10xgenomics.com/samples/xenium/3.0.0/",
#'             "Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny/",
#'             "Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip"
#'         ),
#'         destfile =
#'             "~/data/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip"
#'     )
#'     unzip(
#'         zipfile =
#'             "~/data/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip",
#'         exdir = "~/data/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs",
#'         overwrite = FALSE
#'     )
#'     TENxXenium(
#'        xeniumOut = "~/data/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs"
#'     ) |> import()
#' }
#' @export
TENxXenium <- function(
    resources,
    xeniumOut,
    sample_id = "sample01",
    format = c("mtx", "h5"),
    boundaries_format = c("parquet", "csv.gz"),
    spatialCoordsNames = c("x_centroid", "y_centroid"),
    ...
) {
    format <- match.arg(format)
    boundaries_format <- match.arg(boundaries_format)

    if (!missing(xeniumOut)) {
        if (isScalarCharacter(xeniumOut) && !dir.exists(xeniumOut))
            stop(
                "The '", xeniumOut, "' directory was not found.",
                "\n Verify 'xeniumOut' input directory.",
                call. = FALSE
            )
        resources <- .file_for_format(xeniumOut, format, ...)
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
        xeniumOut <- dirname(path(resources))
    }

    xeniumfile <- .filter_xenium_file(xeniumOut)
    parqfile <- .filter_parquet_file(xeniumOut)
    bounds <- .boundaries_for_format(xeniumOut, boundaries_format)

    .TENxXenium(
        resources = resources,
        boundaries = bounds,
        coordNames = spatialCoordsNames,
        sampleId = sample_id,
        colData = parqfile,
        metadata = xeniumfile
    )
}

.validTENxXenium <- function(object) {
    isFL <- is(object@resources, "TENxFileList_OR_TENxH5")
    isSP <- is(object@colData, "TENxSpatialParquet")
    isXF <- is(object@metadata, "XeniumFile")
    if (all(isFL, isSP, isXF))
        TRUE
    else if (!isFL)
        "'resources' component is not of TENxFileList or TENxH5 class"
    else if (!isSP)
        "'colData' component is not of TENxSpatialParquet class"
    else if (!isXF)
        "'metadata' component is not of XeniumFile class"
}

S4Vectors::setValidity2("TENxXenium", .validTENxXenium)

# import TENxXenium method ------------------------------------------------

#' @describeIn TENxXenium-class Import Xenium Analyzer data
#'
#' @inheritParams BiocIO::import
#'
#' @importFrom BiocIO import
#' @importFrom methods as
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SummarizedExperiment assay assays rowData colData
#' @importFrom SingleCellExperiment mainExpName altExps
#'
#' @exportMethod import
setMethod("import", "TENxXenium", function(con, format, text, ...) {
    sce <- import(con@resources, ...)
    metadata <- import(con@metadata)
    coldata <- import(con@colData)

    ## TODO: mainExpName and altExpNames are lost when SCE sent to constructor
    SpatialExperiment::SpatialExperiment(
        assays = list(counts = assay(sce)),
        rowData = rowData(sce),
        mainExpName = mainExpName(sce),
        altExps = altExps(sce),
        sample_id = con@sampleId,
        colData = as(coldata, "DataFrame"),
        spatialCoordsNames = con@coordNames,
        metadata = list(
            experiment.xenium = metadata,
            polygons = import(con@boundaries)
        )
    )
})
