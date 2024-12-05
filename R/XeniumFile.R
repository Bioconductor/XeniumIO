#' @docType class
#'
#' @title A minimal class to represent Xenium metadata
#'
#' @description This class is a minimal class to represent Xenium metadata.
#' It is dedicated to importing `experiment.xenium` metadata files. It uses
#' the `jsonlite` package to import the metadata.
#'
#' @return A [XeniumFile] object
#'
#' @importClassesFrom TENxIO TENxFile
#' @examples
#' showClass("XeniumFile")
#'
#' @exportClass XeniumFile
.XeniumFile <- setClass("XeniumFile",  contains = "TENxFile")

#' @rdname XeniumFile-class
#'
#' @param resource `character(1)` The path to the Xenium metadata file.
#'
#' @export
XeniumFile <- function(resource) {
    stopifnot(isScalarCharacter(resource))
    .XeniumFile(resource = resource)
}

#' @describeIn XeniumFile-class Import Xenium metadata
#'
#' @inheritParams BiocIO::import
#'
#' @importFrom BiocGenerics path
#' @exportMethod import
setMethod("import", "XeniumFile", function(con, format, text, ...) {
    jsonlite::fromJSON(path(con))
})
