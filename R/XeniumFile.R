#' @importClassesFrom TENxIO TENxFile
#' @exportClass XeniumFile
.XeniumFile <- setClass("XeniumFile",  contains = "TENxFile")

#' @export
XeniumFile <- function(resource) {
    stopifnot(isScalarCharacter(resource))
    .XeniumFile(resource = resource)
}

#' @importFrom BiocBaseUtils checkInstalled
#' @exportMethod import
setMethod("import", "XeniumFile", function(con, format, text, ...) {
    checkInstalled("jsonlite")
    jsonlite::fromJSON(path(con))
})
