.filter_sort_mtx_files <- function(namesvec) {
    files <- .FEATURE_BC_MATRIX_FILES
    names(files) <- files
    res <- lapply(files, function(file) {
        namesvec[startsWith(namesvec, file)]
    })
    unlist(res)
}

.check_filter_mtx <- function(filelist) {
    afiles <- .filter_sort_mtx_files(names(filelist))
    if (!identical(names(afiles), .FEATURE_BC_MATRIX_FILES))
        stop(
            "'TENxFileList' does not contain the expected files:\n  ",
            .FEATURE_BC_MATRIX_FILES_PRINT
        )
    filelist[afiles]
}

.file_for_format <- function(fdir, format, ...) {
    if (identical(format, "h5")) {
        h5f <- file.path(fdir, paste0(.FEATURE_MATRIX_FILE_STUB, ".", format))
        path <- TENxH5(h5f, ...)
    } else if (identical(format, "mtx")) {
        mtxf <- file.path(
            fdir, paste0(.FEATURE_MATRIX_FILE_STUB, ".tar.gz")
        )
        path <- TENxFileList(mtxf)
    }
    path
}

.filter_h5_files <- function(path, format) {
    fname <- file.path(path, paste0(.FEATURE_MATRIX_FILE_STUB, ".", format))
    ish5file <- endsWith(names(path), fname)
    if (!any(ish5file))
        stop("The '", fname, "' file was not found.")
    path(path[ish5file])
}

#' @importFrom TENxIO TENxFileList TENxH5
.find_convert_resources <- function(path, format, ...) {
    if (!is(path, "TENxFileList"))
        path <- .file_for_format(path, format, ...)
    if (identical(format, "h5"))
        path <- .filter_h5_files(path, format)
    else if (identical(format, "mtx"))
        path <- .check_filter_mtx(path)

    path
}

#' @importFrom BiocIO FileForFormat
.filter_xenium_file <- function(path) {
    xf <- list.files(path, pattern = "\\.xenium$", full.names = TRUE)
    FileForFormat(xf)
}
