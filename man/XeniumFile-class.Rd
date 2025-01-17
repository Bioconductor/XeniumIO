% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XeniumFile.R
\docType{class}
\name{XeniumFile-class}
\alias{XeniumFile-class}
\alias{.XeniumFile}
\alias{XeniumFile}
\alias{import,XeniumFile,ANY,ANY-method}
\title{A minimal class to represent Xenium metadata}
\usage{
XeniumFile(resource)

\S4method{import}{XeniumFile,ANY,ANY}(con, format, text, ...)
}
\arguments{
\item{resource}{\code{character(1)} The path to the Xenium metadata file.}

\item{con}{The connection from which data is loaded or to which data is
saved. If this is a \code{character} vector, it is assumed to be a file name and
a corresponding file connection is created and then closed after exporting
the object. If it is a \link[BiocIO:BiocFile-class]{BiocFile} derivative, the
data is loaded from or saved to the underlying resource. If missing, the
function will return the output as a character vector, rather than writing
to a connection.}

\item{format}{The format of the output. If missing and \code{con} is a
file name, the format is derived from the file extension. This argument is
unnecessary when \code{con} is a derivative of \link[BiocIO:BiocFile-class]{BiocFile}.}

\item{text}{If \code{con} is missing, this can be a character vector directly
providing the string data to import.}

\item{...}{Parameters to pass to the format-specific method.}
}
\value{
A \link{XeniumFile} object
}
\description{
This class is a minimal class to represent Xenium metadata.
It is dedicated to importing \code{experiment.xenium} metadata files. It uses
the \code{jsonlite} package to import the metadata.
}
\section{Functions}{
\itemize{
\item \code{import(con = XeniumFile, format = ANY, text = ANY)}: Import Xenium metadata

}}
\examples{
showClass("XeniumFile")

}
