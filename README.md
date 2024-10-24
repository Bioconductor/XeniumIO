
# Introduction

The `XeniumIO` package provides functions to import 10X Genomics Xenium
Analyzer data into R. The package is designed to work with the output of
the Xenium Analyzer, which is a software tool that processes Visium
spatial gene expression data. The package provides functions to import
the output of the Xenium Analyzer into R, and to create a `TENxXenium`
object that can be used with other Bioconductor packages.

# Supported Formats

## TENxIO

The 10X suite of packages support multiple file formats. The following
table lists the supported file formats and the corresponding classes
that are imported into R.

| **Extension**       | **Class**     | **Imported as**                    |
|---------------------|---------------|------------------------------------|
| .h5                 | TENxH5        | SingleCellExperiment w/ TENxMatrix |
| .mtx / .mtx.gz      | TENxMTX       | SummarizedExperiment w/ dgCMatrix  |
| .tar.gz             | TENxFileList  | SingleCellExperiment w/ dgCMatrix  |
| peak_annotation.tsv | TENxPeaks     | GRanges                            |
| fragments.tsv.gz    | TENxFragments | RaggedExperiment                   |
| .tsv / .tsv.gz      | TENxTSV       | tibble                             |

## VisiumIO

| **Extension**  | **Class**          | **Imported as**   |
|----------------|--------------------|-------------------|
| spatial.tar.gz | TENxSpatialList    | DataFrame list \* |
| .parquet       | TENxSpatialParquet | tibble \*         |

## XeniumIO

| **Extension** | **Class** | **Imported as** |
|---------------|-----------|-----------------|
| .zarr.zip     | TENxZarr  | (TBD)           |

# GitHub Installation

``` r
BiocManager::install("Bioconductor/XeniumIO")
```

# Loading package

``` r
library(XeniumIO)
```

# XeniumIO

The `TENxXenium` class has a `metadata` slot for the `experiment.xenium`
file. The `resources` slot is a `TENxFileList` or `TENxH5` object
containing the cell feature matrix. The `coordNames` slot is a vector
specifying the names of the columns in the spatial data containing the
spatial coordinates. The `sampleId` slot is a scalar specifying the
sample identifier.

``` r
TENxXenium(
    resources = "path/to/matrix/folder/or/file",
    xeniumOut = "path/to/xeniumOut/folder",
    sample_id = "sample01",
    format = c("mtx", "h5"),
    boundaries_format = c("parquet", "csv.gz"),
    spatialCoordsNames = c("x_centroid", "y_centroid"),
    ...
)
```

The `format` argument specifies the format of the `resources` object,
either “mtx” or “h5”. The `boundaries_format` allows the user to choose
whether to read in the data using the `parquet` or `csv.gz` format.

# Example Folder Structure

Note that the `xeniumOut` unzipped folder must contain the following
files:

        *outs
        ├── cell_feature_matrix.h5
        ├── cell_feature_matrix.tar.gz
        |   ├── barcodes.tsv*
        |   ├── features.tsv*
        |   └── matrix.mtx*
        ├── cell_feature_matrix.zarr.zip
        ├── experiment.xenium
        ├── cells.csv.gz
        ├── cells.parquet
        ├── cells.zarr.zip
        [...]

Note that currently the `zarr` format is not supported as the
infrastructure is currently under development.

## Xenium class

The `resources` slot should either be the `TENxFileList` from the `mtx`
format or a `TENxH5` instance from an `h5` file. The boundaries can
either be a `TENxSpatialParquet` instance or a `TENxSpatialCSV`. These
classes are automatically instantiated by the constructor function.

``` r
showClass("TENxXenium")
#> Class "TENxXenium" [package "XeniumIO"]
#> 
#> Slots:
#>                                            
#> Name:                             resources
#> Class:               TENxFileList_OR_TENxH5
#>                                            
#> Name:                            boundaries
#> Class: TENxSpatialParquet_OR_TENxSpatialCSV
#>                                            
#> Name:                            coordNames
#> Class:                            character
#>                                            
#> Name:                              sampleId
#> Class:                            character
#>                                            
#> Name:                               colData
#> Class:                   TENxSpatialParquet
#>                                            
#> Name:                              metadata
#> Class:                           XeniumFile
```

## `import` method

The `import` method for a `TENxXenium` instance returns a
`SpatialExperiment` class object. Dispatch is only done on the `con`
argument. See `?BiocIO::import` for details on the generic. The `import`
function call is meant to be a simple call without much input. For more
details in the package, see `?TENxXenium`.

``` r
getMethod("import", c(con = "TENxXenium"))
#> Method Definition:
#> 
#> function (con, format, text, ...) 
#> {
#>     sce <- import(con@resources, ...)
#>     metadata <- import(con@metadata)
#>     coldata <- import(con@colData)
#>     SpatialExperiment::SpatialExperiment(assays = list(counts = assay(sce)), 
#>         rowData = rowData(sce), mainExpName = mainExpName(sce), 
#>         altExps = altExps(sce), sample_id = con@sampleId, colData = as(coldata, 
#>             "DataFrame"), spatialCoordsNames = con@coordNames, 
#>         metadata = list(experiment.xenium = metadata, polygons = import(con@boundaries)))
#> }
#> <bytecode: 0x5cf84e6ff3b8>
#> <environment: namespace:XeniumIO>
#> 
#> Signatures:
#>         con          format text 
#> target  "TENxXenium" "ANY"  "ANY"
#> defined "TENxXenium" "ANY"  "ANY"
```

# Importing an Example Xenium Dataset

The following code snippet demonstrates how to import a Xenium Analyzer
output into R. The `TENxXenium` object is created by specifying the path
to the `xeniumOut` folder. The `TENxXenium` object is then imported into
R using the `import` method for the `TENxXenium` class.

``` r
download.file(
    url = paste0(
        "https://cf.10xgenomics.com/samples/xenium/3.0.0/",
        "Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny/",
        "Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip"
    ),
    destfile =
        "~/data/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip"
)
unzip(
    zipfile =
        "~/data/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip",
    exdir = "~/data/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs",
    overwrite = FALSE
)
TENxXenium(
   xeniumOut = "~/data/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs"
) |> import()
#' class: SpatialExperiment 
#' dim: 8 36 
#' metadata(2): experiment.xenium polygons
#' assays(1): counts
#' rownames(8): DeprecatedCodeword_0321 DeprecatedCodeword_6781 ...
#'   DeprecatedCodeword_16059 DeprecatedCodeword_18533
#' rowData names(3): ID Symbol Type
#' colnames(36): aaamobki-1 aaclkaod-1 ... olbjkpjc-1 omjmdimk-1
#' colData names(13): cell_id transcript_counts ... segmentation_method
#'   sample_id
#' reducedDimNames(0):
#' mainExpName: Deprecated Codeword
#' altExpNames(5): Gene Expression Genomic Control Negative Control
#'   Codeword Negative Control Probe Unassigned Codeword
#' spatialCoords names(2) : x_centroid y_centroid
#' imgData names(0):
```

The dataset was obtained from the 10X Genomics website under the [X0A
v3.0
section](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/resources/xenium-example-data#test-data-v3-0)
and is a subset of the Xenium Prime 5K Mouse Pan Tissue & Pathways
Panel. The link to the data can be seen as the `url` input above and
shown below for completeness.

<https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip>

# Session Info

``` r
sessionInfo()
#> R version 4.4.1 Patched (2024-08-13 r87005)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.1 LTS
#> 
#> Matrix products: default
#> BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: America/New_York
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] XeniumIO_0.99.0             TENxIO_1.7.8               
#>  [3] SingleCellExperiment_1.27.2 SummarizedExperiment_1.35.4
#>  [5] Biobase_2.65.1              GenomicRanges_1.57.2       
#>  [7] GenomeInfoDb_1.41.2         IRanges_2.39.2             
#>  [9] S4Vectors_0.43.2            BiocGenerics_0.51.3        
#> [11] MatrixGenerics_1.17.0       matrixStats_1.4.1          
#> 
#> loaded via a namespace (and not attached):
#>  [1] utf8_1.2.4               SparseArray_1.5.44       lattice_0.22-6          
#>  [4] hms_1.1.3                digest_0.6.37            magrittr_2.0.3          
#>  [7] evaluate_1.0.1           grid_4.4.1               fastmap_1.2.0           
#> [10] jsonlite_1.8.9           Matrix_1.7-0             httr_1.4.7              
#> [13] fansi_1.0.6              VisiumIO_1.1.10          UCSC.utils_1.1.0        
#> [16] codetools_0.2-20         abind_1.4-8              cli_3.6.3               
#> [19] rlang_1.1.4              crayon_1.5.3             XVector_0.45.0          
#> [22] DelayedArray_0.31.14     yaml_2.3.10              BiocBaseUtils_1.7.3     
#> [25] S4Arrays_1.5.11          tools_4.4.1              tzdb_0.4.0              
#> [28] SpatialExperiment_1.15.1 GenomeInfoDbData_1.2.13  vctrs_0.6.5             
#> [31] R6_2.5.1                 magick_2.8.5             BiocIO_1.15.2           
#> [34] lifecycle_1.0.4          zlibbioc_1.51.1          pkgconfig_2.0.3         
#> [37] pillar_1.9.0             Rcpp_1.0.13              glue_1.8.0              
#> [40] xfun_0.48                tibble_3.2.1             knitr_1.48              
#> [43] rjson_0.2.23             htmltools_0.5.8.1        rmarkdown_2.28          
#> [46] readr_2.1.5              compiler_4.4.1
```
