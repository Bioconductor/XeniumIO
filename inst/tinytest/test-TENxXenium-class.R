destfile <- XeniumIO:::.cache_url_file(
    url = paste0(
        "https://cf.10xgenomics.com/samples/xenium/3.0.0/",
        "Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny/",
        "Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip"
    )
)
outfile <- file.path(
    tempdir(), "Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs"
)
if (!dir.exists(outfile))
    dir.create(outfile, recursive = TRUE)
unzip(
    zipfile = destfile, exdir = outfile, overwrite = FALSE
)
semtx <- TENxXenium(xeniumOut = outfile) |>
    import()

expect_true(
    is(assay(semtx), "dgCMatrix")
)

expect_identical(
    colData(semtx)[["sample_id"]] |> unique(),
    "sample01"
)

expect_warning(
    seh5 <- TENxXenium(xeniumOut = outfile, format = "h5") |> import()
)

expect_true(
    is(assay(seh5), "DelayedMatrix")
)

expect_identical(
    colData(seh5)[["sample_id"]] |> unique(),
    "sample01"
)

expect_silent(
    TENxXenium(xeniumOut = outfile, boundaries_format = "csv.gz") |> import()
)

cellsnames <- TENxSpatialCSV(file.path(outfile, "cells.csv.gz")) |>
    import() |>
    names() |>
    tail(-2L)

colDataNames <- colData(
    TENxXenium(
        xeniumOut = outfile, boundaries_format = "csv.gz"
    ) |> import()
) |>
    names()

expect_true(
    all(
        cellsnames %in% colDataNames
    )
)
