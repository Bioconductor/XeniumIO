zipfile <- paste0(
    "https://mghp.osn.xsede.org/bir190004-bucket01/BiocXenDemo/",
    "Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip"
)
destfile <- XeniumIO:::.cache_url_file(zipfile)
outfold <- file.path(
    tempdir(), tools::file_path_sans_ext(basename(zipfile))
)
if (!dir.exists(outfold))
    dir.create(outfold, recursive = TRUE)
unzip(
    zipfile = destfile, exdir = outfold, overwrite = FALSE
)
semtx <- TENxXenium(xeniumOut = outfold) |>
    import()

expect_true(
    is(assay(semtx), "dgCMatrix")
)

expect_identical(
    colData(semtx)[["sample_id"]] |> unique(),
    "sample01"
)

expect_warning(
    seh5 <- TENxXenium(xeniumOut = outfold, format = "h5") |> import()
)

expect_true(
    is(assay(seh5), "DelayedMatrix")
)

expect_identical(
    colData(seh5)[["sample_id"]] |> unique(),
    "sample01"
)

expect_silent(
    TENxXenium(xeniumOut = outfold, boundaries_format = "csv.gz") |> import()
)

cellsnames <- VisiumIO::TENxSpatialCSV(file.path(outfold, "cells.csv.gz")) |>
    import() |>
    names() |>
    tail(-2L)

colDataNames <- colData(
    TENxXenium(
        xeniumOut = outfold, boundaries_format = "csv.gz"
    ) |> import()
) |>
    names()

expect_true(
    all(
        cellsnames %in% colDataNames
    )
)
