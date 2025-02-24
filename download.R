download.url <- function(url, dir) {
    destfile <- file.path(dir, basename(url))
    if (!file.exists(destfile))
        download.file(url, destfile=desfile)
}

dir.create(output.dir <- ".")

full.dir <- file.path(output.dir, "data-full")
data.dir <- file.path(output.dir, "data")

## download the dataset from GEO
if (!file.exists(full.dir)) {
    dir.create(full.dir)
    data.url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53261/matrix/"
    download.url(file.path(data.url, "GSE53261-GPL11154_series_matrix.txt.gz"), full.dir)
    download.url(file.path(data.url, "GSE53261-GPL13534_series_matrix.txt.gz"), full.dir)
    download.url(file.path(data.url, "GSE53261-GPL6883_series_matrix.txt.gz"), full.dir)
    download.url(file.path(data.url, "GSE53261-GPL6984_series_matrix.txt.gz"), full.dir)
}