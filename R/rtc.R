rtc.tree.clustering <- function(ds, nsamples=20, limits=NULL) {
    if (!is.data.frame(ds)) {
        stop("ds must be a data frame")
    }
    if (!is.integer(nsamples) && !is.numeric(nsamples) || nsamples <= 0) {
        stop("nsamples must be a positive integer")
    }
    if (!is.null(limits) && !is.list(limits)) {
        stop("limits must be a list")
    }
    for (range in limits) {
        if (!is.numeric(range)) {
            stop("invalid limits")
        }
    }
    opts <- list()
    opts$nsamples <- as.integer(nsamples)
    .Call(
        "rtc_tree_clustering",
        ds,
        dim(ds)[1],
        dim(ds)[2],
        limits,
        opts,
        PACKAGE="rtc"
    )
}
