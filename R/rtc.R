rtc.tree.clustering <- function(ds, nsamples=20, limits=NULL) {
    gctorture(TRUE)
    if (!is.data.frame(ds)) {
        stop("ds must be a data frame")
    }
    if (!is.integer(nsamples) && !is.numeric(nsamples) || nsamples <= 0) {
        stop("nsamples must be a positive integer")
    }
    if (!is.null(limits) && !is.list(limits)) {
        stop("limits must be a list")
    }
    if (dim(ds)[1] > 0 && !is.null(limits) && length(limits) != dim(ds)[2]) {
        stop("limits must have the same length as the number of parameters")
    }
    for (range in limits) {
        if (!is.numeric(range)) {
            stop("invalid limits")
        }
    }
    N <- dim(ds)[1]
    K <- dim(ds)[2]

    opts <- list()
    opts$nsamples <- as.integer(nsamples)
    segmentation <- .Call(
        "rtc_tree_clustering",
        ds,
        N,
        K,
        limits,
        opts,
        PACKAGE="rtc"
    )
    gctorture(FALSE)
    segmentation
}
