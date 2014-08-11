rtc.tree.clustering <- function(
    ds,
    nsamples=20,
    burnin=0,
    maxiter=0,
    limits=NULL,
    fragment.size=NULL,
    max.segments=0,
    nruns=1
) {
    # gctorture(TRUE)

    N <- dim(ds)[1]
    K <- dim(ds)[2]

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
    if (!is.null(fragment.size) && !is.numeric(fragment.size)) {
        stop("fragment.size must be a numeric vector")
    }
    if (!is.null(fragment.size) && length(fragment.size) != K) {
        stop("fragment.size must have the same length as the number of parameters")
    }
    if (!is.numeric(max.segments) || length(max.segments) != 1) {
        stop("max.segments must be an integer")
    }
    if (!is.numeric(burnin) || length(burnin) != 1) {
        stop("burnin must be an integer")
    }
    if (!is.numeric(nruns) || length(nruns) != 1 || nruns <= 0) {
        stop("nruns must be a positive integer")
    }
    if (!is.numeric(maxiter) || length(maxiter) != 1 || maxiter < 0) {
        stop("maxiter must be a non-negative integer")
    }

    opts <- list()
    opts$nsamples <- as.integer(nsamples)
    opts$max.segments <- as.integer(max.segments)
    opts$maxiter <- as.integer(maxiter)
    opts$burnin <- as.integer(burnin)

    results <- lapply(seq(nruns), function(run) {
        .Call(
            "rtc_tree_clustering",
            ds,
            N,
            K,
            limits,
            as.numeric(fragment.size),
            opts,
            PACKAGE="rtc"
        )
    })
    # gctorture(FALSE)
    unlist(results, recursive=FALSE)
}
