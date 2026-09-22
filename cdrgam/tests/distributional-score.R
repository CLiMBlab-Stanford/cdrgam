library(cdrgam)

set.seed(5291)
impulses <- data.frame(
    time=sort(stats::runif(100, 0, 30)),
    location_signal=stats::rnorm(100),
    scale_signal=stats::rnorm(100)
)
responses <- data.frame(
    time=sort(stats::runif(140, 0.5, 30)),
    response=stats::rnorm(140)
)
design <- prepare_cdrgam(
    list(
        location=response ~
            irf(location_signal, window=c(0, 1.5), k_l=5) - irf(1),
        scale=~ irf(scale_signal, window=c(0, 1.5), k_l=5) - irf(1)
    ),
    impulses,
    responses,
    history='ragged',
    chunk_size=31,
    quiet=TRUE
)
assemblies <- lapply(design$parameters, function(parameter_design) {
    cdrgam:::.cdrgam_distributional_sparse_setup(
        parameter_design,
        list(crossprod_chunk_size=31)
    )
})
family <- cdrgam_family('gaulss')
log_sp <- log(unlist(lapply(
    assemblies,
    cdrgam:::.cdrgam_sparse_initial_sp
), use.names=FALSE))

criterion <- function(parameters, retain=FALSE) {
    sp <- exp(parameters)
    solution <- cdrgam:::.cdrgam_gaulss_sparse_fixed(
        assemblies,
        sp,
        cdrgam:::.cdrgam_gaulss_b(family)
    )
    counts <- vapply(
        assemblies,
        function(assembly) length(assembly$penalty_components),
        integer(1)
    )
    penalty_determinant <- 0
    offset <- 0L
    for (parameter in names(assemblies)) {
        count <- counts[[parameter]]
        indices <- offset + seq_len(count)
        penalty_determinant <- penalty_determinant +
            cdrgam:::.sparse_penalty_logdet(
                assemblies[[parameter]]$blocks,
                sp[indices]
            )$value
        offset <- offset + count
    }
    value <- 2 * solution$objective +
        cdrgam:::.cdr_factor_logdet(solution$factor) - penalty_determinant
    if (retain) list(value=value, solution=solution) else value
}

retained <- criterion(log_sp, retain=TRUE)
analytic <- cdrgam:::.cdrgam_gaulss_sparse_score(
    assemblies,
    retained$solution,
    exp(log_sp),
    family
)
step <- 1e-4
finite <- vapply(seq_along(log_sp), function(i) {
    plus <- minus <- log_sp
    plus[[i]] <- plus[[i]] + step
    minus[[i]] <- minus[[i]] - step
    (criterion(plus) - criterion(minus)) / (2 * step)
}, numeric(1))
stopifnot(max(abs(analytic - finite)) < 2e-4)

cat('Sparse distributional score checks passed.\n')
