# Scalable synthetic mixed-CDR benchmark for the sparse backend.
# Override the defaults with CDRGAM_SCALE_GROUPS and
# CDRGAM_SCALE_RESPONSES_PER_GROUP environment variables.

library(cdrgam)

group_count <- as.integer(Sys.getenv('CDRGAM_SCALE_GROUPS', '100'))
responses_per_group <- as.integer(Sys.getenv(
    'CDRGAM_SCALE_RESPONSES_PER_GROUP',
    '80'
))
item_count <- as.integer(Sys.getenv('CDRGAM_SCALE_ITEMS', '0'))
impulses_per_group <- max(50L, as.integer(responses_per_group * 0.8))
groups <- sprintf('subject_%04d', seq_len(group_count))
population_irf <- function(lag) 0.8 * exp(-1.8 * lag)
deviation_shape <- function(lag) {
    exp(-((lag - 0.7) / 0.3)^2) - 0.3 * exp(-1.5 * lag)
}
amplitudes <- 0.45 * sin(seq(0, 4 * pi, length.out=group_count))
intercepts <- 0.7 * cos(seq(0, 3 * pi, length.out=group_count))

impulse_streams <- vector('list', group_count)
response_streams <- vector('list', group_count)
for (i in seq_len(group_count)) {
    amplitude <- amplitudes[[i]]
    simulation <- simulate_cdr(
        list(x=function(lag) {
            population_irf(lag) + amplitude * deviation_shape(lag)
        }),
        n_impulses=impulses_per_group,
        n_responses=responses_per_group,
        duration=20,
        window=2,
        intercept=20 + intercepts[[i]],
        noise_sd=0.6,
        seed=8100 + i
    )
    simulation$impulses$subject <- groups[[i]]
    simulation$responses$subject <- groups[[i]]
    impulse_streams[[i]] <- simulation$impulses
    response_streams[[i]] <- simulation$responses
}
impulses <- do.call(rbind, impulse_streams)
responses <- do.call(rbind, response_streams)
responses$subject <- factor(responses$subject, levels=groups)
if (item_count > 0L) {
    set.seed(9101)
    items <- sprintf('item_%05d', seq_len(item_count))
    responses$item <- factor(
        sample(items, nrow(responses), replace=TRUE),
        levels=items
    )
    item_effects <- 0.4 * sin(seq(0, 6 * pi, length.out=item_count))
    responses$response <- responses$response +
        item_effects[as.integer(responses$item)]
}

model_formula <- if (item_count > 0L) {
    response ~ s(subject, bs='re') + s(item, bs='re') +
        irf(x, window=c(0, 2), k=10) +
        irf(x, window=c(0, 2), k=10, group=subject)
} else {
    response ~ s(subject, bs='re') +
        irf(x, window=c(0, 2), k=10) +
        irf(x, window=c(0, 2), k=10, group=subject)
}

design <- prepare_cdrgam(
    model_formula,
    impulses,
    responses,
    series='subject',
    history='auto',
    chunk_size=1000,
    quiet=FALSE
)
stopifnot(ncol(design$terms[[2]]$X) == 10L)

elapsed <- system.time({
    trace_optimizer <- identical(Sys.getenv('CDRGAM_SCALE_TRACE'), '1')
    gradient_method <- Sys.getenv('CDRGAM_SCALE_GRADIENT', 'auto')
    schur_method <- Sys.getenv('CDRGAM_SCALE_SCHUR', 'never')
    crossprod_chunk_size <- as.numeric(Sys.getenv(
        'CDRGAM_SCALE_CROSSPROD_CHUNK',
        '10000'
    ))
    supernodal_text <- Sys.getenv('CDRGAM_SCALE_SUPERNODAL', '')
    supernodal <- if (!nzchar(supernodal_text)) NULL else switch(
        tolower(supernodal_text),
        true=TRUE,
        false=FALSE,
        cholmod=NA,
        stop('CDRGAM_SCALE_SUPERNODAL must be true, false, or cholmod')
    )
    fit <- fit_cdrgam(
        design,
        backend='sparse',
        method='REML',
        solver_trace=trace_optimizer,
        sparse_control=c(
            list(
                gradient=gradient_method,
                schur=schur_method,
                crossprod_chunk_size=crossprod_chunk_size
            ),
            if (is.null(supernodal)) list() else list(supernodal=supernodal)
        )
    )
})

dense_grouped_bytes <- nrow(responses) * group_count * 10 * 8
results <- data.frame(
    groups=group_count,
    items=item_count,
    responses=nrow(responses),
    coefficients=length(coef(fit)),
    smoothing_parameters=length(fit$sp),
    exact_trace_rhs=sum(vapply(
        cdrgam:::.sparse_penalty_supports(fit$sparse$penalty_components),
        length,
        integer(1)
    )),
    gradient=fit$sparse$gradient,
    reml=fit$reml,
    objective_evaluations=fit$sparse$convergence$total_objective_evaluations,
    factor_class=fit$sparse$factor_class,
    factor_nonzeros=fit$sparse$factor_nonzeros,
    schur_group=if (is.null(fit$sparse$schur_group)) {
        NA_character_
    } else {
        fit$sparse$schur_group
    },
    schur_core=fit$sparse$schur_core_dimension,
    design_nonzeros=fit$sparse$nnzero,
    system_nonzeros=fit$sparse$system_nnzero,
    numeric_updates=fit$sparse$numeric_updates,
    full_factorizations=fit$sparse$full_factorizations,
    crossprod_chunks=fit$sparse$crossprod_chunks,
    streamed_grouped_terms=fit$sparse$streamed_grouped_terms,
    avoided_dense_grouped_gib=dense_grouped_bytes / 1024^3,
    prepared_design_mib=as.numeric(object.size(design)) / 1024^2,
    fitted_object_mib=as.numeric(object.size(fit)) / 1024^2,
    elapsed_seconds=unname(elapsed[['elapsed']])
)
print(results, row.names=FALSE)

lag <- seq(0, 2, length.out=101)
selected <- groups[unique(round(seq(1, group_count, length.out=8)))]
estimate <- estimate_irf(
    fit,
    term='x|subject',
    lag=lag,
    level=selected
)
stopifnot(all(is.finite(estimate$estimate)), all(is.finite(estimate$se)))
