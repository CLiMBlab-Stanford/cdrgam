library(cdrgam)

# Structural sparse keys must remain exact beyond the 32-bit overflow point
# reached by realistic crossed nonlinear random-IRF models.
large_keys <- cdrgam:::.sparse_matrix_keys(
    c(0L, 49999L, 0L, 49999L),
    c(0L, 0L, 49999L, 49999L),
    50000L
)
stopifnot(
    identical(large_keys, c(0, 49999, 2499950000, 2499999999)),
    length(unique(large_keys)) == length(large_keys)
)

# Nearly collinear impulse predictors, deliberately shuffled responses, and
# very different accumulation chunk sizes must lead to the same optimum.
simulation <- simulate_cdr(
    list(
        x1=function(lag) 0.8 * exp(-1.6 * lag),
        x2=function(lag) -0.4 * exp(-((lag - 0.7) / 0.3)^2)
    ),
    n_impulses=250,
    n_responses=320,
    duration=35,
    window=2,
    noise_sd=0.6,
    seed=7711
)
simulation$impulses$x2 <- simulation$impulses$x1 +
    0.02 * simulation$impulses$x2
formula <- response ~
    irf(x1, window=c(0, 2), k=7) +
    irf(x2, window=c(0, 2), k=7) - irf(1)
design <- prepare_cdrgam(
    formula,
    simulation$impulses,
    simulation$responses,
    history='ragged',
    quiet=TRUE
)
small_chunks <- fit_cdrgam(
    design,
    backend='sparse',
    method='REML',
    sparse_control=list(
        gradient='exact',
        crossprod_chunk_size=17,
        restarts=2
    )
)
one_chunk <- fit_cdrgam(
    design,
    backend='sparse',
    method='REML',
    sparse_control=list(
        gradient='exact',
        crossprod_chunk_size=100000,
        restarts=2
    )
)
native <- fit_cdrgam(design, backend='mgcv', engine='gam', method='REML')
block <- fit_cdrgam(design, backend='block', method='REML')
stopifnot(isTRUE(small_chunks$converged), isTRUE(one_chunk$converged))
stopifnot(small_chunks$sparse$convergence$restart_count == 2L)
stopifnot(length(small_chunks$sparse$convergence$restart_objectives) == 3L)
stopifnot(!small_chunks$sparse$convergence$global_optimum_certified)
stopifnot(max(abs(coef(small_chunks) - coef(one_chunk))) < 1e-6)
stopifnot(max(abs(fitted(small_chunks) - fitted(one_chunk))) < 1e-6)
# Under this severe collinearity mgcv settles at a secondary REML stationary
# point. The independent dense objective agrees with the sparse optimum.
stopifnot(max(abs(fitted(block) - fitted(small_chunks))) < 1e-6)
stopifnot(all(is.finite(fitted(native))))

set.seed(7712)
permutation <- sample.int(nrow(simulation$responses))
shuffled_design <- prepare_cdrgam(
    formula,
    simulation$impulses,
    simulation$responses[permutation, ],
    history='ragged',
    quiet=TRUE
)
shuffled <- fit_cdrgam(
    shuffled_design,
    backend='sparse',
    method='REML',
    sparse_control=list(gradient='exact', crossprod_chunk_size=19)
)
stopifnot(isTRUE(shuffled$converged))
stopifnot(max(abs(coef(small_chunks) - coef(shuffled))) < 1e-5)
stopifnot(max(abs(
    fitted(small_chunks)[permutation] - fitted(shuffled)
)) < 1e-5)

# Strongly unbalanced response-aligned groups exercise sparse random
# intercepts and compact grouped IRFs, including a group with very few rows.
group_sizes <- c(120L, 55L, 22L, 10L)
groups <- paste0('g', seq_along(group_sizes))
impulse_streams <- response_streams <- vector('list', length(groups))
for (i in seq_along(groups)) {
    generated <- simulate_cdr(
        list(x=function(lag) {
            (0.7 + 0.08 * i) * exp(-1.5 * lag)
        }),
        n_impulses=max(30L, group_sizes[[i]]),
        n_responses=group_sizes[[i]],
        duration=25,
        window=2,
        intercept=10 + i / 3,
        noise_sd=0.5,
        seed=7800 + i
    )
    generated$impulses$group <- groups[[i]]
    generated$responses$group <- groups[[i]]
    impulse_streams[[i]] <- generated$impulses
    response_streams[[i]] <- generated$responses
}
impulses <- do.call(rbind, impulse_streams)
responses <- do.call(rbind, response_streams)
responses$group <- factor(responses$group, levels=groups)
unbalanced <- prepare_cdrgam(
    response ~ s(group, bs='re') +
        irf(x, window=c(0, 2), k=6) +
        irf(x, window=c(0, 2), k=6, group=group) - irf(1),
    impulses,
    responses,
    series='group',
    history='ragged',
    quiet=TRUE
)
unbalanced_native <- fit_cdrgam(
    unbalanced,
    backend='mgcv',
    engine='gam',
    method='REML'
)
unbalanced_sparse <- fit_cdrgam(
    unbalanced,
    backend='sparse',
    method='REML',
    sparse_control=list(crossprod_chunk_size=23)
)
stopifnot(isTRUE(unbalanced_sparse$converged))
stopifnot(unbalanced_sparse$sparse$streamed_grouped_terms == 1L)
stopifnot(max(abs(
    fitted(unbalanced_native) - fitted(unbalanced_sparse)
)) < 1e-3)

# Perfectly redundant ordinary parametric columns are benign: they are
# deterministically removed, recorded, and remain valid for prediction.
rank_responses <- base::transform(
    simulation$responses,
    z=seq_len(nrow(simulation$responses))
)
rank_design <- prepare_cdrgam(
    response ~ z + I(z) + irf(x1, window=c(0, 2), k=7) - irf(1),
    simulation$impulses,
    rank_responses,
    history='ragged',
    quiet=TRUE
)
rank_warning <- NULL
rank_fit <- withCallingHandlers(
    fit_cdrgam(rank_design, backend='sparse', method='REML'),
    warning=function(w) {
        if (grepl('aliased ordinary parametric', conditionMessage(w))) {
            rank_warning <<- conditionMessage(w)
            invokeRestart('muffleWarning')
        }
    }
)
rank_reference <- prepare_cdrgam(
    response ~ z + irf(x1, window=c(0, 2), k=7) - irf(1),
    simulation$impulses,
    rank_responses,
    history='ragged',
    quiet=TRUE
)
rank_reference <- fit_cdrgam(rank_reference, backend='sparse', method='REML')
rank_native <- suppressWarnings(fit_cdrgam(
    rank_design,
    backend='mgcv',
    engine='gam',
    method='REML'
))
rank_block <- suppressWarnings(fit_cdrgam(
    rank_design,
    backend='block',
    method='REML'
))
stopifnot(grepl('I\\(z\\)', rank_warning))
stopifnot(identical(rank_fit$cdrgam$rank$parametric$dropped, 'I(z)'))
stopifnot(identical(rank_block$cdrgam$rank$parametric$dropped, 'I(z)'))
stopifnot(identical(
    rank_native$cdrgam$identifiability$global$parametric$dropped,
    'I(z)'
))
stopifnot(inherits(summary(rank_native), 'summary.gam'))
stopifnot(max(abs(fitted(rank_fit) - fitted(rank_reference))) < 1e-7)
rank_prediction <- list(
    impulses=simulation$impulses,
    responses=rank_responses[names(rank_responses) != 'response']
)
stopifnot(max(abs(predict(rank_fit, rank_prediction) - fitted(rank_fit))) < 1e-8)
stopifnot(max(abs(
    predict(rank_native, rank_prediction) - fitted(rank_native)
)) < 1e-8)

# Confounding that spans two complete IRF terms is scientifically meaningful.
# It errors by default, while explicit unattended-run policies are honored.
confounded_impulses <- simulation$impulses
confounded_impulses$x_duplicate <- confounded_impulses$x1
confounded_design <- prepare_cdrgam(
    response ~ irf(x1, window=c(0, 2), k=7) +
        irf(x_duplicate, window=c(0, 2), k=7) - irf(1),
    confounded_impulses,
    simulation$responses,
    history='ragged',
    quiet=TRUE
)
rank_error <- tryCatch(
    fit_cdrgam(confounded_design, backend='sparse', method='REML'),
    error=function(e) conditionMessage(e)
)
drop_error <- tryCatch(
    fit_cdrgam(
        confounded_design,
        backend='sparse',
        method='REML',
        rank_action='drop'
    ),
    error=function(e) conditionMessage(e)
)
minimum_norm <- suppressWarnings(fit_cdrgam(
    confounded_design,
    backend='sparse',
    method='REML',
    rank_action='minimum_norm'
))
penalized <- suppressWarnings(fit_cdrgam(
    confounded_design,
    backend='sparse',
    method='REML',
    rank_action='penalize',
    rank_penalty=1e-5
))
stopifnot(grepl('rank_action="minimum_norm"', rank_error, fixed=TRUE))
stopifnot(grepl('spans smooth or IRF terms', drop_error, fixed=TRUE))
stopifnot(identical(
    minimum_norm$cdrgam$rank$resolution,
    'minimum-norm-ridge-approximation'
))
stopifnot(identical(
    penalized$cdrgam$rank$resolution,
    'fixed-ridge-penalty'
))
stopifnot(all(is.finite(fitted(minimum_norm))))
stopifnot(all(is.finite(fitted(penalized))))
minimum_indices <- lapply(
    minimum_norm$cdrgam$terms,
    function(term) term$coefficient_index
)
stopifnot(max(abs(
    coef(minimum_norm)[minimum_indices[[1L]]] -
        coef(minimum_norm)[minimum_indices[[2L]]]
)) < 1e-5)

block_minimum_norm <- suppressWarnings(fit_cdrgam(
    confounded_design,
    backend='block',
    method='REML',
    rank_action='minimum_norm'
))
stopifnot(identical(
    block_minimum_norm$cdrgam$rank$resolution,
    'minimum-norm-ridge-approximation'
))
stopifnot(all(is.finite(fitted(block_minimum_norm))))
