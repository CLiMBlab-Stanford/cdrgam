library(cdrgam)

set.seed(90210)
impulses <- data.frame(
    time=sort(stats::runif(80L, 0, 20)),
    location_signal=stats::rnorm(80L),
    scale_signal=stats::rnorm(80L)
)
responses <- data.frame(
    time=sort(stats::runif(110L, 0.5, 20)),
    response=stats::rnorm(110L)
)
formulas <- list(
    location=response ~
        irf(location_signal, window=c(0, 1), k_l=4) - irf(1),
    scale=~ irf(scale_signal, window=c(0, 1), k_l=4) - irf(1)
)
design <- prepare_cdrgam(
    formulas,
    impulses,
    responses,
    history='ragged',
    chunk_size=31L,
    quiet=TRUE
)

checkpoint <- tempfile('cdrgam-distributional-', fileext='.rds')
interrupted <- tryCatch(
    cdrgam.fit(
        design,
        family='gaulss',
        backend='sparse',
        checkpoint=checkpoint,
        solver_trace=function(event) {
            if (identical(event$event, 'outer initial')) {
                stop('simulated exact interruption')
            }
        },
        sparse_control=list(gradient='exact', optimizer_maxit=2L)
    ),
    error=function(error) conditionMessage(error)
)
if (!file.exists(checkpoint)) stop(interrupted)
partial <- readRDS(checkpoint)
stopifnot(
    identical(interrupted, 'simulated exact interruption'),
    identical(partial$backend, 'distributional-sparse'),
    identical(partial$stage, 'optimization'),
    identical(partial$phase, 'exact'),
    partial$score_evaluations >= 1L,
    is.finite(partial$best_criterion),
    identical(
        partial$optimizer_state$optimizer,
        'safeguarded_outer_bfgs'
    )
)

resume_events <- list()
resumed <- cdrgam.fit(
    design,
    family='gaulss',
    backend='sparse',
    checkpoint=checkpoint,
    solver_trace=function(event) {
        resume_events[[length(resume_events) + 1L]] <<- event
    },
    sparse_control=list(gradient='exact', optimizer_maxit=2L)
)
complete <- readRDS(checkpoint)
recovered <- cdrgam.fit(
    design,
    family='gaulss',
    backend='sparse',
    checkpoint=checkpoint,
    sparse_control=list(gradient='exact', optimizer_maxit=2L)
)
mismatch <- tryCatch(
    {
        cdrgam.fit(
            design,
            family='gaulss',
            backend='sparse',
            checkpoint=checkpoint,
            sparse_control=list(gradient='exact', optimizer_maxit=3L)
        )
        NA_character_
    },
    error=function(error) conditionMessage(error)
)
stopifnot(
    inherits(resumed, 'cdrgam_distributional_sparse'),
    identical(complete$stage, 'complete'),
    !is.null(complete$optimization),
    any(vapply(resume_events, function(event) {
        identical(event$event, 'checkpoint resumed')
    }, logical(1))),
    any(vapply(resume_events, function(event) {
        identical(event$event, 'outer resumed')
    }, logical(1))),
    identical(
        complete$optimization$optimizer_state$optimizer,
        'safeguarded_outer_bfgs'
    ),
    max(abs(stats::coef(recovered) - stats::coef(resumed))) < 1e-4,
    abs(recovered$reml - resumed$reml) < 1e-6,
    !is.na(mismatch),
    grepl('does not match this model', mismatch, fixed=TRUE)
)

hybrid_checkpoint <- tempfile('cdrgam-distributional-hybrid-', fileext='.rds')
hybrid_interrupted <- tryCatch(
    cdrgam.fit(
        design,
        family='gaulss',
        backend='sparse',
        checkpoint=hybrid_checkpoint,
        solver_trace=function(event) {
            if (identical(event$event, 'outer stochastic score complete')) {
                stop('simulated stochastic interruption')
            }
        },
        sparse_control=list(
            gradient='hybrid',
            gradient_probes=4L,
            optimizer_maxit=1L,
            optimizer_gradient_tolerance=1e6
        )
    ),
    error=function(error) conditionMessage(error)
)
if (!file.exists(hybrid_checkpoint)) stop(hybrid_interrupted)
hybrid_partial <- readRDS(hybrid_checkpoint)
stopifnot(
    identical(hybrid_interrupted, 'simulated stochastic interruption'),
    identical(hybrid_partial$phase, 'stochastic'),
    hybrid_partial$stochastic_score_evaluations >= 1L,
    is.finite(hybrid_partial$best_criterion)
)
hybrid_resumed <- cdrgam.fit(
    design,
    family='gaulss',
    backend='sparse',
    checkpoint=hybrid_checkpoint,
    sparse_control=list(
        gradient='hybrid',
        gradient_probes=4L,
        optimizer_maxit=1L,
        optimizer_gradient_tolerance=1e6
    )
)
stopifnot(
    inherits(hybrid_resumed, 'cdrgam_distributional_sparse'),
    identical(readRDS(hybrid_checkpoint)$stage, 'complete')
)

unlink(c(checkpoint, hybrid_checkpoint))
