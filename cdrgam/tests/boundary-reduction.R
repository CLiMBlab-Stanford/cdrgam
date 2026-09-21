library(cdrgam)

practical <- cdrgam:::.analytic_outer_convergence_assessment(
    hessian=diag(c(1, 1e-9)),
    gradient=c(1e-3, 1e-6),
    criterion=590532,
    gradient_tolerance=2e-4,
    initial_radius=2
)
stopifnot(
    isTRUE(practical$converged),
    practical$diagnostics$predicted_improvement <
        practical$diagnostics$objective_floor,
    practical$diagnostics$flat_directions == 1L
)

noise_certified <- cdrgam:::.analytic_outer_convergence_assessment(
    hessian=diag(c(1, 1e-9)),
    gradient=c(0.02, 1e-6),
    criterion=3779611.693,
    gradient_tolerance=2e-4,
    initial_radius=2,
    objective_noise=0.005
)
stopifnot(
    isTRUE(noise_certified$converged),
    identical(noise_certified$diagnostics$observed_objective_noise, 0.005),
    noise_certified$diagnostics$objective_floor == 0.005,
    noise_certified$diagnostics$predicted_improvement < 0.005
)

noise_records <- list()
noise_assessments <- numeric()
noise_stagnation <- cdrgam:::.safeguarded_outer_bfgs(
    par=1,
    fn=function(x) {
        100 + x^2 / 2 + if (x < 1) 0.005 else 0
    },
    gr=function(x) x,
    lower=-2,
    upper=2,
    maxit=20L,
    initial_radius=1e-5,
    minimum_radius=1e-10,
    maximum_recovery_resets=0L,
    progress=function(record) {
        noise_records[[length(noise_records) + 1L]] <<- record
    },
    convergence_assessment=function(..., objective_noise=0) {
        noise_assessments <<- c(noise_assessments, objective_noise)
        list(
            converged=objective_noise > 0,
            message='test objective-noise certificate',
            diagnostics=list(observed_objective_noise=objective_noise)
        )
    }
)
stopifnot(
    noise_stagnation$convergence == 0L,
    length(noise_assessments) == 1L,
    noise_assessments[[1L]] >= 0.0049,
    any(vapply(noise_records, function(record) {
        identical(record$event, 'convergence_certified') &&
            record$assessment$observed_objective_noise >= 0.0049
    }, logical(1)))
)

noise_state <- list(
    version=1L,
    optimizer='safeguarded_outer_bfgs',
    parameters=1,
    criterion=10000.5,
    gradient=1,
    hessian=matrix(1e6),
    trust_radius=2,
    iteration=0L,
    function_evaluations=1L,
    gradient_evaluations=1L,
    accepted_steps=0L,
    rejected_steps=0L,
    consecutive_rejections=0L,
    curvature_resets=0L,
    consecutive_small_steps=0L,
    recovery_resets=0L
)
noise_fit <- cdrgam:::.safeguarded_outer_bfgs(
    par=1,
    fn=function(x) 10000 + x^2 / 2 + if (x < 1) 1.5e-6 else 0,
    gr=function(x) x,
    lower=-2,
    upper=2,
    maxit=1L,
    gradient_tolerance=1e-8,
    state=noise_state
)
stopifnot(
    noise_fit$history[[2L]]$accepted,
    grepl('noise_limited', noise_fit$history[[2L]]$step_type, fixed=TRUE)
)

# Two marginal curvature penalties leave a one-dimensional joint null space;
# the full-rank group shrinkage penalty survives in the reduced term.
term <- list(
    X=matrix(seq_len(40), nrow=10, ncol=4),
    S=list(
        diag(c(0, 0, 1, 1)),
        diag(c(0, 1, 0, 1)),
        diag(4)
    ),
    S.scale=c(1, 1, 1),
    transform=diag(4),
    rank=c(2L, 2L, 4L),
    null.space.dim=0L,
    constraints=character(),
    group='subject',
    group_levels=c('a', 'b'),
    base_dimension=4L,
    expanded_dimension=8L
)
class(term) <- c('cdrgam_term', 'cdr_compressed_term')
design <- list(
    terms=stats::setNames(list(term), 'x|subject'),
    simplifications=data.frame(
        term=character(), axis=character(), action=character(),
        requested=character(), effective=character(), reason=character(),
        stringsAsFactors=FALSE
    )
)
fit <- list(
    sp=stats::setNames(
        c(exp(14), exp(13), 0.2),
        c('s(cdr_term_1)1', 's(cdr_term_1)2', 's(cdr_term_1)3')
    ),
    optimizer=list(gradient=c(1e-7, -2e-7, 0.1)),
    sparse=list(gradient='exact')
)
class(fit) <- c('cdrgam_sparse', 'cdrgam')
plan <- cdrgam:::.cdr_boundary_reduction_plan(fit, design)
reduced <- cdrgam:::.cdr_apply_boundary_reduction(design, plan)
reduced_initial <- cdrgam:::.cdr_boundary_warm_start(fit, reduced)
untransformed <- design
untransformed$terms[[1L]]$transform <- NULL
untransformed_reduced <- cdrgam:::.cdr_apply_boundary_reduction(
    untransformed,
    plan
)
stopifnot(
    length(plan$entries) == 1L,
    identical(plan$entries[[1L]]$penalty_indices, 1:2),
    ncol(reduced$terms[[1L]]$X) == 1L,
    identical(dim(reduced$terms[[1L]]$transform), c(4L, 1L)),
    length(reduced$terms[[1L]]$S) == 1L,
    isTRUE(all.equal(reduced$terms[[1L]]$S[[1L]], matrix(1))),
    isTRUE(all.equal(
        unname(reduced_initial[['s(cdr_term_1)']]),
        log(0.2)
    )),
    reduced$terms[[1L]]$base_dimension == 1L,
    reduced$terms[[1L]]$expanded_dimension == 2L,
    isTRUE(all.equal(
        untransformed_reduced$terms[[1L]]$transform,
        plan$entries[[1L]]$basis
    )),
    nrow(reduced$simplifications) == 1L,
    reduced$simplifications$action == 'boundary_nullspace_reduction'
)

# A reduced term can retain an unpenalized null space but no smoothing
# parameters. Keep its coefficient block in the model without inventing a
# smoothing-parameter name or index.
zero_term <- term
zero_term$group <- NULL
zero_term$group_levels <- NULL
zero_term$base_dimension <- NULL
zero_term$expanded_dimension <- NULL
zero_term$S <- list(diag(c(0, 0, 1, 1)))
zero_term$S.scale <- 1
zero_term$rank <- 2L
zero_term$null.space.dim <- 2L
zero_design <- design
zero_design$terms <- stats::setNames(list(zero_term), 'x')
zero_plan <- list(
    entries=list(list(
        term_index=1L,
        penalty_indices=1L,
        basis=cbind(c(1, 0, 0, 0), c(0, 1, 0, 0)),
        original_dimension=4L,
        effective_dimension=2L,
        smoothing_parameters='s(cdr_term_1)',
        log_sp=14,
        score=1e-7
    )),
    table=data.frame(
        term='x', term_index=1L, penalties='s(cdr_term_1)',
        components='lag curvature', original_dimension=4L,
        effective_dimension=2L, maximum_log_sp=14,
        maximum_score=1e-7, status='certified_boundary_candidate',
        stringsAsFactors=FALSE
    )
)
zero_reduced <- cdrgam:::.cdr_apply_boundary_reduction(zero_design, zero_plan)
stopifnot(
    length(zero_reduced$terms[[1L]]$S) == 0L,
    length(cdrgam:::.sparse_sp_names('s(cdr_term_1)', 0L)) == 0L
)

# Exercise automatic reduction and reoptimization through the public sparse
# interface. A deliberately permissive gradient tolerance stops the pilot at
# its initial point; the low boundary threshold then makes both marginal
# penalties eligible while the full-rank group penalty remains explicit.
set.seed(20260917)
impulses <- data.frame(
    series=rep(c('a', 'b'), each=30),
    time=rep(seq(0, 29), 2),
    x=rnorm(60),
    z=rnorm(60)
)
responses <- data.frame(
    series=rep(c('a', 'b'), each=25),
    time=rep(seq(4, 28), 2),
    subject=factor(rep(c('a', 'b'), each=25)),
    response=rnorm(50)
)
prepared <- prepare_cdrgam(
    response ~ irf(
        x,
        window=c(0, 4),
        k=c(4, 4),
        nonlinear=TRUE,
        group=subject
    ) - irf(1),
    impulses,
    responses,
    series='series',
    history='ragged',
    quiet=TRUE
)

# Exercise sparse setup when one reduced term has no remaining penalties and
# another term still supplies a smoothing parameter.
zero_fit_design <- prepare_cdrgam(
    response ~ irf(x, window=c(0, 4), k=4) +
        irf(z, window=c(0, 4), k=4) - irf(1),
    impulses,
    responses,
    series='series',
    history='ragged',
    quiet=TRUE
)
zero_fit_term <- zero_fit_design$terms[[1L]]
zero_fit_basis <- cdrgam:::.cdr_penalty_null_basis(zero_fit_term$S)
zero_fit_plan <- list(
    entries=list(list(
        term_index=1L,
        penalty_indices=1L,
        basis=zero_fit_basis,
        original_dimension=ncol(zero_fit_term$X),
        effective_dimension=ncol(zero_fit_basis),
        smoothing_parameters='s(cdr_term_1)',
        log_sp=14,
        score=1e-7
    )),
    table=data.frame(
        term=names(zero_fit_design$terms)[[1L]],
        term_index=1L,
        penalties='s(cdr_term_1)',
        components='lag curvature',
        original_dimension=ncol(zero_fit_term$X),
        effective_dimension=ncol(zero_fit_basis),
        maximum_log_sp=14,
        maximum_score=1e-7,
        status='certified_boundary_candidate',
        stringsAsFactors=FALSE
    )
)
zero_fit_design <- cdrgam:::.cdr_apply_boundary_reduction(
    zero_fit_design,
    zero_fit_plan
)
zero_penalty_fit <- cdrgam.fit(
    zero_fit_design,
    backend='sparse',
    sparse_control=list(
        gradient='exact',
        outer_optimizer='bfgs_trust',
        optimizer_maxit=2L,
        optimizer_gradient_tolerance=1e6,
        boundary_action='report',
        hessian='none'
    ),
    rank_action='minimum_norm'
)
stopifnot(
    length(zero_fit_design$terms[[1L]]$S) == 0L,
    length(zero_penalty_fit$sp) == 1L,
    all(is.finite(zero_penalty_fit$coefficients))
)
boundary_events <- list()
automatic <- cdrgam.fit(
    prepared,
    backend='sparse',
    solver_trace=function(record) {
        boundary_events[[length(boundary_events) + 1L]] <<- record
    },
    sparse_control=list(
        gradient='exact',
        outer_optimizer='bfgs_trust',
        optimizer_gradient_tolerance=1e6,
        boundary_action='reduce',
        boundary_log_sp=-1,
        hessian='none'
    ),
    rank_action='minimum_norm'
)
stopifnot(
    isTRUE(automatic$converged),
    isTRUE(automatic$cdrgam$conditional_on_boundary_reduction),
    sum(automatic$cdrgam$boundary_reductions$status ==
        'applied_and_reoptimized') == 1L,
    sum(automatic$cdrgam$boundary_reductions$status ==
        'full_term_boundary_requires_confirmation') >= 1L,
    all(c('pilot', 'reoptimized') %in%
        automatic$cdrgam$boundary_reductions$fit_stage),
    any(automatic$cdrgam$preparation$simplifications$action ==
        'boundary_nullspace_reduction'),
    automatic$cdrgam$terms[[1L]]$base_dimension <
        prepared$terms[[1L]]$base_dimension,
    all(c(
        'boundary reduction selected',
        'reduced model refit started',
        'reduced model refit complete'
    ) %in% vapply(boundary_events, `[[`, character(1), 'event')),
    any(vapply(boundary_events, function(record) {
        identical(record$event, 'optimizer initialization') &&
            identical(
                record$source,
                'mapped converged boundary-pilot smoothing parameters'
            )
    }, logical(1)))
)
prediction_data <- responses[names(responses) != 'response']
automatic_prediction <- predict(
    automatic,
    newdata=list(impulses=impulses, responses=prediction_data)
)
stopifnot(max(abs(automatic_prediction - fitted(automatic))) < 1e-7)
automatic_irf <- estimate_irf(automatic, term=1, n=7, n_predictor=5, se=FALSE)
stopifnot(
    nrow(automatic_irf) == 35L * nlevels(responses$subject),
    all(is.finite(automatic_irf$estimate))
)
