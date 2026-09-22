# Full Brown rate-model comparison for outer smoothing-parameter optimizers.
# The mixed-effects design is full by subject and random-intercept-only by
# item (docid x docpos).

library(cdrgam)

config <- yaml::read_yaml(Sys.getenv('CDRGAM_BROWN_CONFIG', 'brown.yml'))
read_stream <- function(path) {
    utils::read.csv(path, sep=config$data$sep, header=TRUE)
}
apply_filters <- function(data) {
    keep <- rep.int(TRUE, nrow(data))
    for (filter in config$data$filters) {
        if (!is.null(filter$column)) {
            keep <- keep & do.call(
                match.fun(filter$fun),
                c(list(data[[filter$column]]), list(filter$args))
            )
        } else {
            eligible <- names(which(table(data[[filter$factor]][keep]) >=
                filter$min))
            keep <- keep & as.character(data[[filter$factor]]) %in% eligible
        }
    }
    data[keep, , drop=FALSE]
}
prepare_response <- function(data) {
    data$docpos <- data$Word_Number
    data$item <- interaction(data$docid, data$docpos, drop=TRUE, sep=':')
    data$subject <- droplevels(factor(data$subject))
    data$item <- droplevels(factor(data$item))
    data
}

impulses <- read_stream(config$data$X_train)
responses <- prepare_response(apply_filters(read_stream(config$data$Y_train)))
impulses$subject <- factor(impulses$subject, levels=levels(responses$subject))
window <- c(0, as.numeric(config$data$t_delta_cutoff))
lag_k <- as.integer(Sys.getenv('CDRGAM_BROWN_LAG_K', '6'))
formula <- stats::as.formula(sprintf(
    paste0(
        "fdur ~ s(subject, bs='re') + s(item, bs='re') + ",
        'irf(1, window=c(%s,%s), k=%d) + ',
        'irf(1, window=c(%s,%s), k=%d, group=subject)'
    ),
    window[[1L]], window[[2L]], lag_k,
    window[[1L]], window[[2L]], lag_k
))
design <- prepare_cdrgam(
    formula,
    impulses,
    responses,
    series=c('subject', 'docid'),
    history='auto',
    history_length=as.integer(config$data$history_length),
    chunk_size=10000L,
    quiet=FALSE
)
stopifnot(
    !any(vapply(
        design$terms,
        function(term) identical(term$group, 'item'),
        logical(1)
    )),
    any(vapply(
        design$terms,
        function(term) identical(term$group, 'subject'),
        logical(1)
    ))
)

specifications <- list(
    finite_lbfgsb=list(gradient='finite', outer_optimizer='lbfgsb'),
    exact_lbfgsb=list(gradient='exact', outer_optimizer='lbfgsb'),
    exact_bfgs_trust=list(gradient='exact', outer_optimizer='bfgs_trust'),
    finite_warm_bfgs=list(
        gradient='finite',
        outer_optimizer='lbfgsb',
        warm_start='exact_bfgs_trust'
    )
)
selected <- strsplit(Sys.getenv(
    'CDRGAM_BROWN_OPTIMIZERS',
    paste(names(specifications), collapse=',')
), ',', fixed=TRUE)[[1L]]
selected <- trimws(selected)
stopifnot(length(selected) > 0L, all(selected %in% names(specifications)))

rows <- list()
fits <- list()
for (name in selected) {
    specification <- specifications[[name]]
    message('Brown outer-optimizer benchmark: ', name)
    checkpoint <- NULL
    if (!is.null(specification$warm_start)) {
        warm_fits <- readRDS(Sys.getenv(
            'CDRGAM_BROWN_WARM_FITS',
            'validation/output/brown_outer_optimizer_benchmark_fits.rds'
        ))
        warm <- warm_fits[[specification$warm_start]]
        stopifnot(!is.null(warm), all(is.finite(warm$sp)), is.finite(warm$reml))
        checkpoint <- tempfile('cdrgam-warm-', fileext='.rds')
        saveRDS(list(log_sp=log(warm$sp), criterion=warm$reml), checkpoint)
    }
    elapsed <- system.time({
        fit <- cdrgam.fit(
            design,
            backend='sparse',
            method='REML',
            solver_trace=1,
            checkpoint=checkpoint,
            sparse_control=list(
                gradient=specification$gradient,
                outer_optimizer=specification$outer_optimizer,
                optimizer_maxit=150L,
                optimizer_gradient_tolerance=1e-4,
                hessian='profiled',
                schur='always',
                crossprod_chunk_size=5000L
            )
        )
    })[['elapsed']]
    fits[[name]] <- fit
    rows[[name]] <- data.frame(
        specification=name,
        gradient=specification$gradient,
        outer_optimizer=specification$outer_optimizer,
        observations=nrow(responses),
        coefficients=length(coef(fit)),
        smoothing_parameters=length(fit$sp),
        reml=fit$reml,
        elapsed_seconds=unname(elapsed),
        objective_evaluations=
            fit$sparse$convergence$total_objective_evaluations,
        hessian_evaluations=fit$sparse$hessian_evaluations,
        optimization_factorizations=fit$sparse$numeric_updates -
            fit$sparse$hessian_evaluations,
        total_factorizations=fit$sparse$numeric_updates,
        optimizer_code=fit$optimizer$convergence,
        optimizer_message=if (is.null(fit$optimizer$message)) '' else
            fit$optimizer$message,
        optimizer_iterations=if (is.null(fit$optimizer$iterations)) {
            NA_integer_
        } else {
            fit$optimizer$iterations
        },
        gradient_max=if (is.null(fit$optimizer$gradient)) NA_real_ else
            max(abs(fit$optimizer$gradient)),
        converged=fit$converged
    )
}
results <- do.call(rbind, rows)
output <- Sys.getenv(
    'CDRGAM_BROWN_OPTIMIZER_OUTPUT',
    'validation/output/brown_outer_optimizer_benchmark.csv'
)
utils::write.csv(results, output, row.names=FALSE)
saveRDS(
    lapply(fits, function(fit) list(
        reml=fit$reml,
        sp=fit$sp,
        coefficients=coef(fit),
        optimizer=fit$optimizer,
        convergence=fit$sparse$convergence
    )),
    sub('\\.csv$', '_fits.rds', output),
    compress=FALSE
)
print(results, row.names=FALSE)
