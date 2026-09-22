library(cdrgam)

stopifnot(
    identical(cdrgam_family('poisson')$link, 'log'),
    identical(cdrgam_family('binomial', 'probit')$link, 'probit'),
    identical(cdrgam_family(stats::Gamma(link='log'))$family, 'Gamma')
)
stopifnot(inherits(
    try(cdrgam_family('negative-binomial'), silent=TRUE),
    'try-error'
))

simulation <- simulate_cdr(
    list(signal=function(lag) exp(-2 * lag)),
    n_impulses=180,
    n_responses=160,
    duration=30,
    window=1.5,
    seed=9021
)
latent <- as.numeric(scale(simulation$responses$response))
set.seed(9020)
simulation$impulses$noise <- rnorm(nrow(simulation$impulses))
streams <- function(response) {
    responses <- simulation$responses
    responses$response <- response
    list(impulses=simulation$impulses, responses=responses)
}
fit_and_check <- function(response, family) {
    data <- streams(response)
    data$responses$base_offset <- seq(-0.1, 0.1, length.out=nrow(data$responses))
    prior_weights <- rep(c(1, 2), length.out=nrow(data$responses))
    design <- prepare_cdrgam(
        response ~ offset(base_offset) +
            irf(signal, window=c(0, 1.5), k_l=6) +
            irf(noise, window=c(0, 1.5), k_l=5) - irf(1),
        data$impulses,
        data$responses,
        history='ragged',
        quiet=TRUE
    )
    fit <- suppressWarnings(cdrgam.fit(
        design,
        family=family,
        backend='mgcv',
        engine='gam',
        method='REML',
        weights=prior_weights
    ))
    setup <- getFromNamespace('.fit_compressed_mgcv', 'cdrgam')(
        y=design$responses[[design$response_name]],
        terms=design$terms,
        family=family,
        method='REML',
        engine='gam',
        base_formula=design$ordinary_formula,
        response_data=design$responses,
        setup_only=TRUE,
        drop.unused.levels=design$configuration$drop.unused.levels,
        weights=prior_weights
    )
    pirls <- getFromNamespace('.cdrgam_dense_pirls', 'cdrgam')(
        setup,
        family,
        fit$sp,
        tolerance=1e-10
    )
    sparse_pirls <- getFromNamespace('.cdrgam_sparse_pirls', 'cdrgam')(
        setup,
        family,
        fit$sp,
        tolerance=1e-10
    )
    assembly <- getFromNamespace('.fit_sparse_gaussian', 'cdrgam')(
        design,
        family,
        method='REML',
        sparse_control=list(crossprod_chunk_size=23L),
        setup_only=TRUE,
        weights=prior_weights
    )
    streamed_pirls <- getFromNamespace(
        '.cdrgam_streamed_sparse_pirls',
        'cdrgam'
    )(
        assembly,
        family,
        fit$sp,
        tolerance=1e-10,
        chunk_size=23L
    )
    canonical <- (identical(family$family, 'binomial') &&
            identical(family$link, 'logit')) ||
        (identical(family$family, 'poisson') &&
            identical(family$link, 'log'))
    estimated_gamma <- identical(family$family, 'Gamma') &&
        identical(family$link, 'log')
    custom_supported <- canonical || estimated_gamma
    laml <- if (custom_supported) {
        getFromNamespace('.cdrgam_dense_laml', 'cdrgam')(
            setup,
            family,
            pirls_tolerance=1e-10
        )
    } else NULL
    sparse_laml <- if (custom_supported) {
        evaluator <- getFromNamespace(
            '.cdrgam_streamed_sparse_laml',
            'cdrgam'
        )
        value <- evaluator(
            assembly,
            family,
            c(log(fit$sp), if (estimated_gamma) log(fit$scale) else NULL),
            tolerance=1e-10,
            chunk_size=23L,
            score=TRUE
        )
        step <- 1e-4
        fitted_parameters <- c(
            log(fit$sp),
            if (estimated_gamma) log(fit$scale) else NULL
        )
        finite_score <- vapply(seq_along(fitted_parameters), function(index) {
            upper <- lower <- fitted_parameters
            upper[[index]] <- upper[[index]] + step
            lower[[index]] <- lower[[index]] - step
            (evaluator(
                assembly, family, upper,
                tolerance=1e-10, chunk_size=23L
            )$criterion - evaluator(
                assembly, family, lower,
                tolerance=1e-10, chunk_size=23L
            )$criterion) / (2 * step)
        }, numeric(1))
        stopifnot(max(abs(value$score - finite_score)) < 2e-4)
        value
    } else NULL
    sparse_optimized <- if (custom_supported) {
        getFromNamespace(
            '.cdrgam_optimize_streamed_sparse_laml',
            'cdrgam'
        )(
            assembly,
            family,
            tolerance=1e-10,
            chunk_size=23L,
            max_iterations=60L,
            gradient_tolerance=if (estimated_gamma) 5e-5 else 1e-5
        )
    } else NULL
    newdata <- list(
        impulses=data$impulses,
        responses=data$responses[c('time', 'base_offset')]
    )
    response_prediction <- predict(fit, newdata, type='response')
    link_prediction <- predict(fit, newdata, type='link')
    response_se <- predict(fit, newdata, type='response', se.fit=TRUE)
    stopifnot(
        inherits(fit, 'cdrgam'),
        identical(fit$family$family, family$family),
        identical(fit$family$link, family$link),
        isTRUE(pirls$converged),
        isTRUE(sparse_pirls$converged),
        isTRUE(streamed_pirls$converged),
        sparse_pirls$numeric_updates > 0L,
        streamed_pirls$numeric_updates > 0L,
        streamed_pirls$chunks > 1L,
        max(abs(pirls$coefficients - coef(fit))) < 2e-5,
        max(abs(sparse_pirls$coefficients - pirls$coefficients)) < 1e-8,
        max(abs(
            streamed_pirls$linear_predictors - pirls$linear_predictors
        )) < 1e-8,
        max(abs(
            sparse_pirls$linear_predictors - pirls$linear_predictors
        )) < 1e-8,
        max(abs(pirls$linear_predictors - fit$linear.predictors)) < 2e-5,
        abs(pirls$deviance - deviance(fit)) < 2e-5,
        max(abs(response_prediction - fitted(fit))) < 1e-8,
        max(abs(response_prediction - family$linkinv(link_prediction))) < 1e-8,
        all(is.finite(response_se$se.fit)),
        all(response_se$se.fit >= 0),
        identical(summary(fit)$family$family, family$family)
    )
    if (!is.null(laml)) {
        if (!estimated_gamma) {
            stopifnot(abs(
                sparse_laml$criterion - laml$solution$criterion
            ) < 1e-3)
        }
        stopifnot(
            sparse_optimized$optimization$convergence == 0L,
            max(abs(
                sparse_optimized$retained$solution$linear_predictors -
                fit$linear.predictors
            )) < 3e-3
        )
        if (estimated_gamma) {
            stopifnot(
                abs(log(laml$solution$scale / fit$scale)) < 0.02,
                abs(log(sparse_optimized$retained$scale / fit$scale)) < 0.02
            )
        }
        stopifnot(
            isTRUE(laml$solution$converged),
            max(abs(
                laml$solution$linear_predictors - fit$linear.predictors
            )) < 3e-3,
            abs(laml$solution$deviance - deviance(fit)) < 2e-2
        )
        if (estimated_gamma) {
            stopifnot(max(abs(
                laml$solution$linear_predictors - fit$linear.predictors
            )) < 3e-3)
        } else if (max(fit$sp) < 1e6) {
            stopifnot(max(abs(log(laml$solution$sp / fit$sp))) < 0.01)
        } else {
            stopifnot(min(laml$solution$sp) > 1e6)
        }
        block <- suppressWarnings(cdrgam.fit(
            design,
            family=family,
            backend='block',
            method='REML',
            weights=prior_weights
        ))
        stopifnot(
            isTRUE(block$converged),
            max(abs(block$linear.predictors - fit$linear.predictors)) < 3e-3,
            if (estimated_gamma) {
                max(abs(fitted(block) / fitted(fit) - 1)) < 3e-3
            } else max(abs(fitted(block) - fitted(fit))) < 1e-3,
            abs(deviance(block) - deviance(fit)) < 2e-2,
            abs(as.numeric(logLik(block)) - as.numeric(logLik(fit))) <
                if (estimated_gamma) 0.2 else 1e-2,
            max(abs(predict(block, newdata, type='link') -
                fit$linear.predictors)) < 3e-3
        )
        if (estimated_gamma) {
            stopifnot(abs(log(block$scale / fit$scale)) < 0.02)
        }
        block_summary <- summary(block)
        stopifnot(
            inherits(block_summary, 'summary.cdrgam_block'),
            if (estimated_gamma) {
                't value' %in% colnames(block_summary$p.table)
            } else 'z value' %in% colnames(block_summary$p.table),
            if (estimated_gamma) {
                'Pr(>|t|)' %in% colnames(block_summary$p.table)
            } else 'Pr(>|z|)' %in% colnames(block_summary$p.table),
            if (estimated_gamma) {
                'F' %in% colnames(block_summary$s.table)
            } else 'Chi.sq' %in% colnames(block_summary$s.table),
            abs(block_summary$dev.expl - summary(fit)$dev.expl) < 1e-3
        )
        sparse_fit <- suppressWarnings(cdrgam.fit(
            design,
            family=family,
            backend='sparse',
            method='REML',
            weights=prior_weights,
            sparse_control=list(
                crossprod_chunk_size=23L,
                optimizer_maxit=60L,
                optimizer_gradient_tolerance=
                    if (estimated_gamma) 5e-5 else 1e-5
            )
        ))
        sparse_summary <- summary(sparse_fit)
        sparse_covariance <- vcov(sparse_fit)
        sparse_irf <- estimate_irf(sparse_fit, n=11L)
        sparse_unconditional <- try(
            vcov(sparse_fit, unconditional=TRUE),
            silent=TRUE
        )
        stopifnot(
            isTRUE(sparse_fit$converged),
            max(abs(
                sparse_fit$linear.predictors - fit$linear.predictors
            )) < 3e-3,
            if (estimated_gamma) {
                max(abs(fitted(sparse_fit) / fitted(fit) - 1)) < 3e-3
            } else max(abs(fitted(sparse_fit) - fitted(fit))) < 1e-3,
            abs(deviance(sparse_fit) - deviance(fit)) < 2e-2,
            abs(as.numeric(logLik(sparse_fit)) -
                as.numeric(logLik(fit))) <
                if (estimated_gamma) 0.2 else 1e-2,
            max(abs(predict(sparse_fit, newdata, type='link') -
                fit$linear.predictors)) < 3e-3,
            if (estimated_gamma) {
                't value' %in% colnames(sparse_summary$p.table)
            } else 'z value' %in% colnames(sparse_summary$p.table),
            if (estimated_gamma) {
                'F' %in% colnames(sparse_summary$s.table)
            } else 'Chi.sq' %in% colnames(sparse_summary$s.table),
            all(is.finite(diag(sparse_covariance))),
            all(diag(sparse_covariance) >= 0),
            all(is.finite(sparse_irf$estimate)),
            all(is.finite(sparse_irf$se)),
            inherits(sparse_unconditional, 'try-error'),
            identical(
                sparse_fit$sparse$control$outer_optimizer,
                if (estimated_gamma) 'lbfgsb' else 'bfgs_trust'
            )
        )
        if (estimated_gamma) {
            stopifnot(abs(log(sparse_fit$scale / fit$scale)) < 0.02)
        }
    }
    for (backend in if (custom_supported) character() else 'sparse') {
        error <- try(cdrgam.fit(design, family=family, backend=backend), silent=TRUE)
        stopifnot(
            inherits(error, 'try-error'),
            grepl('currently supports', as.character(error), fixed=TRUE)
        )
    }
    invisible(fit)
}

set.seed(9022)
poisson <- cdrgam_family('poisson', 'log')
fit_and_check(rpois(length(latent), exp(0.7 + 0.25 * latent)), poisson)

set.seed(9023)
binomial <- cdrgam_family('binomial', 'logit')
fit_and_check(rbinom(length(latent), 1, plogis(-0.2 + 0.6 * latent)), binomial)

set.seed(9024)
gamma <- cdrgam_family('Gamma', 'log')
mean <- exp(1 + 0.2 * latent)
fit_and_check(rgamma(length(latent), shape=5, scale=mean / 5), gamma)

set.seed(9025)
binomial_probit <- cdrgam_family('binomial', 'probit')
fit_and_check(
    rbinom(length(latent), 1, stats::pnorm(-0.2 + 0.6 * latent)),
    binomial_probit
)

set.seed(9026)
random_data <- streams(rpois(length(latent), exp(0.7 + 0.25 * latent)))
random_data$responses$subject <- factor(rep(
    paste0('s', 1:8),
    length.out=nrow(random_data$responses)
))
random_design <- prepare_cdrgam(
    response ~ s(subject, bs='re') +
        irf(signal, window=c(0, 1.5), k_l=6) - irf(1),
    random_data$impulses,
    random_data$responses,
    history='ragged',
    quiet=TRUE
)
random_native <- cdrgam.fit(
    random_design,
    family=poisson,
    backend='mgcv',
    engine='gam',
    method='REML'
)
random_sparse <- cdrgam.fit(
    random_design,
    family=poisson,
    backend='sparse',
    method='REML',
    sparse_control=list(
        crossprod_chunk_size=19L,
        optimizer_maxit=60L,
        optimizer_gradient_tolerance=1e-5
    )
)
stopifnot(
    isTRUE(random_sparse$converged),
    max(abs(
        random_sparse$linear.predictors - random_native$linear.predictors
    )) < 3e-3,
    abs(deviance(random_sparse) - deviance(random_native)) < 2e-2
)
