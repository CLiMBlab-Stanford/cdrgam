# Randomized location-scale equivalence checks across all fitting backends.
# Run from the repository root with:
#   Rscript validation/location_scale_randomized.R

library(cdrgam)

simulate_location_scale_streams <- function(seed, subjects=4L) {
    set.seed(seed)
    impulse_streams <- vector('list', subjects)
    response_streams <- vector('list', subjects)
    subject_offsets <- seq(-18, 18, length.out=subjects)
    for (subject_index in seq_len(subjects)) {
        subject <- paste0('s', subject_index)
        impulse_time <- sort(stats::runif(55, 0, 18))
        response_time <- sort(stats::runif(75, 0.5, 18))
        x <- stats::rnorm(length(impulse_time))
        z <- stats::rnorm(length(impulse_time))
        location <- numeric(length(response_time))
        log_sd <- numeric(length(response_time))
        for (row in seq_along(response_time)) {
            lag <- response_time[[row]] - impulse_time
            linked <- lag >= 0 & lag <= 1.5
            location[[row]] <- 320 + subject_offsets[[subject_index]] +
                sum(7 * x[linked] * exp(-1.4 * lag[linked]))
            log_sd[[row]] <- log(38) +
                sum(0.06 * z[linked] * exp(-2 * lag[linked]))
        }
        impulse_streams[[subject_index]] <- data.frame(
            subject=subject, time=impulse_time, x=x, z=z
        )
        response_streams[[subject_index]] <- data.frame(
            subject=subject,
            time=response_time,
            response=stats::rnorm(
                length(response_time), location, exp(log_sd)
            )
        )
    }
    impulses <- do.call(rbind, impulse_streams)
    responses <- do.call(rbind, response_streams)
    levels <- paste0('s', seq_len(subjects))
    impulses$subject <- factor(impulses$subject, levels=levels)
    responses$subject <- factor(responses$subject, levels=levels)
    list(impulses=impulses, responses=responses)
}

fit_replication <- function(seed) {
    streams <- simulate_location_scale_streams(seed)
    formulas <- list(
        location=response ~ s(subject, bs='re') +
            irf(x, window=c(0, 1.5), k_l=6) - irf(1),
        scale=~ s(subject, bs='re') +
            irf(z, window=c(0, 1.5), k_l=6) - irf(1)
    )
    design <- prepare_cdrgam(
        formulas,
        streams$impulses,
        streams$responses,
        series='subject',
        history='ragged',
        chunk_size=61,
        quiet=TRUE
    )
    fits <- list(
        mgcv=cdrgam.fit(
            design, family='gaulss', backend='mgcv', engine='gam',
            method='REML'
        ),
        sparse=cdrgam.fit(
            design, family='gaulss', backend='sparse', method='REML',
            sparse_control=list(
                optimizer_maxit=500L,
                optimizer_gradient_tolerance=2e-3
            )
        )
    )
    newdata <- list(
        impulses=streams$impulses,
        responses=streams$responses[c('subject', 'time')]
    )
    predictions <- lapply(fits, predict, newdata=newdata, type='response')
    data.frame(
        seed=seed,
        converged=fits$sparse$converged,
        inner_gradient=fits$sparse$distributional$gradient_norm,
        outer_gradient=max(abs(fits$sparse$optimizer$gradient)),
        coefficient_max_difference=max(
            abs(coef(fits$sparse) - coef(fits$mgcv))
        ),
        prediction_max_difference=max(
            abs(predictions$sparse - predictions$mgcv)
        ),
        log_likelihood_difference=abs(
            as.numeric(logLik(fits$sparse)) - as.numeric(logLik(fits$mgcv))
        )
    )
}

seeds <- c(1401L, 1402L, 1403L)
metrics <- do.call(rbind, lapply(seeds, fit_replication))
print(metrics, row.names=FALSE)
stopifnot(
    all(metrics$converged),
    max(metrics$outer_gradient) < 2e-3,
    max(metrics$coefficient_max_difference) < 2e-2,
    max(metrics$prediction_max_difference) < 3e-2,
    max(metrics$log_likelihood_difference) < 1e-2
)
