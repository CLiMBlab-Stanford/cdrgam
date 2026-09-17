library(cdrgam)
library(mgcv)

set.seed(20260916)
n <- 240
width <- 8
k <- 10
t_delta <- matrix(runif(n * width, 0, 2), n, width)
by <- matrix(rnorm(n * width), n, width)
mask <- matrix(runif(n * width) > 0.15, n, width)
t_delta[!mask] <- 0
by[!mask] <- 0
knots <- seq(0, 2, length.out=k)

compressed <- compress_cdr_smooth(
    t_delta,
    by,
    k=k,
    knots=knots,
    chunk_size=37
)

reference <- smoothCon(
    s(t_delta, k=k, bs='cr', by=by),
    data=list(t_delta=t_delta, by=by),
    knots=list(t_delta=knots),
    absorb.cons=TRUE,
    scale.penalty=TRUE,
    n=n
)[[1]]

stopifnot(isTRUE(all.equal(compressed$X, reference$X, tolerance=1e-11)))
stopifnot(isTRUE(all.equal(compressed$S, reference$S, tolerance=1e-11)))

beta <- seq(-0.5, 0.5, length.out=ncol(reference$X))
y <- drop(100 + reference$X %*% beta + rnorm(n, sd=0.25))
reference_fit <- gam(
    y ~ s(t_delta, k=k, bs='cr', by=by),
    data=list(y=y, t_delta=t_delta, by=by),
    knots=list(t_delta=knots),
    method='REML'
)
compressed_fit <- fit_compressed_cdr_gam(
    y,
    compressed,
    engine='gam',
    method='REML'
)

stopifnot(is_cdrgam(compressed_fit))
stopifnot(inherits(compressed_fit, 'gam'))
stopifnot(!is_cdrgam(as_gam(compressed_fit)))
stopifnot(inherits(as_gam(compressed_fit), 'gam'))
stopifnot(identical(compressed_fit$cdrgam$term_labels, 'irf_1'))

stopifnot(isTRUE(all.equal(
    fitted(compressed_fit),
    fitted(reference_fit),
    tolerance=1e-7
)))
stopifnot(isTRUE(all.equal(
    unname(compressed_fit$sp),
    unname(reference_fit$sp),
    tolerance=1e-6
)))
# Penalty scaling metadata is part of mgcv's variance-component definition,
# not merely fitting bookkeeping.
invisible(capture.output(reference_vcomp <- gam.vcomp(reference_fit)))
invisible(capture.output(compressed_vcomp <-
    variance_components(compressed_fit)))
stopifnot(isTRUE(all.equal(
    unname(compressed_vcomp),
    unname(reference_vcomp),
    tolerance=1e-6
)))

# Compression is invariant to response chunking.
single_chunk <- compress_cdr_smooth(
    t_delta,
    by,
    k=k,
    knots=knots,
    chunk_size=n
)
one_row_chunks <- compress_cdr_smooth(
    t_delta,
    by,
    k=k,
    knots=knots,
    chunk_size=1
)
stopifnot(isTRUE(all.equal(compressed$X, single_chunk$X, tolerance=0)))
stopifnot(isTRUE(all.equal(compressed$X, one_row_chunks$X, tolerance=0)))
stopifnot(isTRUE(all.equal(compressed$S, single_chunk$S, tolerance=0)))

# Native prediction, summary, covariance, and serialization continue to work.
compressed_prediction <- predict(
    compressed_fit,
    newdata=list(cdr_term_1=compressed$X)
)
stopifnot(isTRUE(all.equal(
    as.numeric(compressed_prediction),
    as.numeric(fitted(compressed_fit)),
    tolerance=1e-9
)))
stopifnot(inherits(summary(compressed_fit), 'summary.gam'))
stopifnot(all(is.finite(vcov(compressed_fit))))
model_path <- tempfile(fileext='.rds')
saveRDS(compressed_fit, model_path)
restored_fit <- readRDS(model_path)
stopifnot(is_cdrgam(restored_fit))
unlink(model_path)
restored_prediction <- predict(
    restored_fit,
    newdata=list(cdr_term_1=compressed$X)
)
stopifnot(isTRUE(all.equal(
    as.numeric(restored_prediction),
    as.numeric(compressed_prediction),
    tolerance=0
)))

# The BAM path uses the same compressed design and is deterministic.
reference_bam <- bam(
    y ~ s(t_delta, k=k, bs='cr', by=by),
    data=list(y=y, t_delta=t_delta, by=by),
    knots=list(t_delta=knots),
    method='fREML',
    chunk.size=53
)
compressed_bam <- fit_compressed_cdr_gam(
    y,
    compressed,
    engine='bam',
    method='fREML',
    chunk.size=53
)
compressed_bam_2 <- fit_compressed_cdr_gam(
    y,
    compressed,
    engine='bam',
    method='fREML',
    chunk.size=53
)
stopifnot(isTRUE(all.equal(
    fitted(compressed_bam),
    fitted(reference_bam),
    tolerance=1e-7
)))
stopifnot(isTRUE(all.equal(
    coef(compressed_bam),
    coef(compressed_bam_2),
    tolerance=0
)))
stopifnot(isTRUE(all.equal(
    compressed_bam$sp,
    compressed_bam_2$sp,
    tolerance=0
)))

# Multiple compressed IRFs retain separate penalties and smoothing parameters.
t_delta_2 <- matrix(runif(n * width, 0, 2), n, width)
by_2 <- matrix(rnorm(n * width), n, width)
t_delta_2[!mask] <- 0
by_2[!mask] <- 0
compressed_2 <- compress_cdr_smooth(
    t_delta_2,
    by_2,
    k=k,
    knots=knots,
    chunk_size=41
)
y_2 <- drop(
    50 +
    compressed$X %*% seq(-0.2, 0.3, length.out=k) +
    compressed_2$X %*% seq(0.25, -0.1, length.out=k) +
    rnorm(n, sd=0.3)
)
reference_multiple <- gam(
    y_2 ~ s(t_delta, k=k, bs='cr', by=by) +
        s(t_delta_2, k=k, bs='cr', by=by_2),
    data=list(
        y_2=y_2,
        t_delta=t_delta,
        by=by,
        t_delta_2=t_delta_2,
        by_2=by_2
    ),
    knots=list(t_delta=knots, t_delta_2=knots),
    method='REML'
)
compressed_multiple <- fit_cdrgam(
    y_2,
    list(primary=compressed, secondary=compressed_2),
    engine='gam',
    method='REML'
)
stopifnot(is_cdrgam(compressed_multiple))
stopifnot(identical(
    compressed_multiple$cdrgam$term_labels,
    c('primary', 'secondary')
))
stopifnot(length(compressed_multiple$sp) == 2)
stopifnot(isTRUE(all.equal(
    fitted(compressed_multiple),
    fitted(reference_multiple),
    tolerance=1e-7
)))
stopifnot(isTRUE(all.equal(
    unname(compressed_multiple$sp),
    unname(reference_multiple$sp),
    tolerance=1e-6
)))

# A non-Gaussian family uses the same design and penalty semantics.
poisson_mean <- exp(0.2 + drop(compressed$X %*%
    seq(-0.04, 0.04, length.out=k)))
poisson_y <- rpois(n, poisson_mean)
reference_poisson <- gam(
    poisson_y ~ s(t_delta, k=k, bs='cr', by=by),
    data=list(poisson_y=poisson_y, t_delta=t_delta, by=by),
    knots=list(t_delta=knots),
    family=poisson(),
    method='REML'
)
compressed_poisson <- fit_compressed_cdr_gam(
    poisson_y,
    compressed,
    family=poisson(),
    engine='gam',
    method='REML'
)
stopifnot(isTRUE(all.equal(
    fitted(compressed_poisson),
    fitted(reference_poisson),
    tolerance=1e-7
)))

# Constant linear-functional weight sums receive the same centering constraint
# and reduced penalty as native mgcv.
constant_by <- matrix(1, n, width)
constant_compressed <- compress_cdr_smooth(
    t_delta,
    constant_by,
    k=k,
    knots=knots,
    chunk_size=31
)
constant_reference <- smoothCon(
    s(t_delta, k=k, bs='cr', by=constant_by),
    data=list(t_delta=t_delta, constant_by=constant_by),
    knots=list(t_delta=knots),
    absorb.cons=TRUE,
    scale.penalty=TRUE,
    n=n
)[[1L]]
stopifnot(isTRUE(all.equal(
    constant_compressed$X,
    constant_reference$X,
    tolerance=1e-10
)))
stopifnot(isTRUE(all.equal(
    constant_compressed$S,
    constant_reference$S,
    tolerance=1e-10
)))
stopifnot(identical(
    constant_compressed$constraints,
    'linear-functional-centering'
))

# Unsupported or malformed inputs fail before fitting rather than changing the
# model silently.
expect_error <- function(expr, pattern) {
    message <- tryCatch(
        {
            force(expr)
            NA_character_
        },
        error=function(e) conditionMessage(e)
    )
    stopifnot(!is.na(message), grepl(pattern, message))
}
expect_error(
    compress_cdr_smooth(t_delta, by[, -1, drop=FALSE], k=k, knots=knots),
    'identical dimensions'
)
expect_error(
    compress_cdr_smooth(t_delta, by, k=k, bs='tp', knots=knots),
    'only bs="cr"'
)
