library(cdrgam)

# These checks compare designs with native-unit references, so they fit
# without predictor rescaling.

set.seed(7)
subject <- factor(rep(letters[1:5], each=40))
responses <- data.frame(
    subject=subject,
    time=rep(seq(0, 39), times=5),
    x=runif(200),
    z=runif(200)
)
subject_shift <- rnorm(5)
subject_slope <- rnorm(5)
responses$y <- subject_shift[subject] +
    subject_slope[subject] * sin(2 * pi * responses$x) +
    rnorm(200, sd=0.3)
impulses <- responses[c('subject', 'time')]

prepare <- function(formula) {
    prepare_cdrgam(rescale_predictors=FALSE,
        formula,
        impulses=impulses,
        responses=responses,
        series='subject',
        window=c(0, 4)
    )
}
assemble <- function(formula) {
    cdrgam:::.fit_sparse_gaussian(
        prepare(formula),
        family=gaussian(),
        method='GCV.Cp',
        setup_only=TRUE
    )
}

# Compare one sparse ti() term with mgcv's own construction. The sparse layout
# orders columns by factor level, so the designs agree up to a permutation.
check_term <- function(term) {
    formula <- stats::as.formula(paste('y ~', term, '+ irf(1)'))
    assembly <- assemble(formula)
    native <- mgcv::gam(
        stats::as.formula(paste('y ~', term)),
        data=responses, method='GCV.Cp', fit=FALSE
    )
    stopifnot(ncol(assembly$setup$X) == 1L)
    sparse_term <- assembly$matrices[[2L]]
    stopifnot(inherits(sparse_term, 'sparseMatrix'))
    native_term <- native$X[, -1L, drop=FALSE]
    stopifnot(ncol(sparse_term) == ncol(native_term))
    dense_term <- as.matrix(sparse_term)
    permutation <- vapply(seq_len(ncol(dense_term)), function(column) {
        distance <- colSums(abs(native_term - dense_term[, column]))
        which.min(distance)
    }, integer(1))
    stopifnot(!anyDuplicated(permutation))
    stopifnot(max(abs(native_term[, permutation] - dense_term)) < 1e-10)
    stopifnot(identical(
        utils::head(assembly$sp_names, length(native$sp)),
        names(native$sp)
    ))
    term_columns <- 1L + seq_len(ncol(dense_term))
    for (index in seq_along(native$S)) {
        native_penalty <- native$S[[index]][permutation, permutation]
        sparse_penalty <- as.matrix(
            assembly$penalty_components[[index]]
        )[term_columns, term_columns, drop=FALSE]
        stopifnot(max(abs(sparse_penalty - native_penalty)) < 1e-10)
    }
    predicted <- cdrgam:::.cdr_predict_random_effect(
        assembly$random_prediction[[1L]],
        responses
    )
    stopifnot(max(abs(as.matrix(predicted) - dense_term)) < 1e-10)
    invisible(assembly)
}

check_term("ti(x, subject, bs=c('cr', 're'), k=c(5, 5))")
check_term("ti(x, subject, bs='fs', k=5)")
check_term("ti(subject, x, bs=c('re', 'tp'), k=c(5, 6))")
check_term("ti(x, z, subject, bs=c('cr', 'cr', 're'), k=c(4, 4, 5))")

# Unseen levels receive zero subject-specific contributions.
assembly <- assemble(y ~ ti(x, subject, bs=c('cr', 're'), k=c(5, 5)) + irf(1))
new_responses <- responses[1:2, ]
new_responses$subject <- factor(c('a', 'new'))
predicted_new <- suppressWarnings(cdrgam:::.cdr_predict_random_effect(
    assembly$random_prediction[[1L]],
    new_responses
))
stopifnot(any(predicted_new[1L, ] != 0))
stopifnot(all(predicted_new[2L, ] == 0))

# Terms without exactly one random-effect factor margin stay in the dense
# ordinary formula.
responses$block <- factor(rep(c('u', 'v'), times=100))
impulses$block <- responses$block
dense <- assemble(y ~ ti(x, z, k=c(4, 4)) + irf(1))
stopifnot(ncol(dense$setup$X) > 1L, !length(dense$random_prediction))
dense <- assemble(y ~ ti(x, subject, block, bs=c('cr', 're', 're'),
    k=c(4, 5, 2)) + irf(1))
stopifnot(ncol(dense$setup$X) > 1L, !length(dense$random_prediction))

# The sparse fit agrees with the native mgcv backend on the same model.
formula <- y ~ s(x, k=5) + s(subject, bs='re') +
    ti(x, subject, bs=c('cr', 're'), k=c(5, 5)) + irf(1)
fit_sparse <- cdrgam(rescale_predictors=FALSE,
    formula,
    impulses=impulses,
    responses=responses,
    series='subject',
    window=c(0, 4),
    family=gaussian(),
    method='GCV.Cp',
    backend='sparse',
    sparse_control=list(optimizer_maxit=1000L)
)
fit_native <- cdrgam(rescale_predictors=FALSE,
    formula,
    impulses=impulses,
    responses=responses,
    series='subject',
    window=c(0, 4),
    family=gaussian(),
    method='GCV.Cp',
    backend='mgcv'
)
stopifnot(
    abs(fit_sparse$gcv.ubre - fit_native$gcv.ubre) < 1e-5 * fit_native$gcv.ubre,
    max(abs(fitted(fit_sparse) - fitted(fit_native))) < 1e-3
)
training_matrix <- predict(
    fit_sparse,
    newdata=list(impulses=impulses, responses=responses),
    type='lpmatrix'
)
stopifnot(ncol(training_matrix) == length(coef(fit_sparse)))
stopifnot(max(abs(
    as.numeric(training_matrix %*% coef(fit_sparse)) -
        fit_sparse$linear.predictors
)) < 1e-8)

cat('sparse tensor interaction: ok\n')
