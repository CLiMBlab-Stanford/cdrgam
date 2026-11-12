library(cdrgam)

set.seed(42)
subject <- factor(rep(letters[1:4], each=20))
responses <- data.frame(
    subject=subject,
    time=rep(seq(0, 19), times=4),
    x=runif(80)
)
responses$y <- 0.4 * responses$x + rnorm(80, sd=0.2)
impulses <- responses[c('subject', 'time')]

formula <- y ~ s(x, subject, bs='fs', k=5) + irf(1)
design <- prepare_cdrgam(
    formula,
    impulses=impulses,
    responses=responses,
    series='subject',
    window=c(0, 4)
)
assembly <- cdrgam:::.fit_sparse_gaussian(
    design,
    family=gaussian(),
    method='GCV.Cp',
    setup_only=TRUE
)
native <- mgcv::gam(formula=y ~ s(x, subject, bs='fs', k=5),
    data=responses, method='GCV.Cp', fit=FALSE)

assembled_matrix <- do.call(cbind, assembly$matrices[1:2])
stopifnot(inherits(assembly$matrices[[2L]], 'sparseMatrix'))
stopifnot(ncol(assembly$setup$X) == 1L)
stopifnot(ncol(assembled_matrix) == ncol(native$X))
stopifnot(max(abs(as.matrix(assembled_matrix) - native$X)) < 1e-10)
stopifnot(identical(
    utils::head(assembly$sp_names, length(native$sp)),
    names(native$sp)
))

native_penalties <- lapply(seq_along(native$S), function(index) {
    output <- matrix(0, ncol(native$X), ncol(native$X))
    columns <- native$off[[index]] + seq_len(nrow(native$S[[index]])) - 1L
    output[columns, columns] <- native$S[[index]]
    output
})
for (index in seq_along(native_penalties)) {
    stopifnot(max(abs(
        as.matrix(assembly$penalty_components[[index]])[
            seq_len(ncol(native$X)),
            seq_len(ncol(native$X)),
            drop=FALSE
        ] -
            native_penalties[[index]]
    )) < 1e-10)
}

predicted <- cdrgam:::.cdr_predict_random_effect(
    assembly$random_prediction[[1L]],
    responses
)
stopifnot(max(abs(
    as.matrix(predicted) - as.matrix(assembly$matrices[[2L]])
)) < 1e-10)

new_responses <- responses[1:2, ]
new_responses$subject <- factor(c('a', 'new'))
predicted_new <- suppressWarnings(cdrgam:::.cdr_predict_random_effect(
    assembly$random_prediction[[1L]],
    new_responses
))
stopifnot(any(predicted_new[1L, ] != 0))
stopifnot(all(predicted_new[2L, ] == 0))

fit <- cdrgam(
    formula,
    impulses=impulses,
    responses=responses,
    series='subject',
    window=c(0, 4),
    family=gaussian(),
    method='GCV.Cp',
    backend='sparse',
    sparse_control=list(
        gradient='finite',
        outer_optimizer='lbfgsb',
        optimizer_maxit=5L,
        schur='always'
    )
)
training_matrix <- predict(
    fit,
    newdata=list(impulses=impulses, responses=responses),
    type='lpmatrix'
)
stopifnot(ncol(training_matrix) == length(coef(fit)))
stopifnot(max(abs(
    as.numeric(training_matrix %*% coef(fit)) - fit$linear.predictors
)) < 1e-8)
stopifnot(nrow(summary(fit)$s.table) == 2L)

cat('sparse factor smooth: ok\n')
