library(cdrgam)

directions <- cdrgam:::.sparse_negative_curvature_directions
names <- c('a', 'b', 'c', 'log_scale')

# Negative curvature along one parameter is "single"; along a combination of
# parameters it is "mixed".
single <- directions(diag(c(2, -0.15, 3, 1)), 1:4, names, c(1, -4, 2))
rotation <- matrix(c(1, 1, -1, 1), 2) / sqrt(2)
ridge <- diag(4)
ridge[1:2, 1:2] <- rotation %*% diag(c(2, -0.5)) %*% t(rotation)
mixed <- directions(ridge, 1:4, names, c(1, 2, 3))
stopifnot(
    nrow(single) == 1L,
    identical(single$kind, 'single'),
    identical(single$parameter, 'b'),
    single$log_sp == -4,
    nrow(mixed) == 1L,
    identical(mixed$kind, 'mixed'),
    nrow(directions(diag(4), 1:4, names, c(1, 2, 3))) == 0L,
    grepl('stopped before the optimum',
        cdrgam:::.sparse_indefinite_message(single, 'REML'), fixed=TRUE),
    grepl('not positive definite',
        cdrgam:::.sparse_indefinite_message(mixed, 'REML'), fixed=TRUE)
)

# fit_diagnostics() reports stopped-short directions and large predicted
# decreases as optimizer_stopped_short, and keeps indefinite_outer_hessian
# for mixed directions.
fake_fit <- function(...) structure(list(
    sparse=list(
        convergence=c(list(converged=TRUE, code=0L, message='converged'),
            list(...)),
        condition_indicator=1
    ),
    cdrgam=list()
), class=c('cdrgam_sparse', 'cdrgam'))
codes <- function(fit) vapply(fit_diagnostics(fit)$reasons, `[[`, '', 'code')
stopifnot(
    identical(
        codes(fake_fit(hessian_positive_definite=FALSE,
            hessian_negative_directions=single)),
        'optimizer_stopped_short'
    ),
    identical(
        codes(fake_fit(hessian_positive_definite=FALSE,
            hessian_negative_directions=mixed)),
        'indefinite_outer_hessian'
    ),
    identical(
        codes(fake_fit(hessian_positive_definite=FALSE,
            hessian_negative_directions=rbind(single, mixed))),
        c('optimizer_stopped_short', 'indefinite_outer_hessian')
    ),
    identical(
        codes(fake_fit(hessian_positive_definite=TRUE,
            predicted_decrease=3.7)),
        'optimizer_stopped_short'
    ),
    identical(
        codes(fake_fit(hessian_positive_definite=FALSE,
            hessian_negative_directions=single, predicted_decrease=3.7)),
        'optimizer_stopped_short'
    ),
    length(codes(fake_fit(hessian_positive_definite=TRUE,
        predicted_decrease=0.01))) == 0L
)
