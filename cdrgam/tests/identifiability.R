library(cdrgam)

# Every response has the same three constant-predictor delays, but the
# non-trivial predictor histories vary. The implicit deconvolutional intercept
# therefore collapses after centering while the requested IRF remains usable.
impulses <- data.frame(
    time=0:11,
    x=c(-1, 0.5, 1.2, -0.7, 0.3, 1.5, -1.1, 0.8, 0.2, -0.4, 1, -0.2)
)
responses <- data.frame(
    time=2:11,
    response=seq(1, 2, length.out=10)
)

implicit_warning <- NULL
implicit <- withCallingHandlers(
    prepare_cdrgam(
        response ~ irf(x, window=c(0, 2), k=3),
        impulses,
        responses,
        history='ragged',
        quiet=TRUE
    ),
    warning=function(w) {
        implicit_warning <<- conditionMessage(w)
        invokeRestart('muffleWarning')
    }
)
stopifnot(grepl('implicit irf(1)', implicit_warning, fixed=TRUE))
stopifnot(length(implicit$terms) == 1L)
stopifnot(identical(names(implicit$terms), 'x'))
stopifnot(identical(implicit$identifiability$rate$requested, 'implicit'))
stopifnot(identical(implicit$identifiability$rate$status, 'removed'))
stopifnot(identical(
    implicit$identifiability$rate$reason,
    'zero design after centering'
))
stopifnot(grepl(
    'irf\\(1',
    paste(deparse(cdr_formula(implicit, type='normalized')), collapse='')
))
stopifnot(!grepl(
    'irf\\(1',
    paste(deparse(cdr_formula(implicit, type='effective')), collapse='')
))

suppressed <- prepare_cdrgam(
    response ~ irf(x, window=c(0, 2), k=3) - irf(1),
    impulses,
    responses,
    history='ragged',
    quiet=TRUE
)
stopifnot(length(suppressed$terms) == 1L)
stopifnot(identical(suppressed$identifiability$rate$requested, 'suppressed'))
stopifnot(identical(suppressed$identifiability$rate$status, 'suppressed'))
stopifnot(isTRUE(all.equal(
    implicit$terms[[1L]]$X,
    suppressed$terms[[1L]]$X,
    tolerance=0
)))

explicit_error <- tryCatch(
    {
        prepare_cdrgam(
            response ~ irf(1, window=c(0, 2), k=3) +
                irf(x, window=c(0, 2), k=3),
            impulses,
            responses,
            history='ragged',
            quiet=TRUE
        )
        NA_character_
    },
    error=function(e) conditionMessage(e)
)
stopifnot(grepl('Explicit irf(1) is unidentifiable', explicit_error, fixed=TRUE))

# With irregular histories, the implicit rate is retained and prediction does
# not require a physical column named "rate" or "1" in the impulse stream.
irregular_responses <- data.frame(
    time=c(1.2, 2.1, 3.8, 5.1, 7.7, 10.4),
    response=seq_len(6)
)
irregular <- prepare_cdrgam(
    response ~ irf(x, window=c(0, 2), k=3),
    impulses,
    irregular_responses,
    history='ragged',
    quiet=TRUE
)
stopifnot(length(irregular$terms) == 2L)
stopifnot(identical(irregular$specification[[1L]]$predictor, '1'))
stopifnot(isTRUE(irregular$specification[[1L]]$constant))
stopifnot(isTRUE(irregular$specification[[1L]]$implicit))
stopifnot(identical(irregular$identifiability$rate$status, 'retained'))

# A physical all-ones predictor is an explicit duplicate of the implicit rate;
# the structural default yields to the explicitly named scientific term.
ones_impulses <- impulses
ones_impulses$ones <- 1
duplicate_warning <- NULL
duplicate <- withCallingHandlers(
    prepare_cdrgam(
        response ~ irf(ones, window=c(0, 2), k=3),
        ones_impulses,
        irregular_responses,
        history='ragged',
        quiet=TRUE
    ),
    warning=function(w) {
        duplicate_warning <<- conditionMessage(w)
        invokeRestart('muffleWarning')
    }
)
stopifnot(grepl('duplicates explicit term ones', duplicate_warning, fixed=TRUE))
stopifnot(length(duplicate$terms) == 1L)
stopifnot(identical(names(duplicate$terms), 'ones'))

# Constant history counts also trigger the generic mgcv centering constraint
# for tensor IRFs, composed with their predictor-deviation constraint.
tensor <- prepare_cdrgam(
    response ~ irf(x, window=c(0, 2), k=c(3, 3), nonlinear=TRUE) - irf(1),
    impulses,
    responses,
    history='ragged',
    quiet=TRUE
)
stopifnot('predictor-deviation-centering' %in% tensor$terms[[1L]]$constraints)
stopifnot('linear-functional-centering' %in% tensor$terms[[1L]]$constraints)
stopifnot(all(is.finite(tensor$terms[[1L]]$X)))
