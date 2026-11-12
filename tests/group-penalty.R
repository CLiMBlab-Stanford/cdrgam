library(cdrgam)

set.seed(20261014)
groups <- factor(rep(letters[1:3], each=18))
impulses <- data.frame(
    subject=groups,
    time=rep(seq(0, 3.4, length.out=18), 3),
    x=rnorm(54)
)
responses <- data.frame(
    subject=factor(rep(letters[1:3], each=12), levels=levels(groups)),
    time=rep(seq(0.25, 3.6, length.out=12), 3),
    y=rnorm(36)
)

prepare <- function(group_penalty='re', nonlinear=FALSE) {
    formula <- if (nonlinear) {
        y ~ irf(
            x, window=c(0, 0.9), k_l=5, k_p=4,
            group=subject, group_penalty=group_penalty
        ) - irf(1)
    } else {
        y ~ irf(
            x, window=c(0, 0.9), k_l=5,
            group=subject, group_penalty=group_penalty
        ) - irf(1)
    }
    prepare_cdrgam(
        formula,
        impulses=impulses,
        responses=responses,
        series='subject',
        quiet=TRUE
    )
}

default <- prepare_cdrgam(
    y ~ irf(x, window=c(0, 0.9), k_l=5, group=subject) - irf(1),
    impulses=impulses,
    responses=responses,
    series='subject',
    quiet=TRUE
)
re <- prepare('re')
t2 <- prepare('t2')
fs <- prepare('fs')

stopifnot(
    identical(re$terms[[1L]]$group_penalty, 're'),
    identical(t2$terms[[1L]]$group_penalty, 't2'),
    identical(fs$terms[[1L]]$group_penalty, 'fs'),
    length(re$terms[[1L]]$S) == 2L,
    length(t2$terms[[1L]]$S) == 2L,
    length(fs$terms[[1L]]$S) == 3L,
    isTRUE(all.equal(default$terms[[1L]]$S, re$terms[[1L]]$S)),
    grepl(
        'group_penalty = "fs"',
        paste(deparse(formula(fs, type='normalized')), collapse=''),
        fixed=TRUE
    )
)

check_native_penalties <- function(design, full) {
    term <- design$terms[[1L]]
    specification <- design$specification[[1L]]
    axes <- if (!is.null(term$axis)) {
        lapply(term$axis, function(axis) list(
            name=axis$internal,
            values=axis$grid,
            k=axis$k,
            bs=axis$bs,
            knots=axis$knots
        ))
    } else {
        list(list(
            name=term$basis$term[[1L]],
            values=term$knots,
            k=specification$k_l,
            bs=specification$bs_l,
            knots=if (identical(specification$bs_l, 'ps')) NULL else term$knots
        ))
    }
    axis_names <- vapply(axes, `[[`, character(1), 'name')
    reference_group <- 'reference_group'
    construction <- do.call(expand.grid, c(
        stats::setNames(lapply(axes, `[[`, 'values'), axis_names),
        stats::setNames(list(factor(c('a', 'b', 'c'))), reference_group),
        list(KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
    ))
    constructor <- if (is.null(full)) quote(mgcv::te) else quote(mgcv::t2)
    controls <- list(
        k=c(vapply(axes, `[[`, integer(1), 'k'), 3L),
        bs=c(vapply(axes, `[[`, character(1), 'bs'), 're')
    )
    if (!is.null(full)) controls$full <- full
    call <- as.call(c(
        list(constructor),
        lapply(c(axis_names, reference_group), as.name),
        controls
    ))
    knots <- stats::setNames(lapply(axes, `[[`, 'knots'), axis_names)
    knots <- knots[!vapply(knots, is.null, logical(1))]
    native <- mgcv::smoothCon(
        eval(call), construction,
        knots=if (length(knots)) knots else NULL,
        absorb.cons=FALSE, scale.penalty=FALSE, n=nrow(construction)
    )[[1L]]
    base <- mgcv::PredictMat(
        term$basis, construction[axis_names], n=nrow(construction)
    )
    if (!is.null(term$transform)) base <- base %*% term$transform
    dimension <- ncol(base)
    grouped <- matrix(0, nrow(construction), 3L * dimension)
    for (index in 1:3) {
        rows <- which(as.integer(construction[[reference_group]]) == index)
        columns <- (index - 1L) * dimension + seq_len(dimension)
        grouped[rows, columns] <- base[rows, , drop=FALSE]
    }
    mapping <- qr.solve(native$X, grouped)
    expected <- lapply(native$S, function(penalty) {
        value <- crossprod(mapping, penalty %*% mapping)
        value / norm(value, type='I')
    })
    observed <- lapply(
        cdrgam:::.materialize_cdr_term(term)$S,
        function(penalty) as.matrix(penalty) / norm(penalty, type='I')
    )
    stopifnot(
        length(observed) == length(expected),
        all(vapply(seq_along(observed), function(index) {
            isTRUE(all.equal(observed[[index]], expected[[index]], tolerance=2e-8))
        }, logical(1)))
    )
}

check_native_penalties(re, full=NULL)
check_native_penalties(t2, full=FALSE)
check_native_penalties(fs, full=TRUE)

nonlinear_re <- prepare('re', nonlinear=TRUE)
nonlinear_t2 <- prepare('t2', nonlinear=TRUE)
nonlinear_fs <- prepare('fs', nonlinear=TRUE)
stopifnot(
    length(nonlinear_re$terms[[1L]]$S) == 3L,
    length(nonlinear_t2$terms[[1L]]$S) == 4L,
    length(nonlinear_fs$terms[[1L]]$S) == 9L
)
check_native_penalties(nonlinear_t2, full=FALSE)
check_native_penalties(nonlinear_fs, full=TRUE)
check_native_penalties(nonlinear_re, full=NULL)

invalid <- try(irf(x, group_penalty='fs'), silent=TRUE)
unknown <- try(irf(x, group=subject, group_penalty='unknown'), silent=TRUE)
stopifnot(inherits(invalid, 'try-error'), inherits(unknown, 'try-error'))

for (design in list(re, t2, fs)) {
    assembly <- cdrgam:::.fit_sparse_gaussian(
        design,
        family=gaussian(),
        method='REML',
        setup_only=TRUE
    )
    expected <- cdrgam:::.materialize_cdr_term(
        design$terms[[1L]], sparse=TRUE
    )
    assembled <- cdrgam:::.sparse_component_rows(
        assembly$matrices[[2L]],
        seq_len(nrow(responses))
    )
    stopifnot(
        inherits(
            assembly$matrices[[2L]],
            'cdrgam_sparse_grouped_component'
        ),
        isTRUE(all.equal(
            as.matrix(assembled),
            as.matrix(expected$X),
            tolerance=2e-12
        )),
        length(assembly$sp_names) >= length(expected$S)
    )
    fit <- cdrgam.fit(
        design,
        family=gaussian(),
        method='REML',
        backend='mgcv',
        engine='gam'
    )
    stopifnot(all(is.finite(fitted(fit))))
}

cat('Grouped IRF penalty checks passed.\n')
