library(cdrgam)

set.seed(20261012)
subjects <- factor(rep(letters[1:4], each=24))
impulses <- data.frame(
    subject=subjects,
    time=unlist(lapply(seq_len(4), function(index) seq(0, 4.6, length.out=24))),
    x=rnorm(96)
)
responses <- data.frame(
    subject=factor(rep(letters[1:4], each=18), levels=levels(subjects)),
    time=unlist(lapply(seq_len(4), function(index) seq(0.3, 4.9, length.out=18))),
    y=rnorm(72)
)

design <- prepare_cdrgam(
    y ~ s(time, k=5) +
        irf(x, window=c(0, 0.8), k_l=5) +
        irf(x, window=c(0, 0.8), k_l=5, group=subject) -
        irf(1),
    impulses=impulses,
    responses=responses,
    series='subject',
    quiet=TRUE
)
specification <- as_mgcv_spec(
    design,
    family=gaussian(),
    method='REML'
)
stopifnot(identical(formula(design, type='mgcv'), specification$formula))
formula_text <- paste(deparse(specification$formula), collapse='')
stopifnot(
    grepl('irf_x', formula_text, fixed=TRUE),
    grepl('irf_x_by_subject', formula_text, fixed=TRUE),
    !grepl('cdr_term_', formula_text, fixed=TRUE),
    !grepl('bs = "cdr"', formula_text, fixed=TRUE),
    all(names(specification$paraPen) %in% names(specification$data))
)

reference <- do.call(mgcv::gam, specification)
native <- cdrgam.fit(
    design,
    family=gaussian(),
    method='REML',
    backend='mgcv',
    engine='gam'
)
legacy <- native
legacy$cdrgam$formula$mgcv <- legacy$cdrgam$formula$computational
legacy_formula_text <- paste(
    deparse(formula(legacy, type='mgcv')),
    collapse=''
)
stopifnot(
    isTRUE(all.equal(
        unname(fitted(reference)),
        unname(fitted(native)),
        tolerance=2e-6
    )),
    isTRUE(all.equal(
        sort(unname(reference$sp)),
        sort(unname(native$sp)),
        tolerance=2e-5
    )),
    grepl('irf_x', legacy_formula_text, fixed=TRUE),
    !grepl('cdr_term_', legacy_formula_text, fixed=TRUE)
)

bad <- try(as_mgcv_spec(native), silent=TRUE)
stopifnot(inherits(bad, 'try-error'))
