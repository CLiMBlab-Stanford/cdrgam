library(cdrgam)

facts <- methods_facts()
stopifnot(identical(facts$method, 'fREML'), identical(facts$engine, 'bam'),
    identical(facts$citations$criterion, 'wood2015bam'))
stopifnot(identical(methods_facts(engine='gam')$method, 'GCV.Cp'))
stopifnot(identical(methods_facts(engine='gam', family=mgcv::gaulss())$method, 'REML'))
sparse <- methods_facts(backend='sparse')
stopifnot(identical(sparse$method, 'REML'), is.null(sparse$engine),
    identical(sparse$citations, list(criterion='wood2011fast', estimator=character())))
stopifnot(identical(methods_facts(backend='block', method='fREML')$method, 'REML'))
stopifnot(identical(methods_facts(backend='sparse', method='GCV.Cp',
    family=stats::poisson())$criterion, 'the unbiased risk estimator'))
stopifnot(inherits(try(methods_facts(backend='block', method='GCV.Cp'), silent=TRUE),
    'try-error'))
stopifnot(setequal(methods_facts()$term_tests$method,
    c('Wald', 'Wood (2013)', 'Wood (2013b) random effect')))

# The resolved criterion matches what a native fit records.
set.seed(1)
data <- data.frame(x=stats::runif(80))
data$y <- sin(2 * pi * data$x) + stats::rnorm(80, sd=0.2)
for (engine in c('bam', 'gam')) {
    fit <- get(engine, asNamespace('mgcv'))(y ~ s(x), data=data)
    # mgcv records GCV.Cp as GCV when the scale parameter is estimated.
    recorded <- sub('^GCV$', 'GCV.Cp', fit$method)
    stopifnot(identical(recorded, methods_facts(engine=engine)$method))
}
cat('methods-facts tests passed\n')
