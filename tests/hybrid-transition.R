library(cdrgam)

policy <- getFromNamespace('.hybrid_rejection_policy', 'cdrgam')
transition_condition <- getFromNamespace(
    '.hybrid_transition_condition', 'cdrgam'
)

productive <- policy(
    rejected_solves=3L,
    rejected_seconds=600,
    reference_seconds=c(250, 300, 350)
)
time_exhausted <- policy(
    rejected_solves=5L,
    rejected_seconds=1200,
    reference_seconds=c(250, 300, 350)
)
count_exhausted <- policy(
    rejected_solves=12L,
    rejected_seconds=0,
    reference_seconds=1
)
caught <- tryCatch(
    stop(transition_condition(time_exhausted)),
    cdrgam_hybrid_transition=identity
)

stopifnot(
    !productive$transition,
    identical(productive$seconds_limit, 1200),
    time_exhausted$transition,
    grepl('cost budget', time_exhausted$reason, fixed=TRUE),
    count_exhausted$transition,
    grepl('safety limit', count_exhausted$reason, fixed=TRUE),
    inherits(caught, 'cdrgam_hybrid_transition'),
    identical(caught$rejected_solves, 5L)
)

cat('hybrid-transition: ok\n')
