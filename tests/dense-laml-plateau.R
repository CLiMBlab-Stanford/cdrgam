library(cdrgam)

# A stall is accepted only on a flat plateau with a negligible gradient.
plateau <- cdrgam:::.cdrgam_dense_laml_plateau
flat <- function(parameters) 10 + exp(-parameters[[1L]]) + parameters[[2L]]^2
accepted <- plateau(flat, c(20, 0), flat(c(20, 0)), penalty_count=2L)
stopifnot(isTRUE(accepted$accepted), identical(accepted$plateau, 1L))

# A remaining slope is a failure even with a parameter on the plateau.
sloped <- plateau(flat, c(20, 1), flat(c(20, 1)), penalty_count=2L)
stopifnot(!sloped$accepted, abs(sloped$gradient[[2L]] - 2) < 1e-6)

# A flat interior stop is not a plateau.
interior <- plateau(flat, c(3, 0), flat(c(3, 0)), penalty_count=2L)
stopifnot(!interior$accepted, !length(interior$plateau))
