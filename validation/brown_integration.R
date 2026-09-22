# Brown-corpus integration test with full subject-specific IRFs and crossed
# subject/item random intercepts. Items are text positions and do not have
# enough replication to support item-specific IRFs.
#
# By default this uses all subjects from three documents, preserving the
# crossed replication needed to identify subject/item variance components.
# Set CDRGAM_BROWN_FULL=1 to use every filtered training row.
# CDRGAM_BROWN_MODELS can select a comma-separated subset of rate, linear,
# nonlinear, and nonstationary.

library(cdrgam)

script_started <- proc.time()[['elapsed']]

config_path <- Sys.getenv('CDRGAM_BROWN_CONFIG', 'brown.yml')
if (!file.exists(config_path)) stop('Brown config not found: ', config_path)
if (!requireNamespace('yaml', quietly=TRUE)) {
    stop('The Brown integration test requires the yaml package')
}
config <- yaml::read_yaml(config_path)

read_stream <- function(path, separator) {
    if (!file.exists(path)) stop('Brown data file not found: ', path)
    utils::read.csv(path, sep=separator, header=TRUE)
}

apply_brown_filters <- function(data, filters) {
    keep <- rep.int(TRUE, nrow(data))
    for (filter in filters) {
        if (!is.null(filter$column)) {
            keep <- keep & do.call(
                match.fun(filter$fun),
                c(list(data[[filter$column]]), list(filter$args))
            )
        } else if (!is.null(filter$factor)) {
            eligible <- names(which(table(data[[filter$factor]][keep]) >= filter$min))
            keep <- keep & as.character(data[[filter$factor]]) %in% eligible
        } else {
            stop('Unsupported filter in Brown metadata')
        }
    }
    data[keep, , drop=FALSE]
}

prepare_response_stream <- function(data, training_levels=NULL) {
    # Brown calls within-document position Word_Number. Expose the name used
    # by the model specification and fail if the source is unexpectedly absent.
    if (!('Word_Number' %in% names(data))) {
        stop('Brown data do not contain the expected Word_Number/docpos field')
    }
    data$docpos <- data$Word_Number
    data$item <- interaction(data$docid, data$docpos, drop=TRUE, sep=':')
    data$experiment_time <- ave(
        data$time,
        interaction(data$subject, data$docid, drop=TRUE),
        FUN=function(value) {
            span <- diff(range(value))
            if (span > 0) (value - min(value)) / span else rep.int(0, length(value))
        }
    )
    if (is.null(training_levels)) {
        data$subject <- droplevels(factor(data$subject))
        data$item <- droplevels(factor(data$item))
    } else {
        data$subject <- factor(data$subject, levels=training_levels$subject)
        data$item <- factor(data$item, levels=training_levels$item)
    }
    data
}

separator <- config$data$sep
impulses <- read_stream(config$data$X_train, separator)
train <- apply_brown_filters(
    read_stream(config$data$Y_train, separator),
    config$data$filters
)
validation <- apply_brown_filters(
    read_stream(config$data$Y_val, separator),
    config$data$filters
)

full_run <- identical(Sys.getenv('CDRGAM_BROWN_FULL', '0'), '1')
if (!full_run) {
    subject_count <- as.integer(Sys.getenv('CDRGAM_BROWN_SUBJECTS', '0'))
    document_count <- as.integer(Sys.getenv('CDRGAM_BROWN_DOCUMENTS', '3'))
    available_subjects <- sort(unique(as.character(train$subject)))
    selected_subjects <- if (subject_count <= 0L) available_subjects else
        available_subjects[seq_len(min(subject_count, length(available_subjects)))]
    selected_documents <- sort(unique(as.character(train$docid)))[seq_len(
        min(document_count, length(unique(train$docid)))
    )]
    select <- function(data) {
        data[
            as.character(data$subject) %in% selected_subjects &
                as.character(data$docid) %in% selected_documents,
            ,
            drop=FALSE
        ]
    }
    impulses <- select(impulses)
    train <- select(train)
    validation <- select(validation)
}

train <- prepare_response_stream(train)
training_levels <- list(
    subject=levels(train$subject),
    item=levels(train$item)
)
validation <- prepare_response_stream(validation, training_levels)
impulses$subject <- factor(impulses$subject, levels=training_levels$subject)

stopifnot(
    nrow(train) > 0L,
    nlevels(train$subject) >= 2L,
    nlevels(train$item) >= 2L,
    all(is.finite(train$fdur))
)
data_seconds <- proc.time()[['elapsed']] - script_started

predictors <- as.character(config$preds$all)
missing_predictors <- setdiff(predictors, names(impulses))
if (length(missing_predictors)) {
    stop('Brown impulse predictors are missing: ', paste(missing_predictors, collapse=', '))
}
window <- c(0, as.numeric(config$data$t_delta_cutoff))
history_length <- as.integer(config$data$history_length)
lag_k <- as.integer(Sys.getenv('CDRGAM_BROWN_LAG_K', '6'))
value_k <- as.integer(Sys.getenv('CDRGAM_BROWN_VALUE_K', '3'))
varying_k <- as.integer(Sys.getenv('CDRGAM_BROWN_VARYING_K', '4'))

format_irf <- function(
        predictor,
        group=NULL,
        nonlinear=FALSE,
        varying=NULL
) {
    arguments <- c(
        predictor,
        sprintf('window=c(%s,%s)', window[[1L]], window[[2L]]),
        if (nonlinear) sprintf('k=c(%d,%d)', lag_k, value_k) else
            if (!is.null(varying)) sprintf('k=c(%d,%d)', lag_k, varying_k) else
                sprintf('k=%d', lag_k),
        if (nonlinear) 'nonlinear=TRUE',
        if (!is.null(varying)) paste0('varying=', varying),
        if (!is.null(group)) paste0('group=', group)
    )
    paste0('irf(', paste(arguments, collapse=', '), ')')
}

with_subject_random_irfs <- function(term) {
    c(
        term(NULL),
        term('subject')
    )
}

rate_terms <- with_subject_random_irfs(function(group) {
    format_irf('1', group=group)
})
linear_terms <- unlist(lapply(predictors, function(predictor) {
    with_subject_random_irfs(function(group) {
        format_irf(predictor, group=group)
    })
}), use.names=FALSE)

# A nonlinear function on a binary domain has only two estimable values and is
# algebraically equivalent to a rate plus a linear binary IRF. Keep that term
# in its identified parameterization; all non-binary predictors use tensors.
unique_values <- vapply(
    impulses[predictors],
    function(value) length(unique(value[is.finite(value)])),
    integer(1)
)
binary_predictors <- predictors[unique_values < 3L]
continuous_predictors <- setdiff(predictors, binary_predictors)
nonlinear_terms <- c(
    unlist(lapply(binary_predictors, function(predictor) {
        with_subject_random_irfs(function(group) {
            format_irf(predictor, group=group)
        })
    }), use.names=FALSE),
    unlist(lapply(continuous_predictors, function(predictor) {
        with_subject_random_irfs(function(group) {
            format_irf(predictor, group=group, nonlinear=TRUE)
        })
    }), use.names=FALSE)
)
varying_rate_terms <- with_subject_random_irfs(function(group) {
    format_irf('1', group=group, varying='experiment_time')
})

ordinary_terms <- c("s(subject, bs='re')", "s(item, bs='re')")
make_formula <- function(terms) {
    stats::as.formula(paste('fdur ~', paste(c(ordinary_terms, terms), collapse=' + ')))
}
formulas <- list(
    rate=make_formula(rate_terms),
    linear=make_formula(c(rate_terms, linear_terms)),
    nonlinear=make_formula(c(rate_terms, nonlinear_terms)),
    nonstationary=make_formula(c(
        rate_terms,
        varying_rate_terms,
        nonlinear_terms
    ))
)
requested <- strsplit(
    Sys.getenv('CDRGAM_BROWN_MODELS', paste(names(formulas), collapse=',')),
    ',',
    fixed=TRUE
)[[1L]]
requested <- trimws(requested)
if (!length(requested) || any(!requested %in% names(formulas))) {
    stop('CDRGAM_BROWN_MODELS must select: ', paste(names(formulas), collapse=', '))
}

output_directory <- file.path('validation', 'output')
dir.create(output_directory, recursive=TRUE, showWarnings=FALSE)
run_label <- Sys.getenv(
    'CDRGAM_BROWN_RUN_LABEL',
    if (full_run) 'full_subject_irf' else paste0(
        'pilot_subject_irf_', nlevels(train$subject), 'subj_',
        length(unique(train$docid)), 'doc'
    )
)
save_fits <- identical(Sys.getenv('CDRGAM_BROWN_SAVE_FITS', '0'), '1')
use_checkpoints <- identical(Sys.getenv('CDRGAM_BROWN_CHECKPOINT', '1'), '1')
solver_trace <- as.integer(Sys.getenv('CDRGAM_BROWN_SOLVER_TRACE', '1'))
fit_backend <- Sys.getenv('CDRGAM_BROWN_BACKEND', 'sparse')
if (!identical(fit_backend, 'sparse')) {
    stop('CDRGAM_BROWN_BACKEND must be sparse')
}
gradient_default <- 'auto'
outer_optimizer_default <- 'lbfgsb'

fits <- list()
metrics <- vector('list', length(requested))
process_peak_rss_mib <- function() {
    status <- readLines('/proc/self/status', warn=FALSE)
    value <- sub('^VmHWM:\\s*([0-9]+)\\s+kB.*$', '\\1',
        status[grepl('^VmHWM:', status)])
    if (!length(value)) return(NA_real_)
    as.numeric(value[[1L]]) / 1024
}
for (i in seq_along(requested)) {
    model_name <- requested[[i]]
    gc(reset=TRUE)
    model_started <- proc.time()[['elapsed']]
    message('Preparing Brown model: ', model_name)
    preparation_elapsed <- system.time({
        design <- prepare_cdrgam(
            formulas[[model_name]],
            impulses,
            train,
            series=c('subject', 'docid'),
            history='auto',
            history_length=history_length,
            chunk_size=as.integer(Sys.getenv('CDRGAM_BROWN_CHUNK_SIZE', '10000')),
            quiet=FALSE
        )
    })[['elapsed']]
    stopifnot(
        all(vapply(design$plan, `[[`, numeric(1), 'maximum_history') <=
            history_length),
        identical(design$stream$history_length, history_length),
        !any(vapply(
            design$terms,
            function(term) identical(term$group, 'item'),
            logical(1)
        )),
        any(vapply(
            design$terms,
            function(term) identical(term$group, 'subject'),
            logical(1)
        ))
    )
    message('Fitting Brown model: ', model_name)
    checkpoint_path <- if (use_checkpoints) file.path(
        output_directory,
        paste0(
            'brown_', run_label, '_', model_name,
            '_optimizer_checkpoint.rds'
        )
    ) else NULL
    fit_elapsed <- system.time({
        fit <- cdrgam.fit(
            design,
            backend=fit_backend,
            family=stats::gaussian(),
            method='REML',
            checkpoint=checkpoint_path,
            solver_trace=solver_trace,
            sparse_control=list(
                gradient=Sys.getenv(
                    'CDRGAM_BROWN_GRADIENT',
                    gradient_default
                ),
                hessian=Sys.getenv(
                    'CDRGAM_BROWN_HESSIAN',
                    'gradient'
                ),
                outer_optimizer=Sys.getenv(
                    'CDRGAM_BROWN_OUTER_OPTIMIZER',
                    outer_optimizer_default
                ),
                optimizer_maxit=as.integer(Sys.getenv(
                    'CDRGAM_BROWN_OPTIMIZER_MAXIT',
                    '300'
                )),
                optimizer_gradient_tolerance=as.numeric(Sys.getenv(
                    'CDRGAM_BROWN_OPTIMIZER_GRADIENT_TOLERANCE',
                    '1e-4'
                )),
                boundary_action=Sys.getenv(
                    'CDRGAM_BROWN_BOUNDARY_ACTION',
                    'report'
                ),
                boundary_log_sp=as.numeric(Sys.getenv(
                    'CDRGAM_BROWN_BOUNDARY_LOG_SP',
                    '12'
                )),
                gradient_probes=as.integer(Sys.getenv(
                    'CDRGAM_BROWN_GRADIENT_PROBES',
                    '12'
                )),
                gradient_cores=as.integer(Sys.getenv(
                    'CDRGAM_BROWN_GRADIENT_CORES',
                    '1'
                )),
                schur='always',
                crossprod_chunk_size=as.integer(Sys.getenv(
                    'CDRGAM_BROWN_CROSSPROD_CHUNK_SIZE',
                    '5000'
                ))
            )
        )
    })[['elapsed']]
    stopifnot(
        all(is.finite(coef(fit))),
        all(is.finite(fitted(fit))),
        identical(fit$family$family, 'gaussian'),
        identical(fit$family$link, 'identity')
    )
    if (!isTRUE(fit$converged)) {
        warning(
            'Brown model ', model_name, ' did not meet the optimizer ',
            'convergence criterion; exporting the finite fit and diagnostics: ',
            fit$sparse$convergence$message,
            call.=FALSE
        )
    }
    validation_predictors <- validation[names(validation) != 'fdur']
    prediction_elapsed <- system.time({
        prediction <- predict(
            fit,
            list(impulses=impulses, responses=validation_predictors)
        )
    })[['elapsed']]
    residual <- validation$fdur - prediction
    if (save_fits) {
        saveRDS(
            fit,
            file.path(
                output_directory,
                paste0('brown_', run_label, '_', model_name, '_fit.rds')
            ),
            compress=FALSE
        )
    }
    if (identical(model_name, 'rate')) {
        lag_grid <- seq(window[[1L]], window[[2L]], length.out=201L)
        population_curve <- estimate_irf(
            fit,
            term='irf(1)',
            lag=lag_grid
        )
        subject_curve <- estimate_irf(
            fit,
            term='irf(1)|subject',
            lag=lag_grid,
            se=FALSE
        )
        population_at <- function(lag) {
            population_curve$estimate[match(lag, population_curve$lag)]
        }
        population_curve$curve <- 'population'
        population_curve$total_estimate <- population_curve$estimate
        subject_curve$curve <- 'subject total'
        subject_curve$total_estimate <- subject_curve$estimate +
            population_at(subject_curve$lag)
        rate_curves <- rbind(population_curve, subject_curve)
        utils::write.csv(
            rate_curves,
            file.path(
                output_directory,
                paste0('brown_', run_label, '_rate_irfs.csv')
            ),
            row.names=FALSE
        )
        plot_path <- file.path(
            output_directory,
            paste0('brown_', run_label, '_rate_irfs.png')
        )
        grDevices::png(plot_path, width=1100, height=1000)
        graphics::par(mfrow=c(2, 1), mar=c(4.5, 4.5, 3, 1))
        graphics::plot(
            population_curve$lag,
            population_curve$estimate,
            type='l',
            lwd=4,
            xlab='Delay (seconds)',
            ylab='Rate IRF',
            ylim=range(rate_curves$total_estimate),
            main=paste('Brown', run_label, 'rate IRF totals')
        )
        if (all(is.finite(population_curve$se))) {
            graphics::polygon(
                c(population_curve$lag, rev(population_curve$lag)),
                c(
                    population_curve$estimate - 1.96 * population_curve$se,
                    rev(population_curve$estimate + 1.96 * population_curve$se)
                ),
                col=grDevices::adjustcolor('grey60', alpha.f=0.3),
                border=NA
            )
        }
        subject_colors <- grDevices::hcl.colors(
            length(unique(subject_curve$group)),
            'Dynamic'
        )
        for (j in seq_along(unique(subject_curve$group))) {
            subject <- unique(subject_curve$group)[[j]]
            rows <- subject_curve$group == subject
            graphics::lines(
                subject_curve$lag[rows],
                subject_curve$total_estimate[rows],
                col=grDevices::adjustcolor(subject_colors[[j]], alpha.f=0.45)
            )
        }
        graphics::lines(
            population_curve$lag,
            population_curve$estimate,
            lwd=4,
            col='black'
        )
        graphics::legend(
            'topright',
            c('population', 'subject totals'),
            col=c('black', subject_colors[[1L]]),
            lwd=c(4, 1),
            bty='n'
        )
        deviation_range <- range(subject_curve$estimate)
        graphics::plot(
            range(lag_grid),
            deviation_range,
            type='n',
            xlab='Delay (seconds)',
            ylab='Random rate deviation',
            main='Grouped deviations (separate vertical scale)'
        )
        graphics::abline(h=0, lty=3, col='grey40')
        for (j in seq_along(unique(subject_curve$group))) {
            subject <- unique(subject_curve$group)[[j]]
            rows <- subject_curve$group == subject
            graphics::lines(
                subject_curve$lag[rows],
                subject_curve$estimate[rows],
                col=grDevices::adjustcolor(subject_colors[[j]], alpha.f=0.55)
            )
        }
        grDevices::dev.off()
    }
    if (!identical(model_name, 'rate')) {
        labels <- fit$cdrgam$term_labels
        population_labels <- labels[!grepl('|', labels, fixed=TRUE)]
        selected_levels <- function(info) {
            levels <- info$group_levels
            if (identical(info$group, 'subject')) return(levels)
            levels[unique(round(seq(
                1,
                length(levels),
                length.out=min(12L, length(levels))
            )))]
        }
        irf_output <- list()
        population_by_label <- list()
        for (label in population_labels) {
            population <- estimate_irf(
                fit,
                term=label,
                n=201L,
                n_predictor=17L
            )
            population$curve <- 'population'
            population$total_estimate <- population$estimate
            population_by_label[[label]] <- population
            irf_output[[length(irf_output) + 1L]] <- population

            filename <- gsub('[^A-Za-z0-9]+', '_', label)
            plot_path <- file.path(
                output_directory,
                paste0(
                    'brown_', run_label, '_', model_name, '_',
                    filename, '.png'
                )
            )
            grDevices::png(plot_path, width=1050, height=750)
            if (all(is.na(population$predictor))) {
                graphics::plot(
                    population$lag,
                    population$estimate,
                    type='l',
                    lwd=4,
                    xlab='Delay (seconds)',
                    ylab='IRF',
                    main=paste(model_name, label)
                )
            } else {
                lag_values <- sort(unique(population$lag))
                predictor_values <- sort(unique(population$predictor))
                surface <- matrix(
                    NA_real_,
                    nrow=length(lag_values),
                    ncol=length(predictor_values)
                )
                surface[cbind(
                    match(population$lag, lag_values),
                    match(population$predictor, predictor_values)
                )] <- population$estimate
                graphics::image(
                    lag_values,
                    predictor_values,
                    surface,
                    xlab='Delay (seconds)',
                    ylab=if (startsWith(label, 'irf(1)~')) {
                        'Normalized experiment time'
                    } else 'Predictor value',
                    main=paste(model_name, label),
                    col=grDevices::hcl.colors(80L, 'Blue-Red 3')
                )
                graphics::contour(
                    lag_values,
                    predictor_values,
                    surface,
                    add=TRUE,
                    drawlabels=FALSE
                )
            }
            grDevices::dev.off()
        }
        grouped_labels <- labels[grepl('|', labels, fixed=TRUE)]
        for (label in grouped_labels) {
            population_label <- sub('\\|[^|]+$', '', label)
            population <- population_by_label[[population_label]]
            if (is.null(population)) next
            info <- fit$cdrgam$terms[[match(label, labels)]]
            deviation <- estimate_irf(
                fit,
                term=label,
                lag=sort(unique(population$lag)),
                predictor=if (all(is.na(population$predictor))) NULL else
                    sort(unique(population$predictor)),
                group=selected_levels(info),
                se=FALSE
            )
            key <- paste(
                format(population$lag, digits=17),
                format(population$predictor, digits=17),
                sep=':'
            )
            deviation_key <- paste(
                format(deviation$lag, digits=17),
                format(deviation$predictor, digits=17),
                sep=':'
            )
            deviation$curve <- paste0(info$group, ' total')
            deviation$total_estimate <- deviation$estimate +
                population$estimate[match(deviation_key, key)]
            irf_output[[length(irf_output) + 1L]] <- deviation
        }
        utils::write.csv(
            do.call(rbind, irf_output),
            file.path(
                output_directory,
                paste0('brown_', run_label, '_', model_name, '_irfs.csv')
            ),
            row.names=FALSE
        )
    }
    evaluations <- fit$sparse$convergence$evaluations
    metrics[[i]] <- data.frame(
        model=model_name,
        train_rows=nrow(train),
        validation_rows=nrow(validation),
        subjects=nlevels(train$subject),
        items=nlevels(train$item),
        coefficients=length(coef(fit)),
        smoothing_parameters=length(fit$sp),
        backend=fit$cdrgam$backend,
        gradient=fit$sparse$gradient,
        outer_optimizer=fit$sparse$outer_optimizer,
        optimizer_message=fit$sparse$convergence$message,
        projected_gradient=if (is.null(fit$sparse$optimizer_progress)) {
            NA_real_
        } else fit$sparse$optimizer_progress$projected_gradient_max,
        gradient_ratio=if (is.null(fit$sparse$optimizer_progress)) {
            NA_real_
        } else fit$sparse$optimizer_progress$gradient_ratio,
        rejected_steps=if (is.null(fit$sparse$optimizer_progress)) {
            NA_integer_
        } else fit$sparse$optimizer_progress$rejected_steps,
        curvature_resets=if (is.null(fit$sparse$optimizer_progress)) {
            NA_integer_
        } else fit$sparse$optimizer_progress$curvature_resets,
        data_seconds=unname(data_seconds),
        preparation_seconds=unname(preparation_elapsed),
        fit_seconds=unname(fit_elapsed),
        prediction_seconds=unname(prediction_elapsed),
        total_seconds=proc.time()[['elapsed']] - model_started,
        process_peak_rss_mib=process_peak_rss_mib(),
        design_mib=as.numeric(object.size(design)) / 1024^2,
        fit_mib=as.numeric(object.size(fit)) / 1024^2,
        design_nonzeros=fit$sparse$nnzero,
        system_nonzeros=fit$sparse$system_nnzero,
        factor_nonzeros=fit$sparse$factor_nonzeros,
        objective_evaluations=
            fit$sparse$convergence$total_objective_evaluations,
        numeric_factorizations=fit$sparse$numeric_updates,
        hessian_factorizations=fit$sparse$hessian_evaluations,
        optimizer_function_evaluations=
            unname(as.numeric(evaluations[[1L]])),
        gradient_evaluations=if (length(evaluations) >= 2L) {
            unname(as.numeric(evaluations[[2L]]))
        } else NA_real_,
        train_rmse=sqrt(mean(residuals(fit)^2)),
        validation_rmse=sqrt(mean(residual^2)),
        validation_mae=mean(abs(residual)),
        logLik=as.numeric(logLik(fit)),
        AIC=as.numeric(AIC(fit)),
        schur_group=fit$sparse$schur_group,
        converged=fit$converged,
        stringsAsFactors=FALSE
    )
    fits[[model_name]] <- fit
    rm(design)
    gc()
}

metrics <- do.call(rbind, metrics)
print(metrics, row.names=FALSE)
utils::write.csv(
    metrics,
    file.path(
        output_directory,
        if (length(requested) == 1L) {
            paste0('brown_', run_label, '_', requested, '_metrics.csv')
        } else 'brown_integration_metrics.csv'
    ),
    row.names=FALSE
)
saveRDS(
    list(
        formulas=formulas[requested],
        metrics=metrics,
        binary_predictors=binary_predictors,
        data=list(
            full=full_run,
            train_rows=nrow(train),
            validation_rows=nrow(validation),
            subjects=levels(train$subject),
            documents=sort(unique(as.character(train$docid)))
        )
    ),
    file.path(
        output_directory,
        if (length(requested) == 1L) {
            paste0('brown_', run_label, '_', requested, '_summary.rds')
        } else 'brown_integration_summary.rds'
    )
)
