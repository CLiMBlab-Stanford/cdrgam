#!/usr/bin/env Rscript

engine_files <- c(
    'R/block.R', 'R/compressed.R', 'R/distributional-block.R',
    'R/distributional-sparse.R', 'R/distributional.R', 'R/family.R',
    'R/formula.R', 'R/parallelism.R', 'R/pirls.R', 'R/prediction.R',
    'R/scaling.R', 'R/solver-control.R', 'R/sparse.R', 'src/init.c',
    'src/schur.c'
)

support_files <- c(
    'R/diagnostics.R', 'R/effects.R', 'R/marginaleffects.R', 'R/plot.R',
    'R/summary.R'
)

usage <- function(status=0L) {
    output <- c(
        'Usage: Rscript scripts/audit-code-size.R [OPTIONS]', '',
        'Report source-line counts for the active cdrgam package.', '',
        'Options:',
        '  --format FORMAT   text or markdown (default: text)',
        '  --output PATH     write the report to PATH instead of stdout',
        '  --root PATH       repository root (normally inferred)',
        '  --self-test       run parser checks, then exit',
        '  -h, --help        show this help'
    )
    writeLines(output, if (status == 0L) stdout() else stderr())
    quit(save='no', status=status)
}

parse_arguments <- function(arguments) {
    options <- list(format='text', output=NULL, root=NULL, self_test=FALSE)
    take_value <- function(index, name) {
        if (index == length(arguments)) {
            stop(name, ' requires a value', call.=FALSE)
        }
        arguments[[index + 1L]]
    }
    index <- 1L
    while (index <= length(arguments)) {
        argument <- arguments[[index]]
        if (argument %in% c('-h', '--help')) usage()
        if (identical(argument, '--self-test')) {
            options$self_test <- TRUE
        } else if (grepl('^--format=', argument)) {
            options$format <- sub('^--format=', '', argument)
        } else if (identical(argument, '--format')) {
            options$format <- take_value(index, '--format')
            index <- index + 1L
        } else if (grepl('^--output=', argument)) {
            options$output <- sub('^--output=', '', argument)
        } else if (identical(argument, '--output')) {
            options$output <- take_value(index, '--output')
            index <- index + 1L
        } else if (grepl('^--root=', argument)) {
            options$root <- sub('^--root=', '', argument)
        } else if (identical(argument, '--root')) {
            options$root <- take_value(index, '--root')
            index <- index + 1L
        } else {
            stop('Unknown argument: ', argument, call.=FALSE)
        }
        index <- index + 1L
    }
    options$format <- match.arg(options$format, c('text', 'markdown'))
    options
}

script_path <- function() {
    file_argument <- grep('^--file=', commandArgs(FALSE), value=TRUE)
    if (!length(file_argument)) return(NULL)
    normalizePath(sub('^--file=', '', file_argument[[1L]]), mustWork=TRUE)
}

repository_root <- function(explicit=NULL) {
    if (!is.null(explicit)) return(normalizePath(explicit, mustWork=TRUE))
    path <- script_path()
    if (!is.null(path)) return(dirname(dirname(path)))
    normalizePath('.', mustWork=TRUE)
}

line_counts <- function(classification) {
    levels <- c('source', 'comment', 'blank')
    values <- table(factor(classification, levels=levels))
    c(
        source=unname(values[['source']]),
        comment=unname(values[['comment']]),
        blank=unname(values[['blank']]),
        physical=length(classification)
    )
}

classify_r <- function(lines) {
    output <- rep.int('source', length(lines))
    output[grepl('^[[:space:]]*$', lines)] <- 'blank'
    output[grepl('^[[:space:]]*#', lines)] <- 'comment'
    output
}

classify_c <- function(lines) {
    output <- character(length(lines))
    in_block <- FALSE
    for (line_number in seq_along(lines)) {
        characters <- strsplit(lines[[line_number]], '', fixed=TRUE)[[1L]]
        position <- 1L
        has_source <- FALSE
        has_comment <- FALSE
        quote <- NULL
        escaped <- FALSE
        while (position <= length(characters)) {
            current <- characters[[position]]
            following <- if (position < length(characters)) {
                characters[[position + 1L]]
            } else ''
            if (in_block) {
                has_comment <- TRUE
                if (identical(current, '*') && identical(following, '/')) {
                    in_block <- FALSE
                    position <- position + 2L
                } else position <- position + 1L
                next
            }
            if (!is.null(quote)) {
                has_source <- TRUE
                if (escaped) {
                    escaped <- FALSE
                } else if (identical(current, '\\')) {
                    escaped <- TRUE
                } else if (identical(current, quote)) quote <- NULL
                position <- position + 1L
                next
            }
            if (current %in% c('"', "'")) {
                has_source <- TRUE
                quote <- current
                position <- position + 1L
            } else if (identical(current, '/') && identical(following, '/')) {
                has_comment <- TRUE
                break
            } else if (identical(current, '/') && identical(following, '*')) {
                has_comment <- TRUE
                in_block <- TRUE
                position <- position + 2L
            } else {
                if (!grepl('[[:space:]]', current)) has_source <- TRUE
                position <- position + 1L
            }
        }
        output[[line_number]] <- if (has_source) {
            'source'
        } else if (has_comment || in_block) {
            'comment'
        } else 'blank'
    }
    output
}

classify_file <- function(path) {
    lines <- readLines(path, warn=FALSE)
    extension <- tolower(tools::file_ext(path))
    classification <- switch(
        extension,
        r=classify_r(lines),
        c=classify_c(lines),
        h=classify_c(lines),
        stop('Unsupported source extension: ', path, call.=FALSE)
    )
    counts <- line_counts(classification)
    data.frame(
        file=path,
        language=if (identical(extension, 'r')) 'R' else 'C',
        source=counts[['source']], comment=counts[['comment']],
        blank=counts[['blank']], physical=counts[['physical']],
        stringsAsFactors=FALSE
    )
}

discover_implementation <- function(package_root) {
    r_files <- file.path('R', list.files(
        file.path(package_root, 'R'), pattern='\\.[Rr]$', full.names=FALSE
    ))
    native_files <- file.path('src', list.files(
        file.path(package_root, 'src'), pattern='\\.[cChH]$', full.names=FALSE
    ))
    sort(c(r_files, native_files))
}

discover_tests <- function(package_root) {
    sort(file.path('tests', list.files(
        file.path(package_root, 'tests'), pattern='\\.[Rr]$', full.names=FALSE
    )))
}

validate_manifest <- function(implementation) {
    classified <- c(engine_files, support_files)
    missing <- setdiff(classified, implementation)
    unclassified <- setdiff(implementation, classified)
    duplicated <- unique(classified[duplicated(classified)])
    messages <- character()
    if (length(missing)) messages <- c(messages, paste0(
        'Manifest entries do not exist: ', paste(missing, collapse=', ')
    ))
    if (length(unclassified)) messages <- c(messages, paste0(
        'Implementation files need an engine/support classification: ',
        paste(unclassified, collapse=', ')
    ))
    if (length(duplicated)) messages <- c(messages, paste0(
        'Manifest entries are duplicated: ', paste(duplicated, collapse=', ')
    ))
    if (length(messages)) stop(paste(messages, collapse='\n'), call.=FALSE)
    invisible(TRUE)
}

count_files <- function(package_root, paths) {
    output <- do.call(rbind, lapply(file.path(package_root, paths), classify_file))
    output$file <- paths
    rownames(output) <- NULL
    output
}

sum_counts <- function(table) {
    columns <- c('source', 'comment', 'blank', 'physical')
    if (!nrow(table)) return(stats::setNames(rep.int(0, length(columns)), columns))
    colSums(table[, columns, drop=FALSE])
}

scope_rows <- function(name, table, include_languages=TRUE) {
    counts <- sum_counts(table)
    output <- data.frame(
        scope=name, source=counts[['source']], comment=counts[['comment']],
        blank=counts[['blank']], physical=counts[['physical']],
        stringsAsFactors=FALSE
    )
    if (!include_languages) return(output)
    for (language in c('R', 'C')) {
        selected <- table[table$language == language, , drop=FALSE]
        if (!nrow(selected)) next
        counts <- sum_counts(selected)
        output <- rbind(output, data.frame(
            scope=paste0('  ', language), source=counts[['source']],
            comment=counts[['comment']], blank=counts[['blank']],
            physical=counts[['physical']], stringsAsFactors=FALSE
        ))
    }
    output
}

git_metadata <- function(root) {
    if (!nzchar(Sys.which('git'))) {
        return(list(revision='unavailable', state='unavailable'))
    }
    run_git <- function(arguments) suppressWarnings(system2(
        'git', c('-C', shQuote(root), arguments), stdout=TRUE, stderr=FALSE
    ))
    revision <- run_git(c('rev-parse', '--short', 'HEAD'))
    status <- attr(revision, 'status')
    if (!length(revision) || (!is.null(status) && status != 0L)) {
        return(list(revision='unavailable', state='unavailable'))
    }
    changes <- run_git(c('status', '--porcelain'))
    list(
        revision=revision[[1L]],
        state=if (length(changes)) 'dirty' else 'clean'
    )
}

format_integer <- function(value) format(
    as.integer(value), big.mark=',', scientific=FALSE, trim=TRUE
)

definitions <- c(
    'Definitions:',
    '- Core engine contains model preparation, fitting, optimization,',
    '  prediction, scaling, parallelism, and native numerical routines.',
    '- Supporting APIs and presentation contain diagnostics, effect',
    '  extraction, marginaleffects integration, plotting, and summaries.',
    '- Source counts exclude blank and comment-only lines. Lines containing',
    '  both source and an inline comment count as source.',
    '- Documentation, archived code, generated validation output, and the',
    '  separate cdrgam.cli repository are excluded.'
)

text_report <- function(root, metadata, scopes, counts) {
    formatted <- data.frame(
        scope=scopes$scope,
        source=format_integer(scopes$source),
        comment=format_integer(scopes$comment),
        blank=format_integer(scopes$blank),
        physical=format_integer(scopes$physical),
        stringsAsFactors=FALSE
    )
    widths <- c(
        scope=max(nchar(c('Scope', formatted$scope))),
        source=max(nchar(c('Source', formatted$source))),
        comment=max(nchar(c('Comment-only', formatted$comment))),
        blank=max(nchar(c('Blank', formatted$blank))),
        physical=max(nchar(c('Physical', formatted$physical)))
    )
    row <- function(scope, source, comment, blank, physical) sprintf(
        paste0('%-', widths[['scope']], 's  %', widths[['source']],
            's  %', widths[['comment']], 's  %', widths[['blank']],
            's  %', widths[['physical']], 's'),
        scope, source, comment, blank, physical
    )
    table <- row('Scope', 'Source', 'Comment-only', 'Blank', 'Physical')
    table <- c(table, paste(rep.int('-', nchar(table)), collapse=''))
    for (index in seq_len(nrow(formatted))) table <- c(table, row(
        formatted$scope[[index]], formatted$source[[index]],
        formatted$comment[[index]], formatted$blank[[index]],
        formatted$physical[[index]]
    ))
    c(
        'CDR-GAM code size audit',
        paste0('Repository: ', root),
        paste0('Revision: ', metadata$revision, ' (', metadata$state, ')'),
        paste0('Generated: ', format(Sys.time(), tz='UTC', usetz=TRUE)), '',
        table, '',
        paste0(
            'Core engine: ', format_integer(counts$engine[['physical']]),
            ' physical lines; ', format_integer(counts$engine[['source']]),
            ' nonblank, noncomment source lines.'
        ),
        paste0(
            'Active package implementation: ',
            format_integer(counts$implementation[['physical']]),
            ' physical lines; ',
            format_integer(counts$implementation[['source']]),
            ' nonblank, noncomment source lines.'
        ),
        paste0(
            'Standalone regression tests: ',
            format_integer(counts$tests[['physical']]), ' physical lines.'
        ), '', definitions
    )
}

markdown_report <- function(root, metadata, scopes, counts) {
    table <- c(
        '| Scope | Source | Comment-only | Blank | Physical |',
        '| --- | ---: | ---: | ---: | ---: |'
    )
    for (index in seq_len(nrow(scopes))) {
        label <- sub('^  ', '&nbsp;&nbsp;', scopes$scope[[index]])
        table <- c(table, paste0(
            '| ', label, ' | ', format_integer(scopes$source[[index]]),
            ' | ', format_integer(scopes$comment[[index]]),
            ' | ', format_integer(scopes$blank[[index]]),
            ' | ', format_integer(scopes$physical[[index]]), ' |'
        ))
    }
    markdown_definitions <- definitions
    markdown_definitions[[1L]] <- '## Definitions'
    markdown_definitions <- append(markdown_definitions, '', after=1L)
    c(
        '# CDR-GAM code size audit', '',
        paste0('- Repository: `', root, '`'),
        paste0('- Revision: `', metadata$revision, '` (', metadata$state, ')'),
        paste0('- Generated: ', format(Sys.time(), tz='UTC', usetz=TRUE)), '',
        table, '',
        paste0(
            'The core engine contains **',
            format_integer(counts$engine[['physical']]),
            ' physical lines** and **',
            format_integer(counts$engine[['source']]),
            ' nonblank, noncomment source lines**.'
        ), '',
        paste0(
            'The complete active package implementation contains **',
            format_integer(counts$implementation[['physical']]),
            ' physical lines** and **',
            format_integer(counts$implementation[['source']]),
            ' nonblank, noncomment source lines**. The standalone regression ',
            'suite adds **', format_integer(counts$tests[['physical']]),
            ' physical lines**.'
        ), '', markdown_definitions
    )
}

self_test <- function() {
    stopifnot(identical(
        classify_r(c('', '  # comment', 'x <- 1 # inline')),
        c('blank', 'comment', 'source')
    ))
    c_lines <- c(
        '', '// comment', 'int x = 1; // inline', '/* block', ' * comment */',
        'const char *s = "// source";'
    )
    stopifnot(identical(
        classify_c(c_lines),
        c('blank', 'comment', 'source', 'comment', 'comment', 'source')
    ))
    stopifnot(identical(
        unname(line_counts(classify_c(c_lines))), c(2L, 3L, 1L, 6L)
    ))
    invisible(TRUE)
}

main <- function() {
    options <- parse_arguments(commandArgs(trailingOnly=TRUE))
    if (isTRUE(options$self_test)) {
        self_test()
        writeLines('Code-size audit self-test passed')
        return(invisible(NULL))
    }
    root <- repository_root(options$root)
    package_root <- file.path(root, 'cdrgam')
    if (!file.exists(file.path(package_root, 'DESCRIPTION'))) {
        stop('Could not find the active cdrgam package below ', root, call.=FALSE)
    }
    implementation <- discover_implementation(package_root)
    validate_manifest(implementation)
    tests <- discover_tests(package_root)
    implementation_table <- count_files(package_root, implementation)
    engine_table <- implementation_table[
        implementation_table$file %in% engine_files, , drop=FALSE
    ]
    support_table <- implementation_table[
        implementation_table$file %in% support_files, , drop=FALSE
    ]
    test_table <- count_files(package_root, tests)
    scopes <- rbind(
        scope_rows('Core engine', engine_table),
        scope_rows('Supporting APIs/presentation', support_table),
        scope_rows('Active package implementation', implementation_table),
        scope_rows('Standalone regression tests', test_table, FALSE)
    )
    counts <- list(
        engine=sum_counts(engine_table),
        implementation=sum_counts(implementation_table),
        tests=sum_counts(test_table)
    )
    metadata <- git_metadata(root)
    report <- switch(
        options$format,
        text=text_report(root, metadata, scopes, counts),
        markdown=markdown_report(root, metadata, scopes, counts)
    )
    if (is.null(options$output)) {
        writeLines(report)
    } else {
        output <- path.expand(options$output)
        if (!dir.exists(dirname(output))) {
            stop('Output directory does not exist: ', dirname(output), call.=FALSE)
        }
        writeLines(report, output)
    }
    invisible(report)
}

tryCatch(
    main(),
    error=function(error) {
        message('Error: ', conditionMessage(error))
        quit(save='no', status=1L)
    }
)
