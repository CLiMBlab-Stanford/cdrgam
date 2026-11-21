arguments <- commandArgs(trailingOnly=TRUE)
if (length(arguments) != 1L) {
    stop('usage: check-release-tag.R TAG', call.=FALSE)
}

git_output <- function(...) {
    arguments <- vapply(list(...), shQuote, character(1L))
    output <- system2('git', arguments, stdout=TRUE, stderr=TRUE)
    if (!is.null(attr(output, 'status'))) {
        stop(paste(output, collapse='\n'), call.=FALSE)
    }
    paste(output, collapse='\n')
}

tag <- arguments[[1L]]
if (!grepl('^v[0-9]+\\.[0-9]+\\.[0-9]+$', tag)) {
    stop('release tags must have the form vMAJOR.MINOR.PATCH', call.=FALSE)
}
package_version <- unname(read.dcf('DESCRIPTION', fields='Version')[1L, 1L])
if (!identical(substring(tag, 2L), package_version)) {
    stop('tag ', tag, ' does not match package version ', package_version,
        call.=FALSE)
}
if (!identical(git_output('cat-file', '-t', paste0('refs/tags/', tag)), 'tag')) {
    stop('release tag ', tag, ' must be annotated', call.=FALSE)
}
annotation <- trimws(git_output(
    'for-each-ref', '--format=%(contents)', paste0('refs/tags/', tag)
))
if (nchar(annotation) < 20L) {
    stop('release tag ', tag, ' must include a meaningful release summary',
        call.=FALSE)
}
head <- git_output('rev-parse', 'HEAD')
target <- git_output('rev-parse', paste0('refs/tags/', tag, '^{commit}'))
if (!identical(target, head)) {
    stop('release tag ', tag, ' does not identify the checked-out commit',
        call.=FALSE)
}
ancestry <- system2(
    'git', c('merge-base', '--is-ancestor', head, 'refs/remotes/origin/main'),
    stdout=FALSE, stderr=FALSE
)
if (!identical(ancestry, 0L)) {
    stop('release tag ', tag, ' does not identify a commit on main',
        call.=FALSE)
}
cat('validated release tag', tag, 'at', head, '\n')
