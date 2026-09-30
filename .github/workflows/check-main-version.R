arguments <- commandArgs(trailingOnly=TRUE)
if (length(arguments) != 1L) {
    stop('usage: check-main-version.R BASE_COMMIT', call.=FALSE)
}

read_version <- function(lines) {
    connection <- textConnection(lines)
    on.exit(close(connection))
    unname(read.dcf(connection, fields='Version')[1L, 1L])
}

release_version <- function(version) {
    if (!grepl('^[0-9]+\\.[0-9]+\\.[0-9]+$', version)) {
        stop('release versions must have the form MAJOR.MINOR.PATCH: ', version,
            call.=FALSE)
    }
    numeric_version(version)
}

base_description <- system2(
    'git', c('show', paste0(arguments[[1L]], ':DESCRIPTION')),
    stdout=TRUE, stderr=TRUE
)
if (!is.null(attr(base_description, 'status'))) {
    stop(paste(base_description, collapse='\n'), call.=FALSE)
}

previous <- read_version(base_description)
proposed <- read_version(readLines('DESCRIPTION', warn=FALSE))
if (release_version(proposed) <= release_version(previous)) {
    stop('main release version must advance: ', previous, ' -> ', proposed,
        call.=FALSE)
}
cat('main release version advances:', previous, '->', proposed, '\n')
