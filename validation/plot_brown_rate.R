arguments <- commandArgs(trailingOnly=TRUE)
run_label <- if (length(arguments)) arguments[[1L]] else 'full'
output_directory <- file.path('validation', 'output')
curves <- utils::read.csv(file.path(
    output_directory,
    paste0('brown_', run_label, '_rate_irfs.csv')
))
population <- curves[curves$curve == 'population', ]
subjects <- curves[curves$curve == 'subject total', ]
items <- curves[curves$curve == 'item total', ]
subject_levels <- unique(subjects$group)
item_levels <- unique(items$group)
subject_colors <- grDevices::hcl.colors(length(subject_levels), 'Dynamic')

grDevices::png(
    file.path(output_directory, paste0('brown_', run_label, '_rate_irfs.png')),
    width=1100,
    height=1000
)
graphics::par(mfrow=c(2, 1), mar=c(4.5, 4.5, 3, 1))
graphics::plot(
    population$lag,
    population$estimate,
    type='l',
    lwd=4,
    xlab='Delay (seconds)',
    ylab='Rate IRF',
    ylim=range(curves$total_estimate),
    main=paste('Brown', run_label, 'rate IRF totals')
)
if (all(is.finite(population$se))) {
    graphics::polygon(
        c(population$lag, rev(population$lag)),
        c(
            population$estimate - 1.96 * population$se,
            rev(population$estimate + 1.96 * population$se)
        ),
        col=grDevices::adjustcolor('grey60', alpha.f=0.3),
        border=NA
    )
}
if (length(item_levels)) for (item in item_levels) {
    rows <- items$group == item
    graphics::lines(
        items$lag[rows],
        items$total_estimate[rows],
        col=grDevices::adjustcolor('grey50', alpha.f=0.35)
    )
}
for (i in seq_along(subject_levels)) {
    rows <- subjects$group == subject_levels[[i]]
    graphics::lines(
        subjects$lag[rows],
        subjects$total_estimate[rows],
        col=grDevices::adjustcolor(subject_colors[[i]], alpha.f=0.45)
    )
}
graphics::lines(population$lag, population$estimate, lwd=4)
legend_labels <- c('population', 'subject totals')
legend_colors <- c('black', subject_colors[[1L]])
if (length(item_levels)) {
    legend_labels <- c(legend_labels, 'sampled item totals')
    legend_colors <- c(legend_colors, 'grey50')
}
graphics::legend(
    'topright',
    legend_labels,
    col=legend_colors,
    lwd=rep.int(1, length(legend_labels)) + c(3, rep.int(0, length(legend_labels) - 1L)),
    bty='n'
)

deviations <- subjects$estimate
if (nrow(items)) deviations <- c(deviations, items$estimate)
graphics::plot(
    range(population$lag),
    range(deviations),
    type='n',
    xlab='Delay (seconds)',
    ylab='Random rate deviation',
    main='Grouped deviations (separate vertical scale)'
)
graphics::abline(h=0, lty=3, col='grey40')
if (length(item_levels)) for (item in item_levels) {
    rows <- items$group == item
    graphics::lines(
        items$lag[rows],
        items$estimate[rows],
        col=grDevices::adjustcolor('grey50', alpha.f=0.5)
    )
}
for (i in seq_along(subject_levels)) {
    rows <- subjects$group == subject_levels[[i]]
    graphics::lines(
        subjects$lag[rows],
        subjects$estimate[rows],
        col=grDevices::adjustcolor(subject_colors[[i]], alpha.f=0.55)
    )
}
grDevices::dev.off()
