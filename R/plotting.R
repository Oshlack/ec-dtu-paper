plot_bottomly_boxplot <- function(results, cols, title='', toplot=c('FDR', 'TPR'), hline=NA) {
    results$feature <- factor(results$feature, levels=names(cols))
    p <- ggplot(results, aes_string('feature', toplot, group='feature', colour='feature')) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width=0.2, size=0.4) +
            theme_bw() +
            theme(legend.position = 'none',
                  plot.title = element_text(hjust=0.5),
                  text = element_text(size = 18)) +
            ylim(0,1) +
            ylab('') +
            xlab('') +
            ggtitle(title) +
            geom_hline(yintercept = 0.05, colour='grey',  linetype='dotted') +
            scale_color_manual(values = cols)
    if(!is.na(hline)) {
        p <- p + geom_hline(yintercept = hline, colour='grey',  linetype='dotted')
    }
    return(p)
}
