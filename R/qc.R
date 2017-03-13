library(dplyr)
library(data.table)
library(assertthat)
library(ggplot2)
library(scales)
library(cowplot)

args <- commandArgs(trailingOnly = T)
f_in <- args[1]
pdf_out <- args[2]


# FUNCTIONS

norm_per_sample <- function(dt, minSum = 0) {
    perSample <- dt %>% group_by(sample) %>%
        summarize(sum = sum(w+c), mean = mean(w+c),median = median(w+c)) %>%
        as.data.table
    plt <- ggplot(perSample) + aes(sum) + geom_histogram(binwidth=1e5) +
        ggtitle("total reads per cell") + 
        geom_vline(xintercept = minSum, col="red", linetype="dashed") +
        scale_x_continuous(label=comma) + xlab("number of reads")
    goodSamples <- perSample[sum >= minSum, ]
    if (nrow(goodSamples) < nrow(perSample)) {
        diff <- nrow(perSample) - nrow(goodSamples)
        message("Removing ", diff, " cells because their total coverage is too low")
    }
    d <- merge(dt, goodSamples, by=c("sample"),all = F)
    d[, wn:= w/median]
    d[, cn:= c/median]
    d[, sum := NULL]
    d[, mean := NULL]
    d[, median := NULL]
    return (list("plot" = plt, "data" = d))
}

format_Mb <- function(x) {
    paste(comma(x/1e6), "Mb")
}
    # MAIN

# Read counts & filter chromosomes
d = fread(f_in)

# Filter/modify data:
d = d[grepl('^chr[0-9XY]+$', chrom),]
d[, chrom := factor(chrom, levels=paste0('chr', c(1:22,'X','Y')), ordered = T)]
d = norm_per_sample(d)$data


# Plot all cells
cairo_pdf(pdf_out, width=14, height=10, onefile = T)
for (s in unique(d$sample))
{

    message("Plotting ", s)
    e = d[sample == s,]
    binwidth = e$end[1] - e$start[1]
    n_bins_removed = "??"
    reads_per_bin = round(sum(c(e$w, e$c)) / nrow(e), 1)
    
    # main plot:
    plt <- ggplot(e) +
        aes(x = (start+end)/2) +
        geom_rect(aes(xmin=start, xmax=end, ymin = -Inf, ymax = Inf, fill = class), alpha=0.2) +
        geom_bar(aes(y = -wn), stat='identity', position = 'identity', fill='darkorange', width=binwidth) +
        geom_bar(aes(y = cn), stat='identity', position = 'identity', fill='dodgerblue3', width=binwidth) +
        coord_flip(expand = F, ylim=c(-2,2)) +
        facet_grid(.~chrom, switch="x") +
        ylab("Watson | Crick") +
        xlab(NULL) +
        scale_x_continuous(breaks = pretty_breaks(12), labels = format_Mb) +
        scale_y_continuous(breaks = pretty_breaks(3)) + 
        #scale_fill_manual(values = c("WC" = "darkolivegreen2", "CC" = "lightskyblue1", "WW" = "orange")) +
        theme_classic() +
        theme(panel.margin = unit(0, "lines"),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.background = element_rect(fill = NA, colour=NA))

    # add vertical bars at +/- 1:
    chrom_sizes = group_by(e, chrom) %>% summarize(xend = max(end))
    plt <- plt + geom_segment(data=chrom_sizes, aes(xend = xend, x=0, y=-1, yend=-1),  linetype="dotted", col="darkgrey") +
                 geom_segment(data=chrom_sizes, aes(xend = xend, x=0, y=+1, yend=+1),  linetype="dotted", col="darkgrey")

    # Histogram:
    e.melt <- melt.data.table(e, "class", measure.vars = c("wn","cn"), variable.name = "strand", value.name = "coverage")
    e.classes = table(e.melt$class)
    e.melt$class = factor(paste0(e.melt$class, " (n=", e.classes[e.melt$class], ")"))
    plt_hist <- ggplot(e.melt) + aes(coverage, fill = strand, col = strand,y = ..scaled..) +
        geom_density(adjust=0.5, alpha=0.3) +
        guides(color = FALSE, fill = FALSE) +
        scale_x_continuous(limits = c(0,2)) +
        facet_wrap(~class, nrow=1, scales = "free_y")
    
    s_name = substr(s,1,25)
    if (nchar(s)>25) s_name = paste0(s_name, "...")

    background_rate = sum(c(e[class=="WW",]$cn, e[class=="CC",]$wn)) / sum(e[class=="WW" | class=="CC"]$cn + e[class=="WW" | class=="CC"]$wn)
    
    all <- ggdraw() + draw_plot(plt) + 
        draw_plot(plt_hist, .39, .73, .6,.22) +
        draw_label(paste("Sample:", s_name), x=.2, y=.97, vjust=1, hjust=0, size = 14) +
        draw_label(paste("Binwidth:", format(binwidth, big.mark=",")), x=.2, y=.94, vjust=1, hjust=0, size=10) +
        draw_label(paste("Reads/Bin:", reads_per_bin), x=.2, y=.91, vjust=1, hjust=0, size=10) +
        draw_label(paste("Background rate: ", round(background_rate*100,2), "%"), x=.2, y=.88, vjust=1, hjust=0, size=10)
    print(all)
}
dev.off()


