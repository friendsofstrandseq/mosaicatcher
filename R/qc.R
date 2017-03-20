library(dplyr)
library(data.table)
library(assertthat)
library(ggplot2)
library(scales)
library(cowplot)

args <- commandArgs(trailingOnly = T)
f_in <- args[1]
pdf_out <- args[2]

format_Mb <- function(x) {
    paste(comma(x/1e6), "Mb")
}
    # MAIN

# Read counts & filter chromosomes
d = fread(f_in)
d = d[grepl('^chr[0-9XY]+$', chrom),]
d = d[, chrom := factor(chrom, levels=paste0('chr', c(1:22,'X','Y')), ordered = T)]

perSample <- d %>% group_by(sample) %>%
    summarize(median = median(w+c)) %>%
    as.data.table


# Plot all cells
cairo_pdf(pdf_out, width=14, height=10, onefile = T)
for (s in unique(d$sample))
{
    message("Plotting ", s)
    e = d[sample == s,]
    binwidth = e$end[1] - e$start[1]
    reads_per_bin = perSample[sample == s,]$median
    
    # main plot:
    plt <- ggplot(e) +
        aes(x = (start+end)/2) +
        geom_bar(aes(y = -w), stat='identity', position = 'identity', fill='darkorange', width=binwidth) +
        geom_bar(aes(y = c), stat='identity', position = 'identity', fill='dodgerblue3', width=binwidth) +
        coord_flip(expand = F, ylim=c(-2*reads_per_bin,2*reads_per_bin)) +
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

    # add vertical bars at +/- median:
    chrom_sizes = group_by(e, chrom) %>% summarize(xend = max(end))
    plt <- plt + geom_segment(data=chrom_sizes, aes(xend = xend, x=0, y=-reads_per_bin, yend=-reads_per_bin),  linetype="dotted", col="darkgrey") +
                 geom_segment(data=chrom_sizes, aes(xend = xend, x=0, y=+reads_per_bin, yend=+reads_per_bin),  linetype="dotted", col="darkgrey")

    # Histogram:
    e.melt <- melt.data.table(e, c("chrom","start","end"), measure.vars = c("w","c"), variable.name = "strand", value.name = "coverage")
    plt_hist <- ggplot(e.melt) + 
        aes(coverage, fill = strand, col = strand, y = ..scaled..) +
        geom_density(adjust=0.5, alpha=0.3) +
        scale_x_continuous(limits = c(0,3*reads_per_bin), breaks = pretty_breaks(4), labels = comma)
    
    s_name = substr(s,1,25)
    if (nchar(s)>25) s_name = paste0(s_name, "...")
    
    all <- ggdraw() + draw_plot(plt) + 
        draw_plot(plt_hist, .69, .77, .3,.22) +
        draw_label(paste("Sample:", s_name), x=.2, y=.97, vjust=1, hjust=0, size = 14) +
        draw_label(paste("Binwidth:", format(binwidth, big.mark=",")), x=.2, y=.94, vjust=1, hjust=0, size=10) +
        draw_label(paste("Reads/Bin:", reads_per_bin), x=.2, y=.92, vjust=1, hjust=0, size=10)
    print(all)
}
dev.off()


