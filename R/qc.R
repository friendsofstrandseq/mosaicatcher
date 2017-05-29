library(data.table)
library(ggplot2)
library(scales)
library(cowplot)

args <- commandArgs(trailingOnly = T)
if (length(args) != 2 | !grepl('\\.pdf$', args[2])) {
    warning("Usage: Rscript R/qc.R input-file output-pdf")
    quit(status = 1)
}
f_in <- args[1]
pdf_out <- args[2]

format_Mb <- function(x) {
    paste(comma(x/1e6), "Mb")
}

# Read counts & filter chromosomes (this is human-specific)
d = fread(f_in)
d = d[, chrom := sub('^chr','',chrom)][]
d = d[grepl('^([1-9]|[12][0-9]|X|Y)$', chrom),]
d = d[, chrom := factor(chrom, levels=as.character(c(1:22,'X','Y')), ordered = T)]

# Plot all cells
cairo_pdf(pdf_out, width=14, height=10, onefile = T)
for (s in unique(d$sample))
{
    for (ce in unique(d[sample == s,]$cell))
    {
        message(paste("Plotting sample", s, "cell", ce,"into",pdf_out))
        e = d[sample == s & cell == ce,]
        
        e$rand_bg = rep(c("a","b"), nrow(e))[1:nrow(e)]
        binwidth = median(e$end - e$start)
        reads_per_bin = median(e$w + e$c)
        chrom_sizes = e[, .(xend = max(end)), by = chrom]
        num_bins = nrow(e)
		y_limit = 2*reads_per_bin+1
    
        # main plot:
        plt <- ggplot(e) +
            aes(x = (start+end)/2) +
            # background rectangles so you see how large the bins are
            #geom_rect(aes(xmin = start, xmax=end, ymin=-2*reads_per_bin, ymax=2*reads_per_bin, fill=rand_bg), inherit.aes=F, alpha=0.2) +
            # Watson/Crick bars
            geom_bar(aes(y = -w, width=(end-start)), stat='identity', position = 'identity', fill='sandybrown') +
            geom_bar(aes(y = c, width=(end-start)), stat='identity', position = 'identity', fill='paleturquoise4') +
            # Trim image to 2*median cov
            coord_flip(expand = F, ylim=c(-y_limit, y_limit)) +
            facet_grid(.~chrom, switch="x") +
            ylab("Watson | Crick") + xlab(NULL) +
            scale_x_continuous(breaks = pretty_breaks(12), labels = format_Mb) +
            scale_y_continuous(breaks = pretty_breaks(3)) + 
            theme_classic() +
            theme(panel.margin = unit(0, "lines"),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  strip.background = element_rect(fill = NA, colour=NA)) + 
            guides(fill = FALSE) +
            # Dotted lines at median bin count
            geom_segment(data=chrom_sizes, aes(xend = xend, x=0, y=-reads_per_bin, yend=-reads_per_bin),  linetype="dotted", col="darkgrey") +
            geom_segment(data=chrom_sizes, aes(xend = xend, x=0, y=+reads_per_bin, yend=+reads_per_bin),  linetype="dotted", col="darkgrey") +
			geom_segment(data = chrom_sizes, aes(xend = xend, x=0), y=0, yend=0, size=0.5)
    
        # Histogram:
        e.melt <- melt.data.table(e, c("chrom","start","end"), measure.vars = c("w","c"), variable.name = "strand", value.name = "coverage")
        plt_hist <- ggplot(e.melt) + 
            aes(coverage, fill = strand, col = strand, y = ..scaled..) +
            geom_density(alpha=0.3) +
            scale_x_continuous(limits = c(0,3*reads_per_bin), breaks = pretty_breaks(4), labels = comma) +
            theme(axis.title.y = element_blank())
    
        # count histogram
        plt_covhist <- ggplot(e) + aes(w+c) + geom_histogram(bins = 30) +
            xlab("coverage") + theme(axis.title.y = element_blank()) +
            scale_x_continuous(limits = c(0, 5*reads_per_bin), label=comma, breaks = pretty_breaks(3))
        
        s_name = substr(s,1,25)
        if (nchar(s)>25) s_name = paste0(s_name, "...")
        
        all <- ggdraw() + draw_plot(plt) +
		    draw_plot(plt_hist,    .74, .75, .25, .2) +
            draw_plot(plt_covhist, .56, .75, .18, .2) +
            draw_label(paste("Sample:", s_name), x=.4, y=.97, vjust=1, hjust=0, size = 14) +
            draw_label(paste("Median binwidth:", format(round(binwidth/1000,0), big.mark=",", scientific=F),"kb"), x=.4, y=.94, vjust=1, hjust=0, size=10) +
            draw_label(paste("Number bins:", format(num_bins, big.mark=",")), x=.4, y=.92, vjust=1, hjust=0, size=10) +
            draw_label(paste("Median reads/bin (dotted):", reads_per_bin), x=.4, y=.90, vjust=1, hjust=0, size=10) +
			draw_label(paste0("Plot limits: [-", y_limit, ",", y_limit, "]"), x=.4, y=.88, vjust=1, hjust=0, size=10)
        print(all)
    }
}
dev.off()
