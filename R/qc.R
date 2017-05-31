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
        
        # Calculate some information
        info_binwidth = median(e$end - e$start)
        info_reads_per_bin = median(e$w + e$c)
        info_chrom_sizes = e[, .(xend = max(end)), by = chrom]
        info_num_bins = nrow(e)
        info_y_limit = 2*info_reads_per_bin+1
        info_sample_name = substr(s,1,25)
        if (nchar(s)>25) info_sample_name = paste0(info_sample_name, "...")
        info_cell_name   = substr(ce,1,25)
        if (nchar(ce)>25) info_cell_name = paste0(info_cell_name, "...")

        # start main plot:
        plt <- ggplot(e) +
            aes(x = (start+end)/2)

        # background rectangles
        if ("class" %in% colnames(e)) {

            # prepare consecutive rectangles for a better plotting experience
            consecutive = cumsum(abs(diff(as.numeric(as.factor(e$class)))))
            consecutive = c(consecutive[1], consecutive)
            e$consecutive = consecutive
            f = e[, .(start = min(start), end = max(end), class = class[1]), by = .(consecutive, chrom)][]

            plt <- plt +
                geom_rect(data = f, aes(xmin = start, xmax=end, ymin=-Inf, ymax=Inf, fill=class), inherit.aes=F, alpha=0.2) +
                scale_fill_manual(values = c(WW = "sandybrown", CC = "paleturquoise4", WC = "yellow", None = NA))
        }

        # Watson/Crick bars
        plt <- plt +
            geom_bar(aes(y = -w, width=(end-start)), stat='identity', position = 'identity', fill='sandybrown') +
            geom_bar(aes(y = c, width=(end-start)), stat='identity', position = 'identity', fill='paleturquoise4') +
            # Trim image to 2*median cov
            coord_flip(expand = F, ylim=c(-info_y_limit, info_y_limit)) +
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
            geom_segment(data = info_chrom_sizes, aes(xend = xend, x=0, y=-info_reads_per_bin, yend=-info_reads_per_bin),
                         linetype="dotted", col="darkgrey", size=0.5) +
            geom_segment(data = info_chrom_sizes, aes(xend = xend, x=0, y=+info_reads_per_bin, yend=+info_reads_per_bin),
                         linetype="dotted", col="darkgrey", size=0.5) +
			geom_segment(data = info_chrom_sizes, aes(xend = xend, x=0), y=0, yend=0, size=0.5)
    


        ###
        # Behaviour when strand states are given:
        # Plot a separate histogram for each state.
        if ("class" %in% colnames(e)) {

            e = merge(e, e[,.N,by=class], by = "class")
            e[, class := paste0(class, " (n=", N, ")")][]
            e.melt <- melt.data.table(e, c("chrom","start","end","class"), measure.vars = c("w","c"), variable.name = "strand", value.name = "coverage")
            plt_hist <- ggplot(e.melt) +
                aes(coverage, fill = strand) +
                geom_histogram(binwidth=1) +
                scale_x_continuous(limits = c(-1, 10+3*info_reads_per_bin), breaks = pretty_breaks(5), labels = comma) +
                theme(text = element_text(size=10), axis.text = element_text(size=8), axis.title.y = element_blank()) +
                      #axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
                scale_fill_manual(values = c(w = "sandybrown", c = "paleturquoise4")) +
                guides(fill=FALSE) +
                xlab("bin coverage") + facet_wrap(~class, nrow=1, scales = "free")


            all <- ggdraw() + draw_plot(plt) +
                draw_plot(plt_hist,                            x=.5, y=.76,  width=.49, height=.23) +
                draw_label(paste("Sample:", info_sample_name), x=.3,  y=.97, vjust=1, hjust=0, size=14) +
                draw_label(paste("Cell:", info_cell_name),     x=.3,  y=.94, vjust=1, hjust=0, size=12) +
                draw_label(paste("Median binwidth:", format(round(info_binwidth/1000,0), big.mark=",", scientific=F),"kb"),
                                                               x=.3,  y=.91, vjust=1, hjust=0, size=10) +
                draw_label(paste("Number bins:", format(info_num_bins, big.mark=",")),
                                                               x=.3,  y=.89, vjust=1, hjust=0, size=10) +
                draw_label(paste("Median reads/bin (dotted):", info_reads_per_bin),
                                                               x=.3,  y=.87, vjust=1, hjust=0, size=10) +
                draw_label(paste0("Plot limits: [-", info_y_limit, ",", info_y_limit, "]"),
                                                               x=.3,  y=.85, vjust=1, hjust=0, size=10)

        ###
        # Behaviour when strand states are not given:
        # Plot a coverage histogram and some information
        } else {

            # Coverage histogram:
            e.melt <- melt.data.table(e, c("chrom","start","end"), measure.vars = c("w","c"), variable.name = "strand", value.name = "coverage")
            plt_hist <- ggplot(e.melt) +
                aes(coverage, fill = strand) +
                geom_histogram(binwidth=1) +
                scale_x_continuous(limits = c(-1, 10+3*info_reads_per_bin), breaks = pretty_breaks(5), labels = comma) +
                theme(text = element_text(size=10), axis.text = element_text(size=8), axis.title.y = element_blank()) +
                      #axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
                scale_fill_manual(values = c(w = "sandybrown", c = "paleturquoise4")) +
                guides(fill=FALSE) +
                xlab("bin coverage")


            all <- ggdraw() + draw_plot(plt) +
                draw_plot(plt_hist,                            x=.5, y=.76,  width=.18, height=.23) +
                draw_label(paste("Sample:", info_sample_name), x=.3,  y=.97, vjust=1, hjust=0, size=14) +
                draw_label(paste("Cell:", info_cell_name),     x=.3,  y=.94, vjust=1, hjust=0, size=12) +
                draw_label(paste("Median binwidth:", format(round(info_binwidth/1000,0), big.mark=",", scientific=F),"kb"),
                                                               x=.3,  y=.91, vjust=1, hjust=0, size=10) +
                draw_label(paste("Number bins:", format(info_num_bins, big.mark=",")),
                                                               x=.3,  y=.89, vjust=1, hjust=0, size=10) +
                draw_label(paste("Median reads/bin (dotted):", info_reads_per_bin),
                                                               x=.3,  y=.87, vjust=1, hjust=0, size=10) +
                draw_label(paste0("Plot limits: [-", info_y_limit, ",", info_y_limit, "]"),
                                                               x=.3,  y=.85, vjust=1, hjust=0, size=10)
        }
        print(all)
    }
}
dev.off()
