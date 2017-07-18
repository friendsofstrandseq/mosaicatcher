library(data.table)
library(assertthat)
library(ggplot2)
library(scales)
library(cowplot)

args <- commandArgs(trailingOnly = T)
if (length(args) < 2 | length(args) > 3 | !grepl('\\.pdf$', args[length(args)])) {
    warning("Usage: Rscript R/qc.R input-file [cell-info-file] output-pdf")
    quit(status = 1)
}
f_in <- args[1]
pdf_out <- args[length(args)]
if(length(args)==3) f_info = args[2] else f_info = NULL

format_Mb <- function(x) {
    paste(comma(x/1e6), "Mb")
}

# if gzip
if (substr(f_in,nchar(f_in)-2,nchar(f_in)) == ".gz")
    f_in = paste("zcat",f_in)

# Read counts & filter chromosomes (this is human-specific)
d = fread(f_in)

# Check that correct files are given:
assert_that("chrom" %in% colnames(d))
assert_that("start" %in% colnames(d) && is.integer(d$start))
assert_that("end" %in% colnames(d) && is.integer(d$end))
assert_that("sample" %in% colnames(d))
assert_that("cell" %in% colnames(d))
assert_that("w" %in% colnames(d) && is.integer(d$w))
assert_that("c" %in% colnames(d) && is.integer(d$c))
assert_that("class" %in% colnames(d))

# Re-name and -order chromosomes
d = d[, chrom := sub('^chr','',chrom)][]
d = d[grepl('^([1-9]|[12][0-9]|X|Y)$', chrom),]
d = d[, chrom := factor(chrom, levels=as.character(c(1:22,'X','Y')), ordered = T)]


# Read cell info file for bells and whistles
if (!is.null(f_info)) {
    I = fread(f_info, skip=12)
    assert_that("sample" %in% colnames(I))
    assert_that("cell" %in% colnames(I))
    assert_that("pass1" %in% colnames(I))
    assert_that("dupl" %in% colnames(I))
    assert_that("mapped" %in% colnames(I))
    assert_that("nb_p" %in% colnames(I) && is.numeric(I$nb_p))
    assert_that("nb_n" %in% colnames(I) && is.numeric(I$nb_n))
    assert_that("nb_z" %in% colnames(I) && is.numeric(I$nb_z))
}


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
        info_total_reads = sum(e$c + e$w)
        info_y_limit = 2*info_reads_per_bin+1
        info_sample_name = substr(s,1,25)
        if (nchar(s)>25) info_sample_name = paste0(info_sample_name, "...")
        info_cell_name   = substr(ce,1,25)
        if (nchar(ce)>25) info_cell_name = paste0(info_cell_name, "...")

        # start main plot:
        plt <- ggplot(e) +
            aes(x = (start+end)/2)

    
        # prepare consecutive rectangles for a better plotting experience
        consecutive = cumsum(c(0,abs(diff(as.numeric(as.factor(e$class))))))
        e$consecutive = consecutive
        f = e[, .(start = min(start), end = max(end), class = class[1]), by = .(consecutive, chrom)][]

        plt <- plt +
            geom_rect(data = f, aes(xmin = start, xmax=end, ymin=-Inf, ymax=Inf, fill=class), inherit.aes=F, alpha=0.2) +
            scale_fill_manual(values = c(WW = "sandybrown", CC = "paleturquoise4", WC = "yellow", None = NA))

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
        
        # Rename classes:
        labels = e[,.N,by = class][,label := paste0(class," (n=", N, ")")][]
        
        e[, class := factor(class, levels = labels$class, labels = labels$label)]
        
        # Histogram in upper right corner
        e.melt <- melt.data.table(e, c("chrom","start","end","class"), measure.vars = c("w","c"), variable.name = "strand", value.name = "coverage")
        plt_hist_xlim = 10+3*info_reads_per_bin
        plt_hist <- ggplot(e.melt) +
            aes(coverage, fill = strand) +
            geom_histogram(binwidth=1, position=position_dodge(), alpha=0.9) +
            scale_x_continuous(limits = c(-1, plt_hist_xlim), breaks = pretty_breaks(5), labels = comma) +
            theme(text = element_text(size=10), axis.text = element_text(size=8)) +
            scale_fill_manual(values = c(w = "sandybrown", c = "paleturquoise4")) +
            guides(fill=FALSE,col=FALSE) + ylab("bin count") +
            xlab("reads per bin") + facet_wrap(~class, nrow=1, scales = "free")

        if (!is.null(f_info)) {
            Ie = I[sample == s & cell == ce,]
            if(nrow(Ie) == 1 && Ie$pass1==1) {
                p = Ie$nb_p
                n = Ie$nb_n
                z = Ie$nb_z
                x = seq(0,plt_hist_xlim)
                scale_factors = e.melt[,.N, by = .(class,strand,coverage)][, .(scale = max(N)), by = class]
                nb = data.table(x      = rep(x,6),
                                strand = rep(c(rep("w",length(x)), rep("c",length(x))),3),
                                class  = c(rep(labels[class=="WW",]$label, 2*length(x)), 
                                           rep(labels[class=="WC",]$label ,2*length(x)), 
                                           rep(labels[class=="CC",]$label, 2*length(x))),
                                scale  = c(rep(scale_factors[class == labels[class=="WW",]$label,]$scale, 2*length(x)),
                                           rep(scale_factors[class == labels[class=="WC",]$label,]$scale, 2*length(x)),
                                           rep(scale_factors[class == labels[class=="CC",]$label,]$scale, 2*length(x))),
                                y      = c(dnbinom(x, z,   p),
                                           dnbinom(x, 2*n, p),
                                           dnbinom(x, n,   p),
                                           dnbinom(x, n,   p),
                                           dnbinom(x, 2*n, p),
                                           dnbinom(x, z,   p)  ))
                plt_hist <- plt_hist + geom_line(data = nb, aes(x,y*scale,col=strand))
            }
        }


        plot_hst_width = .03 + .13*length(unique(e$class))
        all <- ggdraw() + draw_plot(plt) +
            draw_plot(plt_hist,                            x=.45, y=.76,  width=plot_hst_width, height=.23) +
            draw_label(paste("Sample:", info_sample_name), x=.29, y=.97, vjust=1, hjust=0, size=14) +
            draw_label(paste("Cell:", info_cell_name),     x=.29, y=.94, vjust=1, hjust=0, size=12) +
            draw_label(paste("Median binwidth:", format(round(info_binwidth/1000,0), big.mark=",", scientific=F),"kb"),
                                                           x=.29, y=.91, vjust=1, hjust=0, size=10) +
            draw_label(paste("Number bins:", format(info_num_bins, big.mark=",")),
                                                           x=.29, y=.89, vjust=1, hjust=0, size=10) +
            draw_label(paste("Total number of reads:", format(info_total_reads, big.mark=",")),
                                                           x=.29, y=.87, vjust=1, hjust=0, size=10) +
            draw_label(paste("Median reads/bin (dotted):", info_reads_per_bin),
                                                           x=.29, y=.85, vjust=1, hjust=0, size=10) +
            draw_label(paste0("Plot limits: [-", info_y_limit, ",", info_y_limit, "]"),
                                                           x=.29, y=.83, vjust=1, hjust=0, size=10)
        
        # If available, add additional info like duplicate rate and NB params!
        if (!is.null(f_info)) {
            Ie = I[sample == s & cell == ce,]
            if(nrow(Ie) == 1) {
                all <- all +
                    draw_label(paste0("Duplicate rate: ", round(Ie$dupl/Ie$mapped,2)*100,"%"),
                               x=.29,  y=.80, vjust=1, hjust=0, size=10)
                if (Ie$pass1 == 1) {
                    all <- all +
                        draw_label(paste0("NB parameters (p,n,z): ", round(Ie$nb_p,2), ",",  round(Ie$nb_n,2), ",", round(Ie$nb_z,2)),    
                                   x=.29,  y=.78, vjust=1, hjust=0, size=10)
                }
            }
        }

        print(all)
    }
}
dev.off()
