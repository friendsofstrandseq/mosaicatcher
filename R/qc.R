library(dplyr)
library(data.table)
library(assertthat)
library(ggplot2)
library(scales)
library(cowplot)

#args <- commandArgs(trailingOnly = T)
args <- c("data/HG00733.200kb.txt", "data/qc.pdf")
f_in <- args[1]
pdf_out <- args[2]


# FUNCTIONS

format_Mb <- function(x) {
    paste(comma(x/1e6), "Mb")
}

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
    d[, w:= w/median]
    d[, c:= c/median]
    d[, sum := NULL]
    d[, mean := NULL]
    d[, median := NULL]
    return (list("plot" = plt, "data" = d))
}

norm_per_bin <- function(dt, minMedian = 1e-5, maxMedian = 5, maxSD = Inf) {
    perBin <- dt %>% group_by(chrom, start, end) %>%
        summarize(mean = mean(w+c), median = median(w+c), sd = sd(w+c)) %>%
        as.data.table
    plt <- ggplot(perBin) + aes(median) + geom_histogram(binwidth=1e5) +
        ggtitle("median norm. cov per bin") + 
        geom_vline(xintercept = c(minMedian,maxMedian), col="red", linetype="dashed") +
        scale_x_continuous(label=comma) + xlab("number of reads per sample")

    goodBins <- perBin[median >= minMedian & sd <= maxSD,]
    diff = nrow(perBin) - nrow(goodBins)
    if (diff > 0) {
        message("Removing ", diff, " bins")
    }
    d <- merge(dt, goodBins, all=T, by=c("chrom","start","end"))
    d[, w:= w/median]
    d[, c:= c/median]
    d[, sd := NULL]
    d[, mean := NULL]
    d[, median := NULL]
    return (list("plot" = plt, "data" = d, "n_bins_removed" = diff))
}


print_single_sample <- function(d, samplename = NULL) 
{
    assert_that(is.data.table(d))
    assert_that(all(c("start","end","chrom","sample","w","c") %in% colnames(d)))
    assert_that(is.character(samplename) || is.null(samplename))
    
    if (is.null(samplename) || !(samplename %in% unique(d$sample))) {
        samplename <- unique(d_s$sample)[1]
        message("Plot for sample ", samplename)
    }
    e <- d[sample == samplename,]
    
    # Remove outliers
    #bad_idx <- with(e, w+c > mean(w+c, na.rm=T) | w+c < mean(w+c, na.rm=T))
    #if (n2 < n1) message("Removed ", n1-n2, " (out of ", n1, ") bins with too high coverage")
    
    
    # To do: What if whole chromosome drops out ? --> Faceting error
    binwidth = e$end[1] - e$start[1] - 1
    plt <- ggplot(e[!is.na(w),]) + 
        aes(x = (start+end)/2) +
        geom_bar(aes(y = -w), stat='identity', position = 'identity', fill='darkorange', width=binwidth) + 
        geom_bar(aes(y = c), stat='identity', position = 'identity', fill='dodgerblue3', width=binwidth) +
        geom_rect(data = e[is.na(w),], aes(xmin=start, xmax=end, ymin = -Inf, ymax = Inf), fill='grey', alpha=0.2) +
        coord_flip(expand = F) + 
        facet_grid(.~chrom) +
        ggtitle(samplename) +
        ylab("Watson | Crick") +
        xlab(NULL) +
        scale_x_continuous(breaks = pretty_breaks(8), labels = format_Mb) +
        scale_y_continuous(breaks = pretty_breaks(3))
    return(plt)
}





# MAIN

# Read counts & filter chromosomes
d = fread(f_in)
d = d[grepl('^chr[0-9X]+$', chrom),]
d[, chrom := factor(chrom, levels=paste0('chr', c(1:22,'X')), ordered = T)]

# Normalize per sample
nps <- norm_per_sample(d, 3e4)

# Normalize by bin
npb <- norm_per_bin(nps$data, minMedian = 0.05, maxMedian = 2.5)

# Plot all cells
cairo_pdf(pdf_out, width=14, height=10, onefile = T)
for (s in unique(npb$data$sample)) 
{
    message("Plotting ", s)
    norm = npb$data[sample == s,]
    orig = d[sample == s,]
    binwidth = norm$end[1] - norm$start[1]
    n_bins_removed = npb$n_bins_removed
    reads_per_bin = round(sum(c(orig$w, orig$c)) / nrow(orig), 1)
    
    # main plot:
    plt <- ggplot(norm[!is.na(w),]) + 
        aes(x = (start+end)/2) +
        geom_bar(aes(y = -w), stat='identity', position = 'identity', fill='darkorange', width=binwidth) + 
        geom_bar(aes(y = c), stat='identity', position = 'identity', fill='dodgerblue3', width=binwidth) +
        geom_rect(data = norm[is.na(w),], aes(xmin=start, xmax=end, ymin = -Inf, ymax = Inf), fill='grey', alpha=0.2) +
        coord_flip(expand = F) + 
        facet_grid(.~chrom, switch="x") +
        ylab("Watson | Crick") +
        xlab(NULL) +
        scale_x_continuous(breaks = pretty_breaks(8), labels = format_Mb) +
        scale_y_continuous(breaks = pretty_breaks(3)) + 
        theme_classic() +
        theme(panel.margin = unit(0, "lines"),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.background = element_rect(fill = NA, colour=NA))
    
    # auxillaries
    plt_h_norm <- ggplot(norm[!is.na(w),]) + aes(w+c) + 
        geom_histogram(binwidth=0.05) + 
        xlab("bin coverage (w+c)") + 
        ggtitle("normalized coverage per bin") +
        theme_classic() + 
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
    plt_h_orig <- plt_h_norm %+% orig + geom_histogram(binwidth=5) +
        ggtitle("raw coverage per bin")
    
    all <- ggdraw() + draw_plot(plt) + 
        draw_plot(plt_h_orig, .4, .68, .3,.3) +
        draw_plot(plt_h_norm, .7, .68, .3,.3) +
        draw_label(paste("Sample:", s), x=.2, y=.97, vjust=1, hjust=0, size = 14) +
        draw_label(paste("Binwidth:", format(binwidth, big.mark=",")), x=.2, y=.92, vjust=1, hjust=0, size=10) +
        draw_label(paste("Removed bins:", n_bins_removed), x=.2, y=.88, vjust=1, hjust=0, size=10) +
        draw_label(paste("Reads/Bin:", reads_per_bin), x=.2, y=.84, vjust=1, hjust=0, size=10)
    print(all)
}    
dev.off()

