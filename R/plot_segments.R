#
# Copyright (C) 2017 Sascha Meiers
# Distributed under the MIT software license, see the accompanying
# file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.

suppressMessages(library(data.table))
suppressMessages(library(scales))
suppressMessages(library(ggplot2))
suppressMessages(library(assertthat))
suppressMessages(library(cowplot))

print_usage <- function() {
    message("                                                                                ")
    message("Usage: Rscript plot_segments.R count-file segment-file chrom out.pdf            ")
    message("                                                                                ")
    message("       Plot different segmentations for a single chromosome.                    ")
    message("                                                                                ")
}


###########
# Arguments
args = commandArgs(trailingOnly = T)
if (length(args) != 4) {
    print_usage()
    stop("Stopped execution on purpose.")
}
      f_counts = args[1]
    f_segments = args[2]
         CHROM = args[3]
         f_out = args[4]
  zcat_command = "zcat"
     format_Mb = function(x) {paste(comma(round(x/1e6,1)), "Mb")}
cells_per_page = 10
invisible(assert_that(grepl('\\.pdf$', f_out)))


#############
# Read counts
message(" * Reading count data ", f_counts, "...")
if (grepl('\\.gz$',f_counts)) {
    counts = fread(paste(zcat_command, f_counts))
} else {
    counts = fread(f_counts)
}
invisible(assert_that(all(c("chrom","start","end","class","sample","cell","w","c") %in% colnames(counts))))
setkey(counts, chrom)


######################################
# subset to  chromosome and good cells
invisible(assert_that(CHROM %in% unique(counts$chrom)))
counts <- counts[chrom == CHROM]

#######################
# Calculate information
n_cells_orig = nrow(unique(counts[,.(sample,cell)]))
n_samples    = nrow(unique(counts[,.(sample)]))
good_cells = unique(counts[class != "None", .(sample, cell)])
counts <- counts[good_cells, on = .(sample,cell), nomatch = 0]
n_cells = nrow(unique(counts[,.(sample,cell)]))
chrom_size = unique(counts[,.(chrom, start, end)])
chrom_size = max(chrom_size$end) - min(chrom_size$start)


#########################
# Select randomly 8 cells
message(" * Randomly picking ", cells_per_page, " cells. If you're unhappy with this selection, just run the script again!")
good_cells <- good_cells[sample(1:nrow(good_cells),8)]
counts <- counts[good_cells, on = .(sample,cell), nomatch = 0]
y_lim = 3*counts[, median(w+c)]



######################
# Prepare None regions
counts[, cnsc := cumsum(c(0,abs(sign(diff(as.numeric(as.factor(class))))))), by = .(sample,cell)]
none_regions = counts[class == 'None'][, .(start = min(start), end = max(end)), by = .(sample,cell, class, cnsc)]


#########################
# Read segmentation table
message(" * Reading segmentation data ", f_segments, "...")
segs = fread(f_segments)
invisible(assert_that(is.data.table(segs),
            "chrom" %in% colnames(segs),
            "start" %in% colnames(segs),
            "end"   %in% colnames(segs),
            "k"     %in% colnames(segs),
            "sse"   %in% colnames(segs),
            "none_regions" %in% colnames(segs)))
segs = segs[chrom == CHROM]
segs[, SV_class := rep(paste0("seg",1:3), ceiling(.N/3))[1:.N], by = .(k)]

bg_cols = NULL
for (i in 1:nrow(good_cells)) {
    bg_cols <- rbind(bg_cols, cbind(segs, sample = good_cells$sample[i], cell = good_cells$cell[i]))
}




############
# Start plot
cairo_pdf(f_out, width=14, height=10, onefile = T)


################################
# Overview page at the beginning
if (TRUE) {

    message(" * Plotting overview page...")
    plt_sse <- ggplot(unique(segs[, .(k,sse)])) +
        aes(k,sse) +
        geom_line() +
        geom_point(shape = 4) +
        scale_y_log10(labels = comma, breaks = pretty_breaks(10)) +
        theme_gray() +
        ylab("Standard squared error") +
        xlab("Number of segments") +
        geom_vline(xintercept = segs$none_regions[1], linetype = "dashed", color = "dodgerblue") +
        geom_vline(xintercept = 2*segs$none_regions[1], linetype = "dashed", color = "darkorange") +
        annotate("text", x = segs$none_regions[1]-1, y = Inf, label = "none", vjust = 1, hjust = 1) +
        annotate("text", x = 1+2* segs$none_regions[1], y = Inf, label = "2 x none", vjust = 1, hjust = 0)
        #geom_text(data = NULL, label = "none", x = segs$none_regions[1], y = Inf, hjust = 0, vjust = 1, inherit.aes = F) +
        #geom_text(data = NULL, label = "2 x none", x = 2*segs$none_regions[1], y = Inf, hjust = 0, vjust = 1, inherit.aes = F)

    plt_title <- ggdraw() + draw_label(paste("Segmentation on chromosome", CHROM, "across", n_cells, "cells from", n_samples, "sample(s)"), fontface='bold')
    plt_content <- ggdraw() +
        draw_plot(plt_sse, x=.4, y=.1,  width=.55, height=.8) +
        draw_label("Segmentation summary",
                   x=.05, y=.9, vjust=1, hjust=0, size=14, fontface = "bold") +
        draw_label(paste("Chromosome:",CHROM),
                   x=.05, y=.85, vjust=1, hjust=0, size=12) +
        draw_label(paste0("Chrom size: ", format_Mb(chrom_size), ", ", segs$bins[1], " bins"),
                   x=.05, y=.82, vjust=1, hjust=0, size=12) +
        draw_label(paste0("Number of cells: ", n_cells, "/", n_cells_orig),
                   x=.05, y=.79, vjust=1, hjust=0, size=12) +

        draw_label(paste0("None regions: ", segs$none_regions[1], " regions, ", segs$none_bins[1], " bins"),
                   x=.05, y=.74, vjust=1, hjust=0, size=12) +
        draw_label(paste0("How to treat None regions: ", segs$action[1]),
                   x=.05, y=.71, vjust=1, hjust=0, size=12) +

        draw_label(paste0("Max. number of segments: ", segs$maxcp[1]),
                   x=.05, y=.66, vjust=1, hjust=0, size=12) +
        draw_label(paste0("Largest segment: ", segs$maxseg[1], " bins"),
                   x=.05, y=.63, vjust=1, hjust=0, size=12)

    plt_all <- plot_grid(plt_title, plt_content, ncol=1, rel_heights=c(0.07, 1))
    print(plt_all)
}


##############################
# Plot a seletion of #segments
quantiles = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
if (grepl('penalize', segs$action[1])) {
  ranges = min(c(max(segs$k)-1,segs$none_regions[1])) : max(segs$k)
} else {
  ranges = 1:max(segs$k)
}
for (num_seg in unique(quantile(ranges, quantiles, type = 3))) {

    message(" * Plotting segmation with ", num_seg, " segments...")
    local_bg_cols <- bg_cols[k == num_seg]
    plt <- ggplot(counts) +
        geom_rect(data = local_bg_cols, alpha = 0.5, size = 0.1,
                  aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = SV_class)) +
        scale_fill_manual(values = c(seg1 = "#edf8b1", seg2 = "#7fcdbb", seg3 = "#2c7fb8")) +
        geom_rect(data = none_regions, fill = "black", alpha = 0.33,
                  aes(xmin = start, xmax = end, ymin = -y_lim/2, ymax = y_lim/2)) +
        geom_rect(data = none_regions, fill = "black", alpha = 0.666,
                  aes(xmin = start, xmax = end, ymin = 0.99*y_lim, ymax = Inf)) +
        geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax = -w), fill='sandybrown') +
        geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax =  c), fill='paleturquoise4') +
        facet_wrap(~ paste(sample, cell), ncol = 1) +
        ylab("Watson | Crick") + xlab(NULL) +
        scale_x_continuous(breaks = pretty_breaks(12), labels = format_Mb) +
        scale_y_continuous(breaks = pretty_breaks(3)) +
        coord_cartesian(ylim = c(-y_lim, y_lim)) +
        theme_minimal() +
        theme(panel.spacing = unit(0, "lines"),
              axis.ticks.x = element_blank(),
              strip.background = element_rect(color = "#eeeeee", fill = "#eeeeee"),
              strip.text = element_text(size = 5),
              legend.position = "bottom") +
        ggtitle(paste0(num_seg, " segments on ", CHROM,  " (", basename(f_counts), ")"))

    print(plt)
}

dev.off()

