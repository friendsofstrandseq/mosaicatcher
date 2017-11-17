#
# Copyright (C) 2017 Sascha Meiers
# Distributed under the MIT software license, see the accompanying
# file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.

library(data.table)
library(ggplot2)
library(scales)
library(assertthat)

args = c("/Volumes/korbel/meiers/projects/20171024_strandseq_hackathon/unified_pipeline/counts/D2Rfb.100000_fixed.txt.gz",
         "/Volumes/korbel/meiers/projects/20171024_strandseq_hackathon/unified_pipeline/nonSVcellsGTprobs.reformatted2.txt",
         "/Volumes/korbel/meiers/projects/20171024_strandseq_hackathon/unified_pipeline/nonSVcellsGTprobs.plot")
args = commandArgs(trailingOnly = T)
if (length(args) != 3) {
    message("Usage: Rscript script.R count-file sv-prob-file output-prefix")
    stop()
}

#################
# Check arguments
zcat_command = "zcat"
f_counts = args[1]
f_prob   = args[2]
f_out    = args[3]
if (grepl('\\.pdf$', f_out)) f_out = substr(f_out, 1, nchar(f_out)-4)


#######################
# Other global settings
chrom_regex <- '^chr[1-9XY][0-9]?$'
format_Mb   <- function(x) {paste(comma(x/1e6), "Mb")}
cells_per_page <- 8


#############
# Read counts
if (grepl('\\.gz$',f_counts)) {
    counts = fread(paste(zcat_command, f_counts))
} else {
    counts = fread(f_counts)
}
assert_that(all(c("chrom","start","end","class","sample","cell","w","c") %in% colnames(counts)))
counts <- counts[ grepl(chrom_regex, chrom),]
setkey(counts, sample, cell, chrom)


##########################
# Read SV probability file
if (grepl('\\.gz$',f_prob)) {
    prob = fread(paste(zcat_command, f_prob))
} else {
    prob = fread(f_prob)
}
assert_that("chrom" %in% colnames(prob),
            "start" %in% colnames(prob),
            "end" %in% colnames(prob),
            "cell" %in% colnames(prob),
            "p_ref" %in% colnames(prob))

# Remove CN columns
if (length(grep("^p_cn", colnames(prob)))>0) {
    prob[, grep("^p_cn", colnames(prob)) := NULL]
}

# Add prior probabilities: this is not the correct way to do it!
prob$p_ref = prob$p_ref * 1.1

# Check that cells are the same as in "counts"
c1 = unique(counts[,.(sample,cell)])
c2 = unique(prob[,.(sample,cell)])
c3 = c1[c2,, on = .(sample,cell), nomatch=0]
if(nrow(c1) != nrow(c3)) {
    message("ERROR: Not all the cells from the count table are also covered in the probabilities")
    message("Cells in the count table:")
    print(paste(c1$sample,c1$cell))
    message("Cells in the probability table:")
    print(paste(c2$sample,c2$cell))
    message("Cells of the count table that are also covered by probs:")
    print(paste(c3$sample,c3$cell))
    stop()
}



#######################
# Reformat prob as long
prob_long <- melt(prob, 
                  id.vars = c("chrom","start","end","sample","cell"), 
                  measure.vars = colnames(prob)[grepl('^p_',colnames(prob))],
                  variable.name = "SV_class",
                  value.name    = "probability")
prob_long[, SV_class := substr(SV_class,3,nchar(as.character(SV_class)))]
# Select 
prob_long <- prob_long[, .SD[probability == max(probability),][order(probability),][1,], by = .(chrom,start,end,sample,cell)]
setkey(prob_long, sample, cell, chrom)





################
# Colors for SVs
manual_colors = c(ref = "#bbbbbb", 
                  del_hom = "dodgerblue4", del_h1 = "dodgerblue1", del_h2 = "dodgerblue3",
                  inv_hom = "seagreen",    inv_h1 = "seagreen1",   inv_h2 = "seagreen3",
                  dup_hom = "red4",        dup_h1 = "red1",        dup_h2 = "red3",
                                           idup_h1 = "yellow1",    idup_h2 = "yellow3")
levels(prob_long$SV_class)

##########################
# Plot, one file per chrom
y_lim = 3 * counts[,median(w+c)]
n_cells = nrow(unique(counts[,.(sample,cell)]))

for (CHROM in unique(counts[, chrom])) {
    
    message("Plotting ", CHROM)
    cairo_pdf(paste0(f_out,".",CHROM,".pdf"), width=14, height=10, onefile = T)

    i = 1
    while (i < n_cells) {
        
        CELLS = unique(counts[,.(sample,cell)])[i:(min(i+cells_per_page-1,n_cells))]    
        setkey(CELLS, sample, cell)
        
        local_counts = counts[CELLS, on = .(sample,cell), nomatch = 0][chrom == CHROM]
        local_probs = prob_long[CELLS, on = .(sample,cell), nomatch=0][chrom == CHROM,]
        
    
        plt <- ggplot(local_counts)
        
        # Add GT colors:
        if(nrow(local_probs)>0) {
            plt <- plt + 
                geom_rect(data = local_probs, alpha = 0.2, size = 0.1,
                          aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = SV_class)) +
                scale_fill_manual(values = manual_colors)
        }
        
        plt <- plt +
            geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax = -w), fill='sandybrown') +
            geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax =  c), fill='paleturquoise4') +
            facet_grid(sample + cell ~ .) +
            ylab("Watson | Crick") + xlab(NULL) +
            scale_x_continuous(breaks = pretty_breaks(12), labels = format_Mb) +
            scale_y_continuous(breaks = pretty_breaks(3)) + 
            coord_cartesian(ylim = c(-y_lim, y_lim)) +
            theme_minimal() +
            theme(panel.spacing = unit(0, "lines"),
                  axis.ticks.x = element_blank(),
                  #strip.background = element_rect(fill = NA, colour=NA), 
                  legend.position = "bottom") +
            ggtitle(CHROM)
        
        print(plt)
        i = i + cells_per_page
    }
    dev.off()
}
 


#plt + geom_rect(data = prob_long, alpha = 0.2, size = 0.1,
#              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = SV_class)) +
#    geom_label(data = PLOT_M, inherit.aes = F, x = 0, y = Inf, vjust = 1,
#               aes(label = paste0(max_class, ", ", round(max_class_perc*100), "%")))

    
