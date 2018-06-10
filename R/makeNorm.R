library(dplyr)
library(data.table)
library(ggplot2)
library(assertthat)

args = commandArgs(trailingOnly = T)
if (length(args)<2) {
  message("Usage: Rscript create_normalization.R test-sample.txt.gz control-sample1.txt.gz [control-sample2.txt.gz [...]]")
  message("       Will write files <norm.txt> and <norm.pdf>")
  stop()
}


# Read test sample
message(" * Read test sample: ", args[1])
TEST = fread(paste("zcat", args[1]))

DF = NULL
for (x in args[2:length(args)]) {
  message(" * Read control sample: ", x)
  DF = rbind(DF,
             cbind(data.table(fread(paste("zcat",x)), 
                              file = sub(paste0("samples/(.*)_[0-9]+\\.txt.gz"), "\\1", x))))
}


# Normalize each cell by its mean
setkey(TEST, sample, cell, chrom, start, end)
setkey(DF, sample, cell, chrom, start, end)
TEST[, cov_norm := (c+w)/mean(c+w), by = .(sample, cell)]
DF[, cov_norm := (c+w)/mean(c+w), by = .(sample, cell)]


# Make sure bin positions are the same for all cells
message(" * Check whether the bins align in all samples")
bins = unique(DF[, .(chrom, start, end)])
DF[, assert_that(all(.SD == bins)), by = .(sample, cell), .SDcols = c("chrom","start","end")] %>% invisible
TEST[, assert_that(all(.SD == bins)), by = .(sample, cell), .SDcols = c("chrom","start","end")] %>% invisible


# Take mean coverage across all cells per sample
message(" * Get mean coverage across all control samples")
CONTROL = DF[, 
             .(mean = mean(cov_norm), sd = sd(cov_norm)), 
             by = .(chrom, start, end)]
setkey(CONTROL, chrom, start, end)
TEST = TEST[, 
            .(mean = mean(cov_norm), sd = sd(cov_norm)), 
            by = .(chrom, start, end)]


# Black-list weird bins:
CONTROL[mean < 1e-3, mean := 1e-3]
CONTROL[, class := ifelse(mean > 3 | mean < 1/3, "None","good")]
message(" * Black-list ", nrow(CONTROL[class == "None"]), "/", nrow(CONTROL), " control bins based on mean")


  # Plot information about black-listing
  p1 = ggplot(CONTROL) + 
    geom_histogram(aes(mean, fill = class), binwidth = 0.05) + 
    scale_x_log10(labels = scales::percent, breaks = c(0.1, 1, 10)) +
    geom_vline(xintercept = 3, linetype = "dotted") +
    geom_vline(xintercept = 1/3, linetype = "dotted") +
    xlab("Normalized coverage per bin") +
    ggtitle(paste("Mean normalized coverage in control samples"))


# Further black-list those bins that appear to have a high standard deviation (scaled to their mean)
  p2 = ggplot(CONTROL[mean>0]) + 
    geom_point(aes(mean, sd/mean, col = class, shape = sd/mean>1)) +
    scale_x_log10() + scale_y_log10() +
    geom_hline(yintercept = 1, linetype = "dotted") +
    ggtitle("Mean vs SD of norm coverage across bins")

none_before = nrow(CONTROL[class == "None"])
CONTROL[sd/mean > 1, class := "None"]
message(" * Black-list another ", nrow(CONTROL[class=="None"])-none_before, " bins based on the SD criterion")


# Look only at the correlation of TEST and CONTROL only in non-blacklisted bins:
X = merge(CONTROL[, .(chrom, start, end, control_cov = mean, class)], 
          TEST[,.(chrom, start, end, test_cov = mean)], 
          by = c("chrom", "start", "end"))
X[class!="None", control_cov := control_cov / mean(control_cov)]
X[class!="None", test_cov    := test_cov / mean(test_cov)]
slope = lm(test_cov ~ control_cov, data = X[class != "None"])$coefficients[2]
message(" * test and control have a linear relationship with slope ",
        round(slope,3), " and a correlation of Pearson's r = ", 
        round(cor(X[class!="None",control_cov], X[class!="None", test_cov]),3))


  # Now plot CONTROL  vs. TEST
  p3 <- ggplot(X[class != "None"]) + 
    aes(control_cov, test_cov, col = chrom) +
    geom_point(alpha = 0.5) +
    geom_smooth(aes(control_cov, test_cov), method = "lm", inherit.aes = F) +
    ggtitle(paste0("Correlation of mean coverages (Slope = ", slope, "; None bins removed)"))


# Apply scaling factor and write table  
X[, scalar := 1.0]
X[class!="None", scalar := 1/(((control_cov - mean(control_cov)) * slope) + mean(control_cov))]


message(" * Writing table to norm.txt")
write.table(X[, .(chrom, start, end, scalar, class)], "norm.txt", quote=F, col.names = T, sep="\t", row.names = F)


# Make plots
message(" * Writing plots to norm.pdf")
cairo_pdf("norm.pdf", width=8, height=6, onefile = T)
print(p1)
print(p2)
print(p3)
dev.off()



