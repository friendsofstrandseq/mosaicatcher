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
    plt <- ggplot(perBin) + aes(median) + geom_histogram(binwidth=1e-1) +
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


