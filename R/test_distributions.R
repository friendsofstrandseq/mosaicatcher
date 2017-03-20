# Test distributions
library(data.table)
library(ggplot2)
MAX=30
MEAN=10
plot_NB <- function(mean, max) {
    f = data.table(x = c(seq(max),
                         seq(max),
                         seq(max),
                         seq(max),
                         seq(max)),
                   y = c(dpois(seq(max),mean),
                         dnbinom(seq(max),mean,0.5),
                         dnbinom(seq(max),mean/4,0.2),
                         dnbinom(seq(max),mean/9,0.1),
                         dnbinom(seq(max),mean/49,0.02)),
                   type = c(rep("Poisson",max),
                            rep("NB(p=0.5)",max),
                            rep("NB(p=0.2)",max),
                            rep("NB(p=0.1)",max),
                            rep("NB(p=0.02)",max)))
    plt <- ggplot(f) +
        aes(x,y,col=type) +
        geom_line() +
        ggtitle(paste("Mean",mean)) +
        theme_classic()
    return(plt)
}

plot_NB(1200,2000)
plot_NB(10,20)
