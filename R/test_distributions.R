# Test distributions
library(data.table)
library(ggplot2)
f = data.table(x = c(seq(100),
                     seq(100),
                     seq(100),
                     seq(100)), 
               y = c(dpois(seq(100),50), 
                     dnbinom(seq(100),50,0.5),
                     dnbinom(seq(100),50/3,0.25),
                     dnbinom(seq(100),50/9,0.1)), 
               type = c(rep("Poisson",100),
                        rep("NB(50,0.5)",100),
                        rep("NB(16.6,0.25)",100),
                        rep("NB(5.5,0.1)",100)))
ggplot(f) + aes(x,y,col=type) + geom_line() + theme_classic()

