# TO DO: Explain here what this script does, and how to use it.
library(data.table)
library(Matrix)

# Load the fragment counts as a 2,953 x 491,437 matrix.
x <- readLines("GSE96769_scATACseq_counts.txt.gz",n = 1)
x <- unlist(strsplit(x,"\t",fixed = TRUE))
x <- unlist(strsplit(x,";",fixed = TRUE))
x <- x[-1]
dat <- fread("GSE96769_scATACseq_counts.txt.gz",sep = "\t",skip = 1)
class(dat) <- "data.frame"
names(dat) <- c("i","j","x")
counts <- sparseMatrix(i = dat$i,j = dat$j,x = dat$x)
counts <- t(counts)
rownames(counts) <- x
rm(x,dat)

# TO DO: Explain here what this code chunk is doing.
y <- table(summary(counts)$x)
x <- names(y)
y <- as.numeric(y)
plot(x,y,pch = 20,log = "y",
     xlab = "number of fragments mapping to peak",
     ylab = "number of peaks")
lines(x,y)

# TO DO: save counts.
