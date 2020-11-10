# Short script to load the Buenrostro et al (2018) data ("Integrated
# single-cell analysis maps the continuous regulatory landscape of
# human hematopoietic differentiation") downloaded from Gene Omnibus,
# accession GSE96769.
library(tools)
library(data.table)
library(Matrix)

# Load the fragment counts as a 2,953 x 491,437 sparse matrix.
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

# PLot the distribution of fragment counts.
y <- table(summary(counts)$x)
x <- names(y)
y <- as.numeric(y)
plot(x,y,pch = 20,log = "y",cex = 0.65,
     xlab = "number of fragments mapping to peak",
     ylab = "number of peaks")
lines(x,y)

# Save the counts matrix.
save(list = "counts",file = "buenrostro2018.RData")
resaveRdaFiles("buenrostro2018.RData")
