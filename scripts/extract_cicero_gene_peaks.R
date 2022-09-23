# Using the Cicero co-accessibility data, generate a data structure
# containing a list of peaks for each gene.
library(Matrix)
library(tools)

# Import the Cicero co-accessibility data.
cicero <- readRDS(file.path("../data/Cusanovich_2018/processed_data",
                            "master_cicero_conns.rds"))
class(cicero) <- "data.frame"
cicero <- cicero[c("Peak1","Peak2","peak1.tss.gene_id",
                   "peak2.tss.gene_id","cluster")]
cicero <- subset(cicero,is.element(cluster,c(11,18,22,25)))
cicero <- transform(cicero,
                    Peak1             = factor(as.character(Peak1)),
                    Peak2             = factor(as.character(Peak2)),
                    peak1.tss.gene_id = factor(peak1.tss.gene_id),
                    peak2.tss.gene_id = factor(peak2.tss.gene_id))

# Here we extract, for each gene, the distal and proximal sites that
# are connected to Peak1. The two columns used in this calculation
# are: (1) "peak1.tss.gene_id", Ensemble gene id of the proximal TSS
# overlapping Peak1; (2) "Peak2", peak id of the second peak in the
# connection.
genes <- levels(cicero$peak1.tss.gene_id)
peaks <- levels(cicero$Peak2)
cicero <- transform(cicero,
                    Peak1 = as.numeric(Peak1),
                    Peak2 = as.numeric(Peak2))
cicero_gene <- tapply(cicero$Peak2,cicero$peak1.tss.gene_id,unique,
                      simplify = FALSE)
# cicero_gene <- tapply(cicero$Peak1,cicero$peak1.tss.gene_id,unique,
#                       simplify = FALSE)

# Encode these data as a sparse matrix cicero_gene_mat such that
# cicero_gene_mat[i,j] = 1 if and only if peak i is connected to gene j.
n <- length(peaks)
m <- length(genes)
j <- cicero_gene
for (t in 1:m)
  j[[t]][] <- t
i <- unlist(cicero_gene)
j <- unlist(j)
cicero_gene_mat <- sparseMatrix(i = i,j = j,dims = c(n,m))
cicero_gene_mat <- as(cicero_gene_mat,"dgCMatrix")
rownames(cicero_gene_mat) <- peaks
colnames(cicero_gene_mat) <- genes

# Save these data to an .RData file.
save(list = "cicero_gene_mat",file = "cicero_gene_kidney.RData")
resaveRdaFiles("cicero_gene_kidney.RData")
