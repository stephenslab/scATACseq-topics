# Using the Cicero co-accessibility data, generate a data structure
# containing a list of peaks for each gene.
library(tools)

# Import the Cicero co-accessibility data.
cicero <- readRDS(file.path("../data/Cusanovich_2018/processed_data",
                            "master_cicero_conns.rds"))
class(cicero) <- "data.frame"
cicero <- cicero[c("Peak1","Peak2","peak1.tss.gene_id",
                   "peak2.tss.gene_id","cluster")]
cicero <- subset(cicero,is.element(cluster,c(11,18,22,25)))
cicero <- transform(cicero,
                    Peak1             = as.character(Peak1),
                    Peak2             = as.character(Peak2),
                    peak1.tss.gene_id = factor(peak1.tss.gene_id),
                    peak2.tss.gene_id = factor(peak2.tss.gene_id))

# Here we extract, for each gene, the distal and proximal sites that
# are connected to Peak1. The two columns used in this calculation
# are: (1) "peak1.tss.gene_id", Ensemble gene id of the proximal TSS
# overlapping Peak1; (2) "Peak2", peak id of the second peak in the
# connection.
cicero_gene <- tapply(cicero$Peak2,cicero$peak1.tss.gene_id,unique,
                      simplify = FALSE)
# cicero_gene <- tapply(cicero$Peak1,cicero$peak1.tss.gene_id,unique,
#                       simplify = FALSE)

# Save these data to an .RData file.
save(list = "cicero_gene",file = "cicero_gene_kidney.RData")
resaveRdaFiles("cicero_gene_kidney.RData")
