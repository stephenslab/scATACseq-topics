# Downloaded seq_gene.md.gz from
# https://ftp.ncbi.nih.gov/genomes/MapView/Mus_musculus/sequence/BUILD.37.2/initial_release/
# gunzip seq_gene.md.gz
# remove "#" from first line
# gzip seq_gene.md
library(tools)
library(data.table)
seq_gene <- fread("seq_gene.md.gz",sep = "\t",quote = "",
                  header = TRUE,stringsAsFactors = FALSE)
class(seq_gene) <- "data.frame"
seq_gene <- subset(seq_gene,
                   feature_type == "GENE" &
                   grepl("MGSCv37",group_label,fixed = TRUE))
seq_gene <- seq_gene[c("chromosome","chr_start","chr_stop",
                       "feature_name","feature_id","group_label")]
seq_gene <-
  transform(seq_gene,
            chromosome  = factor(chromosome),
            group_label = factor(group_label),
            feature_id  = as.numeric(sapply(strsplit(feature_id,":"),"[[",2)))
save(list = "seq_gene",file = "mm9_seq_gene.RData")
resaveRdaFiles("mm9_seq_gene.RData")
