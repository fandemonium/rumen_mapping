library(ggplot2)
library(DESeq2)
library(phyloseq)
library(reshape2)

# read in the parsed coverage file
coverage <- read.delim("all_coverage_parsed.txt", header=F)
# add in header
names(coverage) <- c("sample_name", "contig", "gene_start", "gene_end", "rmg_genome", "protein_id", "gc", "strand", "read_count", "mapped_len", "gene_length", "coverage")

# read in sample metadata
si <- read.delim("../shotgun_metadata.txt")

# save merged coverage for later use
coverage <- merge(coverage, si[, c(3, 11, 14)], by.x="sample_name", by.y="Sample_Name")
saveRDS(coverage, "merged_coverage_w_sample_info.RDS")

# create the gene information table
cov.gene <- coverage[, c("protein_id", "rmg_genome", "contig", "gene_start", "gene_end", "gc")]
# dereplicate
cov.gene <- cov.gene[! duplicated(cov.gene), ]
# make the protein_id the row.names
row.names(cov.gene) <- cov.gene$protein_id

# create the simplified sample information
cov.si <- si[, c("Sample_Name", "PicoGreen_Conc", "Metadata1", "Clustering")]
# make the sample_name the row.names
row.names(cov.si) <- cov.si$Sample_Name
## ad hoc adding more information
cov.si$sample_id <- data.frame(do.call('rbind', strsplit(as.character(cov.si$Sample_Name), "_", fixed=T)))[, 2]
cov.si$time <- data.frame(do.call('rbind', strsplit(as.character(cov.si$Sample_Name), "_", fixed=T)))[, 1]

# run Rscript to loop through different cut off values:
# Rscript eval_diff_coverage_cutoff.R merged_coverage_w_sample_info.RDS ../shotgun_metadata.txt

# 0.7 looks good
 
cov <- subset(coverage, coverage >= 0.7)
# prepare for phyloseq object
# create the wide count table
cov.wide <- dcast(cov[, c("sample_name", "protein_id", "read_count")], protein_id ~ sample_name, value.var = "read_count")
cov.wide[is.na(cov.wide)] <- 0
# make protein_id the row.names
row.names(cov.wide) <- cov.wide$protein_id
# remove the protein_id column
cov.wide$protein_id <- NULL

# make the phyloseq object
cov.phy <- phyloseq(otu_table(as.matrix(cov.wide), taxa_are_rows=T), tax_table(as.matrix(cov.gene)))
sample_data(cov.phy) <- cov.si

# mds with treatment and clustering
## clustering doesn't separate anything
plot_ordination(cov.phy, ordinate(cov.phy, "MDS"), color = "sample_id", shape="Metadata1", label="Sample_Name") + geom_point(size = 3) + theme_bw()


# deseq2
cov.diagdds <- phyloseq_to_deseq2(cov.phy, ~ Metadata1)
# use local to estimate dispersion
cov.deseq <- DESeq(cov.diagdds, test="Wald", fitType="local")
# get result
cov.res <- data.frame(results(cov.deseq))
sig <- subset(cov.res, padj < 0.05)
sig$protein_id <- row.names(sig)
####
## > sig
##               baseMean log2FoldChange    lfcSE     stat       pvalue
##k87_11225348_4 2.613768       20.06864 2.956836 6.787199 1.143318e-11
##k87_16477912_4 2.613768       20.06864 2.956836 6.787199 1.143318e-11
##                       padj     protein_id
##k87_11225348_4 5.794905e-08 k87_11225348_4
##k87_16477912_4 5.794905e-08 k87_16477912_4
####
save.image("coverage_ana.RData")
#savehistory("coverage_ana_history.R")
