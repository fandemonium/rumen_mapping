library(ggplot2)
library(DESeq2)
library(phyloseq)
library(reshape2)

# read in the parsed coverage file
coverage <- read.delim("all_coverages_parsed.txt", header=F)
# add in header
names(coverage) <- c("sample_name", "contig", "gene_start", "gene_end", "rmg_genome", "protein_id", "read_count", "mapped_len", "gene_length", "coverage")

# read in sample metadata
si <- read.delim("../shotgun_metadata.txt")
## ad hoc adding more information
si$sample_id <- data.frame(do.call('rbind', strsplit(as.character(si$Sample_Name), "_", fixed=T)))[, 2]
si$time <- data.frame(do.call('rbind', strsplit(as.character(si$Sample_Name), "_", fixed=T)))[, 1]

# save merged coverage for later use
coverage <- merge(coverage, si[, c(2, 3, 11, 14:16)], by.x="sample_name", by.y="Sample_Name")
coverage$group_id <- paste0(coverage$rmg_genome,  ":", coverage$protein_id, ":", coverage$gene_start, ":", coverage$gene_end)
saveRDS(coverage, "merged_coverage_w_sample_info.RDS")

# create the gene information table
cov.gene <- coverage[, c("group_id", "protein_id", "rmg_genome", "contig", "gene_start", "gene_end")]
# dereplicate
cov.gene <- cov.gene[! duplicated(cov.gene), ]
# make the protein_id the row.names
row.names(cov.gene) <- cov.gene$group_id

# create the simplified sample information
cov.si <- coverage[, c("sample_name", "PicoGreen_Conc", "Metadata1", "Clustering", "sample_id", "time")]
cov.si <- cov.si[! duplicated(cov.si), ]
# make the sample_name the row.names
row.names(cov.si) <- cov.si$sample_name

# run Rscript to loop through different cut off values:
# Rscript ~/Documents/repos/rumen_mapping/scripts/eval_diff_coverage_cutoff.R merged_coverage_w_sample_info.RDS

# 0.7 looks good
# parametric for fitType
 
cov <- subset(coverage, coverage >= 0.7)
# prepare for phyloseq object
# create the wide count table
cov.wide <- dcast(cov[, c("sample_name", "group_id", "read_count")], group_id ~ sample_name, value.var = "read_count")
cov.wide[is.na(cov.wide)] <- 0
# make group_id the row.names
row.names(cov.wide) <- cov.wide$group_id
# remove the group_id column
cov.wide$group_id <- NULL

# make the phyloseq object
cov.phy <- phyloseq(otu_table(as.matrix(cov.wide), taxa_are_rows=T), tax_table(as.matrix(cov.gene)))
sample_data(cov.phy) <- cov.si

# mds with treatment and clustering
## clustering doesn't separate anything
plot_ordination(cov.phy, ordinate(cov.phy, "MDS"), color = "sample_id", shape="Metadata1", label="Sample_Name") + geom_point(size = 3) + theme_bw()


# deseq2
cov.diagdds <- phyloseq_to_deseq2(cov.phy, ~ Metadata1)
# use parametric to estimate dispersion
cov.deseq <- DESeq(cov.diagdds, test="Wald", fitType="parametric")
# get result
cov.res <- data.frame(results(cov.deseq))
sig <- subset(cov.res, padj < 0.05)
sig$group_id <- row.names(sig)
####
## > sig
##               baseMean log2FoldChange    lfcSE     stat       pvalue
##k87_11225348_4 2.613768       20.06864 2.956836 6.787199 1.143318e-11
##k87_16477912_4 2.613768       20.06864 2.956836 6.787199 1.143318e-11
##                       padj     group_id
##k87_11225348_4 5.794905e-08 k87_11225348_4
##k87_16477912_4 5.794905e-08 k87_16477912_4
####
save.image("coverage_ana.RData")
#savehistory("coverage_ana_history.R")
