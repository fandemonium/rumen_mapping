suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(phyloseq))
suppressMessages(library(reshape2))

args <- commandArgs(TRUE)

coverage <- readRDS(args[1])
#coverage <- readRDS( "merged_coverage_w_sample_info.RDS")

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


# subset for different cut offs
for (i in seq(0.4, 0.9, 0.1)){
	print(i)
	cov <- subset(coverage, coverage >= i)
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
	print("plot mds!!!")
	pdf(paste0("plots/cov_", i, "_mds.pdf"))
	p1 <- plot_ordination(cov.phy, ordinate(cov.phy, "MDS"), color = "Clustering", shape="Metadata1", label="sample_name") + geom_point(size = 3) + theme_bw()
	print(p1)
	dev.off()

	# deseq2
	cov.diagdds <- phyloseq_to_deseq2(cov.phy, ~ Metadata1)
	# use parametric to estimate dispersion
	for (j in c("local", "parametric")){
		print(j)
		cov.deseq <- DESeq(cov.diagdds, test="Wald", fitType=j)
		pdf(paste0("plots/cov_", i, "_dispersion_", j, ".pdf"))
		p2 <- plotDispEsts(cov.deseq, main=j)
		print(p2)
		dev.off()
	}
}

