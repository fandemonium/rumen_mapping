suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(phyloseq))
suppressMessages(library(reshape2))

args <- commandArgs(TRUE)

coverage <- readRDS(args[1])
#coverage <- readRDS( "merged_coverage_w_sample_info.RDS")
si <- read.delim(args[2])
#si <- read.delim("../shotgun_metadata.txt")

# create the gene information table
cov.gene <- coverage[, c("protein_id", "rmg_genome", "contig", "gene_start", "gene_end", "gc")]
# dereplicate
cov.gene <- cov.gene[! duplicated(cov.gene), ]
# make the protein_id the row.names
row.names(cov.gene) <- cov.gene$protein_id

# create the simplified sample information
cov.si <- si[, c("Sample_Name", "PicoGreen_Conc", "Metadata1")]
# make the sample_name the row.names
row.names(cov.si) <- cov.si$Sample_Name


# subset for different cut offs
for (i in seq(0.5, 0.9, 0.1)){
	print(i)
	cov <- subset(coverage, coverage >= i)
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
	print("plot mds!!!")
	pdf(paste0("plots/cov_", i, "_mds.pdf"))
	p1 <- plot_ordination(cov.phy, ordinate(cov.phy, "MDS"), color = "Metadata1") + geom_point()
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

