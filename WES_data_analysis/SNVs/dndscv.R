library("devtools")
# install_github("im3sanger/dndscv")
library("dndscv")
library("tidyr")


# By default, dNdScv assumes human data from the GRCh37/hg19 assembly.
# Precomputed reference file for GRCh38/hg38 at:
# https://github.com/im3sanger/dndscv_data/tree/master/data

# load("/home/rania/WM_WES/RefCDS_human_GRCh38_GencodeV18_recommended.rda")

RefCDS <- "/home/rania/WM_WES/RefCDS_human_GRCh38_GencodeV18_recommended.rda"


# Use as input a maf object and keep only the SNVs and INDELs
maf.data <- laml@data

patients.all <- maf.data[maf.data$Variant_Type != "CNV" ,]




MYD88 <- patients.all[patients.all$Hugo_Symbol == "MYD88", c("Hugo_Symbol","Start_Position","End_Position","Strand","Variant_Classification","HGVSp",
                                                             "Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode",
                                                             "t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count","n_alt_count")]


table(MYD88$Variant_Classification)



# dNdScv input prep: data.frame with five columns, snvs and indels
# (sampleID, chromosome, position, reference base and mutant base)

mutations <- patients.all[, c("Tumor_Sample_Barcode","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2")]

# Rename colnames
colnames(mutations) <- c("sampleID", "chr", "pos", "ref", "mut")

# Set all chromosomes to just the chromosome number
mutations$chr <- gsub("chr","", as.character(mutations$chr))
# head(mutations)
# nrow(mutations)

# remove some contigs that are not chrom 1- 22, X, Y
mutations <- mutations[mutations$chr %in% c(as.character(seq.int(1, 22)), "X", "Y") ,]
# unique(mutations$chr)

dndsout <- dndscv::dndscv(mutations, refdb = RefCDS, cv = NULL)


# dndscv outputs: Table of significant genes
# Result of neutrality tests at gene level
# P-values for substitutions are obtained by Likelihood-Ratio Tests.
# Q-values are obtained by Benjamini-Hodgberg’s multiple testing correction.
# Numbers of substitutions of each class observed in each gene & maximum-likelihood estimates (MLEs) of the dN/dS ratios for each gene,
# for missense (wmis), nonsense (wnon), essential splice site mutations (wspl) and indels (wind).
# The global q-value integrating all mutation types are available in the qglobal_cv and qallsubs_cv columns for analyses with and without indels, respectively.
sel_cv = dndsout$sel_cv
print(head(sel_cv, n=10), digits = 3)

# qglobal_cv when you also have indels
signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
rownames(signif_genes) = NULL
print(signif_genes)

# Search for genes of interest:
sel_cv[sel_cv$gene_name == "RSAD2",]

genes.with.qmis.lesser0.05 <- sel_cv[sel_cv$qmis_cv <= 0.05,]
print(genes.with.qmis.lesser0.05 )
genes.with.qtrunc.lesser0.05 <- sel_cv[sel_cv$qtrunc_cv <= 0.05,]
print(genes.with.qtrunc.lesser0.05)


# global MLEs for the dN/dS ratios across all genes
print(dndsout$globaldnds)

# annotated table of coding mutations
head(dndsout$annotmuts)

# θ<1 may reflect problems with the suitability of the dNdScv model for the dataset
print(dndsout$nbreg$theta)

# in very large datasets
# neutrality tests per gene: performed using a more traditional dN/dS model
# in which the local mutation rate for a gene is estimated exclusively from the synonymous mutations observed in the gene (dNdSloc)
signif_genes_localmodel = as.vector(dndsout$sel_loc$gene_name[dndsout$sel_loc$qall_loc<0.1])
print(signif_genes_localmodel)



# Using different substitution models
# In general, the full trinucleotide model is recommended for cancer genomic datasets as it typically provides the least biased dN/dS estimates

# 192 rates (used as default): sm = "192r_3w"
data("submod_192r_3w", package="dndscv")
colnames(substmodel) = c("syn","mis","non","spl")
head(substmodel)

# 12 rates (no context-dependence)
data("submod_12r_3w", package="dndscv")
colnames(substmodel) = c("syn","mis","non","spl")
head(substmodel)

# 2 rates (classic ts/tv model)
data("submod_2r_3w", package="dndscv")
colnames(substmodel) = c("syn","mis","non","spl")
head(substmodel)


dndsout.2 <- dndscv::dndscv(mutations, refdb = RefCDS, cv = NULL, sm = "12r_3w")

dndsout.3 <- dndscv::dndscv(mutations, refdb = RefCDS, cv = NULL, sm = "2r_3w")

# Lower AIC values indicate a better-fit model, dndsout is better
AIC(dndsout$poissmodel)
AIC(dndsout.2$poissmodel)
AIC(dndsout.3$poissmodel)
