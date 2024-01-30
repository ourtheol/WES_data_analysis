# if (!require("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install("maftools")

library("maftools")
library("tidyverse")
library("ggplot2")
library("RColorBrewer")

library("vcfR")
library("GenomicRanges")
library("CNVRanger")
library("AnnotationHub")

library("mclust")



setwd("/home/rania/Documents/Tina/WES/maftools_plots/")


# Path to the somatic SNVs maf files
somatic_snvs_maf_path <- "/home/rania/Documents/Tina/WES/mafs/new_mafs/snvs/SNVs/"

# Path to the somatic INDELSs maf files
somatic_indels_maf_path <- "/home/rania/Documents/Tina/WES/mafs/new_mafs/indels/INDELS/"

# path to the vcfs output from cnv_facets
cnv.vcf.files.dir <- "/home/rania/Documents/Tina/WES/single_samples/cnv_vcfs/"




# function: vcf2df -------------------------------------------------------------
# insert the path of a vcf file output from CNV_FACETS, and
# return a data frame with colnames: chr start end patient_id tcn svtype

vcf2df <- function(vcf.path) {
  
  vcf <- read.vcfR(vcf.path)
  
  chrom <- substring(as.data.frame(getFIX(vcf))$CHROM, 4)   # remove the "chrom"
  
  pos <- as.data.frame(getFIX(vcf))$POS
  
  # access the INFO column of the vcf file
  end <- INFO2df(vcf)$END
  
  svtype <- INFO2df(vcf)$SVTYPE
  
  tcn <- INFO2df(vcf)$TCN_EM
  
  num_mark <- INFO2df(vcf)$NUM_MARK
  
  cnlr_median <- INFO2df(vcf)$CNLR_MEDIAN
  
  # data frame with the CNVs called by cnv_facets
  cnv.df <- data.frame(chrom, pos, end, svtype, tcn, num_mark, cnlr_median)
  
  return(cnv.df)
  
}



# Prep the input to maftools ---- 

## Variant_Classifications
nonSilentPLUSsilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                        "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",
                        
                        "3'UTR", "5'UTR", "3'Flank", "Targeted_Region", "Silent", "Intron", "RNA", "IGR",
                        "Splice_Region", "5'Flank", "lincRNA", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame",
                        "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del")

## Set colors for each Variant_Classification
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
colourCount = length(unique(nonSilentPLUSsilent))
nonSilentPLUSsilent.cols = getPalette(colourCount)
names(nonSilentPLUSsilent.cols) = nonSilentPLUSsilent

# ## Set colors for MYD88L265P
# mutVSwt.cols = c("mediumaquamarine", "mediumpurple2")
# names(mutVSwt.cols) = c("mut","wt")
# 
# ## Set colors for stage of WM
# stage.cols = c("khaki3", "salmon1", "indianred2", "orangered3", "chartreuse4")
# names(stage.cols) = c("IGM-MGUS", "aWM_stable", "aWM_progressed", "sWM", "treated")


## clinical info df for each patient
clinical.info <- read.table("/home/rania/Documents/Tina/WES/single_samples/annotations.txt")
clinical.info$Tumor_Sample_Barcode <- rownames(clinical.info)

# Reorder the patient_IDs according to the input matrix (wide) in the oncoPrint()
clinical.info$State <- factor(clinical.info$State, levels = c("IGM-MGUS", "aWM_stable", "aWM_progressed", "sWM", "treated"))
clinical.info <- clinical.info[order(clinical.info$State, clinical.info$cAst_PCR),]
patient.order <- as.vector(rownames(clinical.info))

print(head(as.data.frame(clinical.info)))


# Patient states 
IGM_MGUS <- clinical.info[clinical.info$State == "IGM-MGUS",]$Tumor_Sample_Barcode
aWM_progressed <- clinical.info[clinical.info$State == "aWM_progressed",]$Tumor_Sample_Barcode
aWM_stable <- clinical.info[clinical.info$State == "aWM_stable",]$Tumor_Sample_Barcode
sWM <- clinical.info[clinical.info$State == "sWM",]$Tumor_Sample_Barcode





## vector of maf files (snv and indel results)
mafs.read <- c()



###### SNVs
somatic_snvs <- list.files(path = somatic_snvs_maf_path)

for (maf in somatic_snvs){
  maf.read = read.maf(maf = paste0(somatic_snvs_maf_path, maf),
                      vc_nonSyn = nonSilentPLUSsilent,  
                      clinicalData = clinical.info)
  mafs.read=append(mafs.read, maf.read)
}



###### INDELs
somatic_indels <- list.files(path = somatic_indels_maf_path)

for (maf in somatic_indels){
  maf.read = read.maf(maf = paste0(somatic_indels_maf_path, maf),
                      vc_nonSyn = nonSilentPLUSsilent) 
  mafs.read=append(mafs.read, maf.read)
}



###### CNVs

# Read all the vcf files (CNV_FACETS output) in a dir 
# and create a data-frame with all CNVs in the cohort

# empty data-frame to host all the cohort cnvs
patients.cnvs.df <- data.frame()

for (vcf.file in list.files(cnv.vcf.files.dir)) {
  
  # get the patient_ID
  patient_ID <- strsplit(vcf.file, split = "\\.")[[1]][1]
  
  # get a df for each patient's cnvs
  cnv.df <- vcf2df(paste0(cnv.vcf.files.dir, vcf.file))
  
  # Remove Neutral CNVs    ########################### read about that
  cnv.df <- cnv.df[cnv.df$svtype != "NEUTR",]
  
  # insert the patient_ID into the patient's df
  cnv.df$patient_ID <- patient_ID
  
  patients.cnvs.df <- rbind(patients.cnvs.df, cnv.df)
  
}




# Rename the column names to the exprected from "makeGRangesListFromDataFrame" function
colnames(patients.cnvs.df) <- c("chr", "start", "end", "svtype", "state", "num_mark", "cnlr_median",  "patient_ID")

# Read the R data-frame and group the calls by patient ID and convert them to a GRangesList
grl <- GenomicRanges::makeGRangesListFromDataFrame(patients.cnvs.df[, c("chr", "start", "end", "svtype", "state", "patient_ID")], 
                                                   split.field="patient_ID", keep.extra.columns=TRUE)




# Map the genomic locations of the CNVs to the genes located there 

# Overlap analysis of CNVs with functional genomic regions
ah <- AnnotationHub::AnnotationHub()

ahDb <- AnnotationHub::query(ah, pattern = c("Homo sapiens", "EnsDb", "GRCh38"))

# retrieve selected record
ahEdb <- ahDb[["AH53211"]]



genes <- ensembldb::genes(ahEdb)

genes <- subset(genes, gene_biotype %in% c("protein_coding",
                                           "pseudogene", "IG_V_pseudogene",
                                           "TR_V_gene","TR_V_pseudogene","TR_C_gene","TR_D_gene","TR_J_gene",
                                           "TR_J_pseudogene","IG_C_gene", "IG_C_pseudogene","IG_J_gene","IG_J_pseudogene","IG_D_gene",
                                           "IG_V_gene","IG_pseudogene")
)



# dataframe input into maftools with colnames: GeneName, sampleName, CNVType ("Amp" or "Del") 
cnTable <- data.frame()

for (i in 1:length(grl)){

  # insert the patient_id as a meta data column
  grl[[i]]$patient_id <- names(grl)[i]
  
  # # Match genomic locations with affected genes
  olaps <- GenomicRanges::findOverlaps(genes, grl[[i]], ignore.strand=TRUE)
  qh <- S4Vectors::queryHits(olaps)
  sh <- S4Vectors::subjectHits(olaps)
  cgenes <- genes[qh]
  # annotate the CNV type (gain/loss) for genes overlapping with CNVs
  cgenes$svtype <- grl[[i]]$svtype[sh]
  cgenes$patient_id <- grl[[i]]$patient_id[sh]
  metadata.df <- as.data.frame(cgenes@elementMetadata[, c("symbol","patient_id","svtype")])
  
  cnTable <- rbind(cnTable, metadata.df)
}


# Set colnames as expected by read.maf()
colnames(cnTable) <- c("GeneName", "sampleName", "CNVType")

# Set the CNVType as "Amp" or "Del"
unique(cnTable$CNVType)
# rename
cnTable[cnTable$CNVType == "DEL", ]$CNVType <- "Del"
cnTable[cnTable$CNVType == "HEMIZYG", ]$CNVType <- "Del"
cnTable[cnTable$CNVType == "DUP", ]$CNVType <- "Amp"
cnTable[cnTable$CNVType == "LOH", ]$CNVType <- "Amp"  ## ??
cnTable[cnTable$CNVType == "DUP-LOH", ]$CNVType <- "Amp"


# Data into maftools ----

laml <- merge_mafs(c(mafs.read),
                   vc_nonSyn = c(nonSilentPLUSsilent, "Amp", "Del"),
                   cnTable = cnTable)   ## here we can add info for genes affected by CNVs using the param: cnTable
# Typing laml shows basic summary of MAF file.
print(laml)


#Shows sample summary.
getSampleSummary(laml)
#Shows gene summary.
getGeneSummary(laml)
#shows clinical data associated with samples
getClinicalData(laml)
# #Shows all fields in MAF
getFields(laml)
# #Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'laml')



# oncoplots ----

pdf(file = "maftools_out.pdf")

# plot the summary of the maf file: number of variants in each sample as a stacked barplot and variant types as a boxplot summarized by Variant_Classification.
plotmafSummary(maf = laml, addStat = 'median', titvRaw = FALSE, top = 25, fs = 0.8, color = nonSilentPLUSsilent.cols, showBarcodes = TRUE, textSize = 0.3)

# oncoplot for top mutated genes
oncoplot(maf = laml,
         top = 50,
         fontSize = 0.3, annotationFontSize = 1, legendFontSize = 1, titleFontSize = 0.8,
         showTumorSampleBarcodes = TRUE,
         SampleNamefontSize = 0.7,
         genesToIgnore = "Unknown",
         clinicalFeatures = c("State","cAst_PCR"),
         sortByAnnotation = TRUE,   ##sorts only according to the 1st clinical feature
         groupAnnotationBySize = FALSE,
         annotationOrder = c("IGM-MGUS", "aWM_stable", "aWM_progressed", "sWM", "treated"),
         includeColBarCN = FALSE
)

# oncoplot with genes Tina had asked me about
oncoplot(maf = laml,
         genes = c("CXCR4","TBL1XR1","PTPN13","MALT1","BCL10","NFKB1","NKFB2","NFKBIB","NFKBIZ","UDRL1F","TP53","ATM","TRRAP"),
         fontSize = 0.3, annotationFontSize = 1, legendFontSize = 1, titleFontSize = 0.8,
         showTumorSampleBarcodes = TRUE,
         clinicalFeatures = c("State","cAst_PCR"),
         sortByAnnotation = TRUE,   ##sorts only according to the 1st clinical feature
         groupAnnotationBySize = FALSE,
         annotationOrder = c("IGM-MGUS", "aWM_stable", "aWM_progressed", "sWM", "treated"),
         includeColBarCN = FALSE
)


# oncoplot with genes dndscv genes
oncoplot(maf = laml,
         genes = c("MYD88", "CXCR4", "FOXD4L1", "GYPB", "RSAD2", "CD79B", "ZYG11B"),
         fontSize = 0.3, annotationFontSize = 1, legendFontSize = 1, titleFontSize = 0.8,
         showTumorSampleBarcodes = TRUE,
         clinicalFeatures = c("State","cAst_PCR"),
         sortByAnnotation = TRUE,   ##sorts only according to the 1st clinical feature
         groupAnnotationBySize = FALSE,
         annotationOrder = c("IGM-MGUS", "aWM_stable", "aWM_progressed", "sWM", "treated"),
         includeColBarCN = FALSE
)


# oncoplot with immunoglobuline genes (only the ones with mutations)
immunoglobuline.genes = read.table(file="/home/rania/Documents/Tina/WES/maftools_plots/immunoglobulin_genes.txt", header = FALSE)
immunoglobuline.genes.vec <- immunoglobuline.genes[ , 1]
oncoplot(maf = laml,
         top = 50,
         # minMut = 8,
         genes = immunoglobuline.genes.vec,
         altered = TRUE,
         SampleNamefontSize = 0.7,
         showTumorSampleBarcodes = TRUE,
         fontSize = 0.05, annotationFontSize = 1, legendFontSize = 1, titleFontSize = 0.8,
         clinicalFeatures = c("State","cAst_PCR"),
         sortByAnnotation = TRUE,   ##sorts only according to the 1st clinical feature
         groupAnnotationBySize = FALSE,  ##orders annotations by their size
         annotationOrder = c("IGM-MGUS", "aWM_stable", "aWM_progressed", "sWM", "treated"),
         includeColBarCN = FALSE,
         fill = FALSE  # to remove the non-mutated genes
)

# oncoplot for top mutated genes, remove the immunoglobulin genes!!
oncoplot(maf = laml,
         top = 85,
         fontSize = 0.3, annotationFontSize = 1, legendFontSize = 1, titleFontSize = 0.8,
         showTumorSampleBarcodes = TRUE,
         SampleNamefontSize = 0.7,
         genesToIgnore = c("Unknown", immunoglobuline.genes.vec, "IGHM", "IGHJ6"),
         clinicalFeatures = c("State","cAst_PCR"),
         sortByAnnotation = TRUE,   ##sorts only according to the 1st clinical feature
         groupAnnotationBySize = FALSE,
         annotationOrder = c("IGM-MGUS", "aWM_stable", "aWM_progressed", "sWM", "treated"),
         includeColBarCN = FALSE
)

# oncoplot with pathways -----


par(oma=c(2,2,2,2)) 
oncoplot(maf = laml,

         drawColBar = FALSE,
         # gene_mar = 14,
         pathways = "sigpw",
         selectedPathways = c("NFKB_signaling", "MAPK_signaling", "Genome_integrity", "Epigenetics_DNA_modifiers", "Immune_signaling"),
         showTumorSampleBarcodes = TRUE,
         fontSize = 0.5, annotationFontSize = 0.8, legendFontSize = 1, titleFontSize = 0.8,
         sortByAnnotation = TRUE,
         annotationOrder = c("IGM-MGUS", "aWM_stable", "aWM_progressed", "sWM", "treated"),
         clinicalFeatures = c("State","cAst_PCR")
)




# Get ALL the genes in selected pathways

library("msigdbr")

# show available species
msigdbr_species()

all_gene_sets = msigdbr(species = "Homo sapiens")
head(all_gene_sets)

# show all the available collections
msigdbr_collections()

all_gene_sets %>%
  dplyr::filter(gs_cat == "C2", gs_subcat == "CP:BIOCARTA") %>%
  head()

# C2: curated gene sets
biocarta_gene_sets <- dplyr::filter(all_gene_sets, gs_cat == "C2", gs_subcat == "CP:BIOCARTA")


# NFKB, ERK/MAPK, DNA  repair pathway, epigenetic signaling
pathways <- biocarta_gene_sets[biocarta_gene_sets$gs_name %in%  c("BIOCARTA_NFKB_PATHWAY", "BIOCARTA_ERK_PATHWAY", "BIOCARTA_MAPK_PATHWAY", "BIOCARTA_RB_PATHWAY", "BIOCARTA_P53_PATHWAY"), c("human_gene_symbol", "gs_name")]
colnames(pathways) <- c("Genes", "Pathway")


par(oma=c(2,2,2,2))
oncoplot(maf = laml,
         drawColBar = FALSE,
         pathways = pathways, gene_mar = 8,
         showTumorSampleBarcodes = TRUE,
         clinicalFeatures = c("State","cAst_PCR"),
         fontSize = 0.2, annotationFontSize = 0.5, legendFontSize = 0.4, titleFontSize = 0.7,
         sortByAnnotation = TRUE,
         annotationOrder = c("IGM-MGUS", "aWM_stable", "aWM_progressed", "sWM", "treated")
         )



# descriptive plots ----

# SNPs classified into Transitions and Transversions, summarization tables
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv, showBarcodes = TRUE, textSize = 0.5)

# lollipop plot for selected gene
## SOS: protein structure must be available!!
lollipopPlot(maf = laml,gene = 'MYD88',AACol ='HGVSp_Short',showMutationRate = TRUE,labelPos = "all",printCount = TRUE)
lollipopPlot(maf = laml,gene = 'IGLL5',AACol ='HGVSp_Short',showMutationRate = TRUE,labelPos = "all",printCount = TRUE)
lollipopPlot(maf = laml,gene = 'MUC2',AACol ='HGVSp_Short',showMutationRate = TRUE,labelPos = "all",printCount = TRUE)
lollipopPlot(maf = laml,gene = 'DDX11',AACol ='HGVSp_Short',showMutationRate = TRUE,labelPos = "all",printCount = TRUE)
lollipopPlot(maf = laml,gene = 'KMT2C',AACol ='HGVSp_Short',showMutationRate = TRUE,labelPos = "all",printCount = TRUE)

# # draw protein domains
# plotProtein(gene = "TP53", refSeqID = "NM_000546")

# genomic loci with localized hyper-mutations
# rainfall plots: plotting inter variant distance on a linear genomic scale
# also highlight regions where potential changes in inter-event distances are located
rainfallPlot(maf = laml, tsb = "T1_608", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")
rainfallPlot(maf = laml, tsb = "T2_608", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")

rainfallPlot(maf = laml, tsb = "T_650_1", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")
rainfallPlot(maf = laml, tsb = "T_650_2", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")

rainfallPlot(maf = laml, tsb = "T_655_1", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")
rainfallPlot(maf = laml, tsb = "T_655_2", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")

rainfallPlot(maf = laml, tsb = "T_790_1", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")
rainfallPlot(maf = laml, tsb = "T_790_2", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")

rainfallPlot(maf = laml, tsb = "T1_805", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")
rainfallPlot(maf = laml, tsb = "T2_805", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")

rainfallPlot(maf = laml, tsb = "T_650_1", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")
rainfallPlot(maf = laml, tsb = "T_650_2", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")

rainfallPlot(maf = laml, tsb = "T1_2308", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")
rainfallPlot(maf = laml, tsb = "T_2308_2", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")

rainfallPlot(maf = laml, tsb = "T_1928", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")
rainfallPlot(maf = laml, tsb = "T2_1928", detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38")


rainfallPlot(maf = laml, tsb = IGM_MGUS, detectChangePoints = FALSE, pointSize = 0.4, ref.build = "hg38", fontSize = 0)
rainfallPlot(maf = laml, tsb = aWM_stable, detectChangePoints = FALSE, pointSize = 0.4, ref.build = "hg38", fontSize = 0)
rainfallPlot(maf = laml, tsb = aWM_progressed, detectChangePoints = FALSE, pointSize = 0.4, ref.build = "hg38", fontSize = 0)
rainfallPlot(maf = laml, tsb = sWM, detectChangePoints = TRUE, pointSize = 0.4, ref.build = "hg38", fontSize = 0)



## compares mutation load in input MAF against all of 33 TCGA cohorts derived from MC3 project
par(oma=c(4,4,4,4))  # outer margins
tcgaCompare(maf = laml, cohortName = 'WM', logscale = TRUE, cohortFontSize = 0.5, axisFontSize = 0.5)
# primary site of cancer as labels
par(oma=c(4,4,4,4))
tcgaCompare(maf = laml, cohortName = 'WM', primarySite = TRUE, logscale = TRUE, cohortFontSize = 0.5, axisFontSize = 0.5)

## plotting VAF
# boxplot which quickly helps to estimate clonal status of top mutated genes
# (clonal genes usually have mean allele frequency around ~50% assuming pure sample)
# no t_vaf column, plotVaf will crate one since we have t_ref_count and t_alt_count
## t_vaf field is missing, but found t_ref_count & t_alt_count columns. Estimating vaf..
par(oma=c(2,2,2,2))
plotVaf(maf = laml, vafCol = NULL, top = 25)
# plotVaf(maf = laml, vafCol = NULL, genes = c("MYD88"))

# # checks for drug–gene interactions and gene druggability information compiled from Drug Gene Interaction DB
par(oma=c(2,2,2,2))
drugInteractions(maf = laml, fontSize = 0.75, top = 40)
 
# detecting cancer driver genes based on positional clustering
# identify cancer genes (driver) from a given MAF
# Most of the variants in cancer causing genes are enriched at few specific loci (hot-spots).
# This method takes advantage of such positions to identify cancer genes.
laml.sig = oncodrive(maf = laml, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
head(laml.sig)
par(oma=c(2,2,2,2))
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
# what domain in given cancer cohort, is most frequently affected
# Summarizes amino acid positions and annotates them with pfam domain information
par(oma=c(2,2,2,2))
laml.pfam = pfamDomains(maf = laml, AACol = 'HGVSp_Short', top = 10)

#exclusive/co-occurance event analysis on top mutated genes.
somaticInteractions(maf = laml, top = 30, pvalue = c(0.05, 0.1))


# checks for enrichment of known Oncogenic Signaling Pathways in TCGA cohorts
OncogenicPathways(maf = laml)

# visualize complete pathway
PlotOncogenicPathways(maf = laml, pathways = "NOTCH", showTumorSampleBarcodes = TRUE)
PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS", showTumorSampleBarcodes = TRUE)
PlotOncogenicPathways(maf = laml, pathways = "WNT", showTumorSampleBarcodes = TRUE)
PlotOncogenicPathways(maf = laml, pathways = "Hippo", showTumorSampleBarcodes = TRUE)
PlotOncogenicPathways(maf = laml, pathways = "PI3K", showTumorSampleBarcodes = TRUE)
PlotOncogenicPathways(maf = laml, pathways = "Cell_Cycle", showTumorSampleBarcodes = TRUE)
PlotOncogenicPathways(maf = laml, pathways = "MYC", showTumorSampleBarcodes = TRUE)
PlotOncogenicPathways(maf = laml, pathways = "NRF2", showTumorSampleBarcodes = TRUE)
PlotOncogenicPathways(maf = laml, pathways = "TP53", showTumorSampleBarcodes = TRUE)


dev.off()



# clinical info used in analysis ----

# Clinical enrichment analysis

pdf(file = "maftools_out_byState.pdf", width = 150, height = 20,)

clinicalEnrichment = clinicalEnrichment(maf = laml, clinicalFeature = 'State')
par(oma=c(2,2,2,2))
plotEnrichmentResults(enrich_res = clinicalEnrichment , pVal = 0.05, geneFontSize = 1, annoFontSize = 1)

dev.off()


pdf(file = "maftools_out_byPCR.pdf")

clinicalEnrichment.2 = clinicalEnrichment(maf = laml, clinicalFeature = 'cAst_PCR')

par(oma=c(2,2,2,2))
plotEnrichmentResults(enrich_res = clinicalEnrichment.2, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)

dev.off()



# Compare the aWM_stable VS aWM_progressed

pdf(file = "maftools_aWM_stableVSprogressed.pdf")

##Subset maf and return output as an MAF object (Default behavior)
aWM_stable.maf <- subsetMaf(maf = laml, tsb = aWM_stable)
aWM_progressed.maf <- subsetMaf(maf = laml, tsb = aWM_progressed)

#Considering only genes which are mutated in at-least in 3 samples in one of the cohort to avoid bias due to genes mutated in single sample.
pt.vs.rt <- mafCompare(m1 = aWM_stable.maf, m2 = aWM_progressed.maf, m1Name = 'aWM_stable', m2Name = 'aWM_progressed', minMut = 3)
print(pt.vs.rt)

forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, geneFontSize = 0.8)



genes = c("CACNA1A", "PKHD1", "CD79B", "DNM1P47", "NBPF1", "IGLL5", "CDH6", "CELF2", "MUC2","MUC4",
          "ADCY1", "KIR3DL1", "KIR3DL3", "KIR3DX1", "LILRB3")


coOncoplot(m1 = aWM_stable.maf, m2 = aWM_progressed.maf, m1Name = 'aWM_stable', m2Name = 'aWM_progressed', removeNonMutated = TRUE, showSampleNames = TRUE, genes = genes,  
           geneNamefont = 0.4, legendFontSize = 0.5, titleFontSize = 0.8)

coBarplot(m1 = aWM_stable.maf, m2 = aWM_progressed.maf, m1Name = "aWM_stable", m2Name = "aWM_progressed", genes = genes)

dev.off()




# Mutational Signatures ----
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)

#  I work only with the SNVs here!!!
laml <- merge_mafs(c(mafs.read),
                   vc_nonSyn = c(nonSilentPLUSsilent)
                   )

laml.tnm = trinucleotideMatrix(maf = laml, prefix = NULL, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

# Estimate APOBEC enrichment scores
# APOBEC induced mutations are more frequent in solid tumors, associated with C>T transition events
# Differences between APOBEC enriched and non-enriched samples
plotApobecDiff(tnm = laml.tnm, maf = laml, pVal = 0.2)
# ---APOBEC related mutations are enriched in  0 % of samples (APOBEC enrichment score > 2 ;  0  of  67  samples

# Signature analysis
library('NMF')
# runs NMF on a range of values and measures the goodness of fit - in terms of Cophenetic correlation.
laml.sign = estimateSignatures(mat = laml.tnm, nTry = 10)

# # elbow plot to visualize and decide optimal number of signatures
plotCophenetic(res = laml.sign)

# set the optimal number of signatures
laml.sig = extractSignatures(mat = laml.tnm, n = 6)



#Compate against original 30 signatures
laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")

#Compate against updated version3 60 signatures
laml.v3.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "SBS")


# comparison of similarities of detected signatures against validated signatures
library('pheatmap')
pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
pheatmap::pheatmap(mat = laml.v3.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")


# plot signatures
maftools::plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = "legacy")
maftools::plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = "SBS")




# to plot the signature contribution per patient
maftools::plotSignatures(nmfRes = laml.sig, contributions = TRUE, title_size = 0.05, sig_db = "legacy")
# par(oma=c(1,1,1,1)) 
maftools::plotSignatures(nmfRes = laml.sig, contributions = TRUE, title_size = 0.05, font_size = 0.6, sig_db = "SBS",
                         patient_order = patient.order,
                         show_barcodes = TRUE, show_title = TRUE,
                         legend.text = c("A", "B"),
                         args.legend = list(x = "topleft")
                         )





# Density plots, tumor heterogeneity and MATH scores ----

# Prepare the segmentation file (must be tab separated file):
# Column names should be: Sample, Chromosome, Start, End, Num_Probes and Segment_Mean (log2 scale)
# !! Num_Probes is not important for the analysis !!

segFile.df <- patients.cnvs.df[, c("patient_ID", "chr", "start", "end", "num_mark", "cnlr_median")]

colnames(segFile.df) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")

write.table(segFile.df, file = "segFile.txt" , quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")



# "T_1558","T_1919","T_2029","T_2065","T_2088","T_2090","T_2245","T_2298","T_2303","T_2306","T_3100","T_3316",T_3540","T_3546","T_3672","T_3721","T_3826","T_3860","T_3905","T_4010","T_4021","T_4162","T_4198","T_4250","T_4269","T_4272","T_4302","T_4366","T_4388","T_4445","T_4496","T_4511","T1_2308","T_2308_2"

samplename <- "T_2308_2"

pdf(file = paste0("./density_plots/", samplename, ".pdf"))

het1 = inferHeterogeneity(maf = laml, tsb = samplename, vafCol = NULL, useSyn = TRUE)

# ignore variants located on copy-number altered regions
het2 = inferHeterogeneity(maf = laml, tsb = samplename, vafCol = NULL, useSyn = TRUE, segFile = "segFile.txt")

print(het1$clusterMeans)

plotClusters(clusters = het1, genes = "all")

print(het2$clusterMeans)

plotClusters(clusters = het2, genes = 'CN_altered', showCNvars = TRUE)

dev.off()


# # Plots ignoring chromosomes with CNVs
# het = inferHeterogeneity(maf = laml, tsb = "T_1962", vafCol = NULL, ignChr = c("chr2","chr6","chr14","chr18","chrX","chrY"), useSyn = TRUE)
# het = inferHeterogeneity(maf = laml, tsb = "T_1558", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_1919", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_2029", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_2088", vafCol = NULL, ignChr = c("chr1","chr2","chr7","chr8","chr9","chr11","chr14","chr16","chr17","chr19","chr20","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_2090", vafCol = NULL, ignChr = c("chr7","chr9","chr12","chr19","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_2245", vafCol = NULL, ignChr = c("chr1","chr4","chr7","chr12","chr14","chr16","chr19","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_2298", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_2303", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_2306", vafCol = NULL, ignChr = c("chr7","chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_3100", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_3100", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_3316", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_3540", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_3546", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_3672", vafCol = NULL, ignChr = c("chr1","chr2","chr5","chr6","chr14","chr16","chr19","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_3721", vafCol = NULL, ignChr = c("chr1","chr2","chr9","chr14","chr15","chr16","chr18","chr19","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_3826", vafCol = NULL, ignChr = c("chr1","chr2","chr5","chr7","chr11","chr14","chr15","chr16","chr19","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_3860", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_3905", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_3967", vafCol = NULL, ignChr = c("chr1","chr5","chr15","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T_4010", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T1_2308", vafCol = NULL, ignChr = c("chr1","chr6","chr12","chr14","chr19","chrX","chrY"))
# het = inferHeterogeneity(maf = laml, tsb = "T2_3164", vafCol = NULL, ignChr = c("chr14","chrX","chrY"))





