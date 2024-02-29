library("maftools")
library("data.table")
library("tidyverse")
library("ggplot2")


somatic_snvs_maf_wgs_14samples <- "/home/rania/Documents/Tina/WES/Fra_SNVs/maf/"


## Variant_Classifications
nonSynonymous = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                  "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")


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




# SNVs from WGS
setwd(somatic_snvs_maf_wgs_14samples)
somatic_snvs_wgs <- list.files(path = somatic_snvs_maf_wgs_14samples)
maf_list_somatic_snvs_wgs <- lapply(somatic_snvs_wgs, function(x) data.table::fread(x))


maf_somatic_snvs_wgs <- rbindlist(c(maf_list_somatic_snvs_wgs, maf_list_somatic_snvs_wgs))

laml <- read.maf(maf_somatic_snvs_wgs,
                 vc_nonSyn = nonSynonymous,
                 clinicalData = clinical.info)





pdf(file = "maftools_out.pdf")

# plot the summary of the maf file: number of variants in each sample as a stacked barplot and variant types as a boxplot summarized by Variant_Classification.
plotmafSummary(maf = laml, addStat = 'median', titvRaw = FALSE, top = 25, fs = 0.8, showBarcodes = TRUE, textSize = 0.3)

# oncoplot for top mutated genes
oncoplot(maf = laml,
         top = 50,
         fontSize = 0.3, annotationFontSize = 1, legendFontSize = 1, titleFontSize = 0.8,
         showTumorSampleBarcodes = TRUE,
         SampleNamefontSize = 0.7,
         clinicalFeatures = c("cAst_PCR", "State"),
         sortByAnnotation = TRUE,   ##sorts only according to the 1st clinical feature
         #groupAnnotationBySize = FALSE,
         #annotationOrder = c("IGM-MGUS", "aWM_stable", "aWM_progressed", "sWM", "treated"),
         includeColBarCN = FALSE
)

# oncoplot with genes Tina had asked me about
oncoplot(maf = laml,
         genes = c("MYD88","CXCR4","TP53","ARID1A","ARID1B","CD79A","CD79B","KMT2D","KDM6A","KMT2C","NFKB1","NFKB2","NFKBIZ","BCL2","BCL6","BCL10","ATM","NOTCH1","NOTCH2","NOTCH3","TRAP","HIST1H1D","EZH2","EGF","TBL1XR1"),
         fontSize = 0.3, annotationFontSize = 1, legendFontSize = 1, titleFontSize = 0.8,
         showTumorSampleBarcodes = TRUE,
         clinicalFeatures = c("State","cAst_PCR"),
         sortByAnnotation = TRUE,   ##sorts only according to the 1st clinical feature
         groupAnnotationBySize = FALSE,
         annotationOrder = c("IGM-MGUS", "aWM_stable", "aWM_progressed", "sWM", "treated"),
         includeColBarCN = FALSE
)

dev.off()





