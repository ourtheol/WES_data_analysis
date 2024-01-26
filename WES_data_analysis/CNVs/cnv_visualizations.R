# Plot the CNVs infered by CNV_FACETS, according to their locations on the short/long arms of each chromosome 


library("GenomicRanges")

library("vcfR")

library("biscuiteer")

library("dplyr")
library("tidyr") # pivot_wider
library("tidyverse") 

library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")


# Colornames for the plot can be found here:
# https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf


# path to the vcfs output from cnv_facets
cnv.vcf.files.dir <- "/home/rania/Documents/Tina/WES/single_samples/cnv_vcfs/"

# hg38 chromosome arms ranges
data(GRCh38.chromArm, package="biscuiteer")

# Add a meta-data column
GRCh38.chromArm$chrom.arm <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13p","13q","14p","14q",
                               "1p5","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21p","21q","22p","22q","Xp","Xq","Yp","Yq","MT")



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
  
  # data frame with the CNVs called by cnv_facets
  cnv.df <- data.frame(chrom, pos, end, svtype, tcn)
  
  return(cnv.df)
  
}



# function: vec2list -----------------------------------------------------------
# insert a vector, and return a list with the elements divided by ";"

vec2list <- function(x) {
  
  paste( unlist(x), collapse=';' )
  
}



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


# # The CNV types detected:
# unique(patients.cnvs.df$svtype)
# # "NEUTR"   "HEMIZYG" "DEL"     "DUP"     "LOH"     "DUP-LOH"


# Rename the column names to the exprected from "makeGRangesListFromDataFrame" function
colnames(patients.cnvs.df) <- c("chr", "start", "end", "svtype", "state", "patient_ID")

# Read the R data-frame and group the calls by patient ID and convert them to a GRangesList
grl <- GenomicRanges::makeGRangesListFromDataFrame(patients.cnvs.df, 
                                                   split.field="patient_ID", keep.extra.columns=TRUE)



# Loop throught the GRangesList and map each cnv to the long/short arm of a chrom
for (i in 1:length(unique(patients.cnvs.df$patient_ID))){
  
  print(i)

  # Match genomic locations with chromosome arm information
  matching_rows <- findOverlaps(grl[[i]], GRCh38.chromArm, select = "arbitrary") # choose this if a cnv spans to multiple chrom arms
  grl[[i]]$chrom.arm <- mcols(GRCh38.chromArm)$chrom.arm[matching_rows]
  
  
}




# Plotting time

# GRanges list into data-frame:
grl.df <- as.data.frame(grl)

# Turn into wide table and insert the chromosome arms as rownames
wide <- pivot_wider(grl.df[c("group_name", "svtype", "chrom.arm")], names_from = group_name, values_from = svtype, names_sep = ";" ) %>% 
  column_to_rownames(var="chrom.arm") %>% 
  replace(.=="NULL", "") # Replace NULL with blank


# Handling chrom.arms with multiple entries within the same patient
# multiple entries are as vectors, turn them into strings
wide <- apply(wide, c(1,2),  vec2list )


# Remove the X chromosome
wide <- wide[!(row.names(wide) %in% c("Xp","Xq")),]



# Read file with the annotations for each patient
annotation.table <- read.table("/home/rania/Documents/Tina/WES/single_samples/annotations.txt", row.names = 1)

# Reorder the patient_IDs according to the input matrix (wide) in the oncoPrint()
annotation.table$State <- factor(annotation.table$State, levels = c("IGM-MGUS", "aWM_stable", "aWM_progressed", "sWM", "treated"))
annotation.table.reordered <- annotation.table[order(annotation.table$State, annotation.table$cAst_PCR),]
patient.order <- rownames(annotation.table.reordered)

# Prep annotations: State, PCR for the oncoplot
anno.df <- annotation.table[colnames(wide), c("State", "cAst_PCR"), drop = FALSE]
mix.ha <- HeatmapAnnotation(df = anno.df, 
                           col = list(State = c("IGM-MGUS" = "khaki3", "aWM_stable" = "salmon1", "aWM_progressed" = "indianred2", "sWM" = "orangered3", "treated" = "chartreuse4"),
                                      cAst_PCR = c("wt" = "mediumpurple2", "mut" = "mediumaquamarine")))





# oncoPrint
col = c(DUP = "orangered2", HEMIZYG = "skyblue3", DEL = "blue", NEUTR = "yellow2", "LOH" = "orange", "DUP-LOH" = "sienna3")


oncoht <- oncoPrint(wide,
                    
                    alter_fun = function(x, y, w, h, v) {
            
                    n = sum(v)  # how many alterations for current gene in current sample
                    h = h*0.9
                    # use `names(which(v))` to correctly map between `v` and `col`
                    if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h, 
                                    gp = gpar(fill = col[names(which(v))], col = NA), just = "top")
                    },
                    
                    col = col,
                    show_column_names = TRUE,
                    column_names_rot = 90,
                    bottom_annotation = mix.ha,
                    
                    column_order = patient.order
                    ) 





################################################################################

library("CNVRanger")
library("Gviz")

library("AnnotationHub")


# Summarize CNV calls to CNV regions
cnvrs <- populationRanges(grl, density = 0.1)


# Distinguish driver from passenger mutations, i.e. to distinguish meaningful events from random background aberrations. 
# The GISTIC method identifies those regions of the genome that are aberrant more often than would be expected by chance, 
# with greater weight given to high amplitude events (high-level copy-number gains or homozygous deletions) that are less likely to represent random aberrations 
cnvrs <- populationRanges(grl, density = 0.1, est.recur = TRUE)
# and filter for recurrent CNVs that exceed a significance threshold of 0.05.
subset(cnvrs, pvalue < 0.05)


# work independently for each group:
# Read the R data-frame and group the calls by patient ID and convert them to a GRangesList
grl.IGM_MGUS <- GenomicRanges::makeGRangesListFromDataFrame(patients.cnvs.df[patients.cnvs.df$patient_ID %in% IGM_MGUS, ],
                                                            split.field="patient_ID", keep.extra.columns=TRUE)
cnvrs.IGM_MGUS <- populationRanges(grl.IGM_MGUS, density = 0.1, est.recur = TRUE)
subset(cnvrs.IGM_MGUS, pvalue < 0.05)


grl.aWM_stable <- GenomicRanges::makeGRangesListFromDataFrame(patients.cnvs.df[patients.cnvs.df$patient_ID %in% aWM_stable, ],
                                                              split.field="patient_ID", keep.extra.columns=TRUE)
cnvrs.aWM_stable <- populationRanges(grl.aWM_stable, density = 0.1, est.recur = TRUE)
subset(cnvrs.aWM_stable, pvalue < 0.05)


grl.aWM_progressed <- GenomicRanges::makeGRangesListFromDataFrame(patients.cnvs.df[patients.cnvs.df$patient_ID %in% aWM_progressed, ],
                                                                  split.field="patient_ID", keep.extra.columns=TRUE)
cnvrs.aWM_progressed <- populationRanges(grl.aWM_progressed, density = 0.1, est.recur = TRUE)
subset(cnvrs.aWM_progressed, pvalue < 0.05)


grl.sWM <- GenomicRanges::makeGRangesListFromDataFrame(patients.cnvs.df[patients.cnvs.df$patient_ID %in% sWM, ],
                                                       split.field="patient_ID", keep.extra.columns=TRUE)
cnvrs.sWM <- populationRanges(grl.sWM, density = 0.1, est.recur = TRUE)
subset(cnvrs.sWM, pvalue < 0.05)



# plotRecurrentRegions(cnvrs, genome = "hg38", chr = "chr14")


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



# chromosome selection
sel.genes <- subset(genes, seqnames == "16")

olaps <- GenomicRanges::findOverlaps(sel.genes, cnvrs, ignore.strand=TRUE)
qh <- S4Vectors::queryHits(olaps)
sh <- S4Vectors::subjectHits(olaps)
cgenes <- sel.genes[qh]
# annotate the CNV type (gain/loss) for genes overlapping with CNVs
cgenes$type <- cnvrs$type[sh]
subset(cgenes, select = "type")

# cnvOncoPrint plot
# the input matrix and the annotation matrix must be in the same order, for the annotations to not get mixed up
# cnvs.chr@column_names_param$labels is the order the sampleIDs are arranged in the plot

# since the samples are reordered first do the plot and then rerun to add the annotations



cnvs.chr <- cnvOncoPrint(grl, 
             cgenes,
             show_column_names = TRUE,
             top.features = 80,
             column_names_rot = 90)



anno.reordered <- anno.df[order(match(rownames(anno.df), cnvs.chr@column_names_param$labels)), , drop = FALSE]

cnvs.chr <- cnvOncoPrint(grl, 
                         cgenes,
                         show_column_names = TRUE,
                         top.features = 80,
                         column_names_rot = 90,
                         bottom_annotation = HeatmapAnnotation(df = anno.reordered,
                                                               col = list(State = c("IGM-MGUS" = "khaki3", "aWM_stable" = "salmon1", "aWM_progressed" = "indianred2", "sWM" = "orangered3", "treated" = "chartreuse4"),
                                                                          cAst_PCR = c("wt" = "mediumpurple2", "mut" = "mediumaquamarine")))
)



png("/home/rania/Documents/Tina/WES/cnv_plots/genes_affected_by_CNVs/chr16.png", width = 1800, height = 1200)
print(cnvs.chr)
dev.off()





################################################################################

# Plot percentage of the cohort having alterations on the selected genomic regions

library("GenVisR")

cnSpec.input <- subset(patients.cnvs.df, select = c("chr", "start", "end", "state", "patient_ID"))
colnames(cnSpec.input) <- c("chromosome", "start", "end", "segmean" ,"sample")


cnFreq(cnSpec.input,
       genome="hg38", 
       plotChr=c(paste0("chr",rep(1:22))))


# plot the number of patients
cnFreq(cnSpec.input,
       genome="hg38", 
       plotChr=c(paste0("chr",rep(1:22))),
       plotType = "frequency")
