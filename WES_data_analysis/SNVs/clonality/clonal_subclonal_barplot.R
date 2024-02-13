library(dplyr)
library(tidyr)
library(ggplot2)
library("cowplot")
# library("xlsx")


# the genes I have selected
genes <- c("MYD88", "CXCR4", "CD79B", "IGLL5", "MUC2", "MUC5AC", "KMT2C", "PKHD1", "CACNA1A", "FOXD4L1")


sample_purity <- c(
"T_1546"=0.92,"T_1558"=0.17,"T_1735"=0.70,"T_1781"=0.78,"T_1831"=0.72,"T_1840"=0.50,"T_1873"=0.28,"T_1876"=0.20,"T_1919"=0.14,"T_1924"=0.45,
"T_1962"=0.80,"T_2029"=0.17,"T_2048"=0.3,"T_2065"=0.15,"T_2088"=0.16,"T_2090"=0.2,"T_2102"=0.20,"T_2126"=0.63,"T_2142"=0.30,"T_2222"=0.90,
"T_2245"=0.16,"T_2251"=0.32,"T_2257"=0.82,"T_2298"=0.12,"T_2299"=0.22,"T_2300"=0.40,"T_2303"=0.34,"T_2306"=0.19,"T_2316"=0.30,"T_2359"=0.22,
"T_2501"=0.56,"T_3100"=0.32,"T_3316"=0.22,"T_3323"=0.32,"T_3373"=0.4,"T_3400"=0.49,"T_3419"=0.40,"T_3540"=0.2,"T_3546"=0.16,"T_3597"=0.54,
"T_3672"=0.4,"T_3721"=0.22,"T_3826"=0.18,"T_3829"=0.50,"T_3860"=0.16,"T_3905"=0.14,"T_3967"=0.06,"T_4010"=0.3,"T_4021"=0.22,"T_4153"=0.12,
"T_4162"=0.18,"T_4198"=0.26,"T_4250"=0.16,"T_4269"=0.12,"T_4272"=0.46,"T_4302"=0.12,"T_4366"=0.15,"T_4384"=0.26,"T_4388"=0.08,"T_4445"=0.2,
"T_4496"=0.17,"T_4511"=0.18,"T_4514"=0.64,"T_4528"=0.32,"T_4641"=0.49,"T_721"=0.26,"T_800"=0.54,"T2_3164"=0.39,"T1_608"=0.58,"T_608"=0.26,
"T_650_1"=0.63,"T_650_2"=0.71,"T_650_3"=0.52,"T_655_1"=0.84,"T_655_2"=0.76,"T_790_1"=0.44,"T_790_2"=0.20,"T1_805"=0.72,"T2_805"=0.70,
"T1_2308"=0.15,"T_2308_2"=0.25,"T_1928"=0.76,"T2_1928"=0.76
)                                              


patient.gene.vaf <- read.csv("/home/rania/Documents/Tina/WES/clonal_subclonal/snv_mafs_combined.tsv", sep=",", header=TRUE)
View(patient.gene.vaf)


# Calculate t_vaf, ccf, set to clonal/subclonal according to ccf
patient.gene.vaf <- patient.gene.vaf %>% 
  dplyr::mutate(t_vaf = t_alt_count / (t_alt_count + t_ref_count)) %>% 
  dplyr::mutate(purity = recode(Tumor_Sample_Barcode, !!!`sample_purity`)) %>% 
  dplyr::mutate(ccf = t_vaf / purity) %>%
  dplyr::mutate(variant_clonality = if_else(ccf < 0.7 , "subclonal", "clonal")) %>%
  # count the number of variants (in all patients) per gene and put the genes in descending order
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::arrange(desc(count), .by_group=FALSE)

# Subset only the genes with more than 5 variants across our samples
patient.gene.vaf <- patient.gene.vaf[patient.gene.vaf$count >= 4,]                     
View(patient.gene.vaf)                     


# make the gene names an ordered factor, so that they do not get alphabetically ordered in the plot
patient.gene.vaf$Hugo_Symbol <- factor(patient.gene.vaf$Hugo_Symbol, levels = unique(patient.gene.vaf$Hugo_Symbol))

dotplot <- ggplot(patient.gene.vaf, aes(x=Hugo_Symbol, y=ccf)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="red") +
  xlab("") +
  ylab("CCF")



View(patient.gene.vaf[patient.gene.vaf$ccf > 1,])
View(patient.gene.vaf[patient.gene.vaf$Hugo_Symbol == "MYD88",])
View(patient.gene.vaf[patient.gene.vaf$Hugo_Symbol == "CXCR4",])

write.xlsx(patient.gene.vaf[patient.gene.vaf$ccf > 1,], file="clonal_subclonal.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = FALSE, append = FALSE)
write.xlsx(patient.gengene_liste.vaf[patient.gene.vaf$Hugo_Symbol == "MYD88",], file="MYD88.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = FALSE, append = FALSE)


# Count the number of clonal and subclonal variants for each gene
plot.prep <- subset(patient.gene.vaf, select = c(Hugo_Symbol, variant_clonality)) %>% 
  dplyr::group_by(Hugo_Symbol, variant_clonality) %>% 
  dplyr::summarise(total_count=n(), .groups = 'drop') 


# Stacked barplot with multiple groups
barplot <- ggplot(data=plot.prep, aes(x=Hugo_Symbol, y=total_count, fill=variant_clonality)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(legend.position="none") +
  ylab("variant count")
  
# Stacked barplot with multiple groups
barplot2 <- ggplot(data=plot.prep, aes(x=Hugo_Symbol, y=total_count, fill=variant_clonality)) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position='top') +
  xlab("") +
  ylab("%")
  
# arrange the two plots into one column
pdf("clonal_subclonal.pdf")
plot_grid(
  dotplot, barplot2, barplot, 
  labels = "AUTO", ncol = 1)
dev.off()






################################################################################
# Work again the same way with correcting for copy number
patient.gene.vaf <- read.csv("/home/rania/Documents/Tina/WES/clonal_subclonal/snv_mafs_combined.tsv", sep=",", header=TRUE)

# Subset only part of the columns
patient.gene.vaf <- patient.gene.vaf[, c("Hugo_Symbol","Tumor_Sample_Barcode","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele1",
                                         "Tumor_Seq_Allele2","t_ref_count","t_alt_count","t_depth",
                                         "copy_num_status","minor_copy_num_status","sv_type","cellular_fraction")] 
View(patient.gene.vaf)


# Calculate t_vaf, ccf, set to clonal/subclonal according to ccf
patient.gene.vaf <- patient.gene.vaf %>% 
  dplyr::mutate(t_vaf = t_alt_count / (t_alt_count + t_ref_count)) %>% 
  dplyr::mutate(purity = recode(Tumor_Sample_Barcode, !!!sample_purity)) %>% 
  # correcting for copy number
  dplyr::mutate(ccf = t_vaf * copy_num_status / purity) %>%   
  dplyr::mutate(variant_clonality = if_else(ccf < 0.75 , "subclonal", "clonal")) %>%
  # count the number of variants (in all patients) per gene and put the genes in descending order
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::arrange(desc(count), .by_group=FALSE)


# Subset only the entries with the aforementioned genes of interest
patient.gene.vaf <- patient.gene.vaf[patient.gene.vaf$Hugo_Symbol %in% genes,]
# # Subset only the genes with more than 5 variants across our samples
# patient.gene.vaf <- patient.gene.vaf[patient.gene.vaf$count >= 4,]
# View(patient.gene.vaf)



View(patient.gene.vaf[patient.gene.vaf$Hugo_Symbol == "MYD88",])

View(patient.gene.vaf[patient.gene.vaf$Hugo_Symbol == "CXCR4",])

View(patient.gene.vaf[patient.gene.vaf$ccf > 1,])
View(patient.gene.vaf[patient.gene.vaf$copy_num_status == 3,])





# make the gene names an ordered factor, so that they do not get alphabetically ordered in the plot
patient.gene.vaf$Hugo_Symbol <- factor(patient.gene.vaf$Hugo_Symbol, levels = unique(patient.gene.vaf$Hugo_Symbol))

dotplot <- ggplot(patient.gene.vaf, aes(x=Hugo_Symbol, y=ccf)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="red") +
  xlab("") +
  ylab("CCF")


# Count the number of clonal and subclonal variants for each gene
plot.prep <- subset(patient.gene.vaf, select = c(Hugo_Symbol, variant_clonality)) %>% 
  dplyr::group_by(Hugo_Symbol, variant_clonality) %>% 
  dplyr::summarise(total_count=n(), .groups = 'drop') 


# Stacked barplot with multiple groups
barplot <- ggplot(data=plot.prep, aes(x=Hugo_Symbol, y=total_count, fill=variant_clonality)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(legend.position="none") +
  ylab("variant count")

# Stacked barplot with multiple groups
barplot2 <- ggplot(data=plot.prep, aes(x=Hugo_Symbol, y=total_count, fill=variant_clonality)) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position='top') +
  xlab("") +
  ylab("%")

# arrange the two plots into one column
pdf("clonal_subclonal.pdf")
plot_grid(
  dotplot, barplot2, barplot, 
  labels = "AUTO", ncol = 1)
dev.off()




# # Alternative plotting option:
# # long data to wide data
# plot.prep.long <- tidyr::spread(plot.prep,  
#                                 key = Hugo_Symbol,
#                                 value = total_count,
#                                 fill = 0)
# 
# # Set variant_clonality into rownames
# plot.prep.long <- data.frame(plot.prep.long, row.names = 1)
# # Into numeric
# plot.prep.long <- sapply(plot.prep.long, as.numeric)
# 
# 
# 
# # Get the stacked barplot
# barplot <- barplot(plot.prep.long, 
#                    col=colors()[c(23,89)] , 
#                    border="white",
#                    space=0.04, 
#                    font.axis=2, 
#                    xlab="gene") 

