

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# setRepositories(ind=c(1:6))
# BiocManager::install("biomaRt")
# BiocManager::install("tximport", force=TRUE)
# BiocManager::install("DESeq2", force=TRUE)
# BiocManager::install('rhdf5')
# BiocManager::install('this.path')
# BiocManager::install('ramify')
# BiocManager::install('crayon')

library(biomaRt)        # required to map transcripts to genes
library(tximport)       # required to read input files
library(DESeq2)         # the library that will call DEGs
library(crayon)         # so the messages are blue
library(this.path)      # necessary to locate where this file is
library(ggplot2)        # useful for plotting
library(ramify)         # necessary for the clip function
library(rhdf5)          # necessary for reading the input files

#
# 0. user-defined variables
#

# set your working directory. This is an option, but you can change as you prefer.
script_path <- "C:/Users/thorb/Documents/OfvirkiFiskurinn/RNASequencinganalysis/020_vik"
script_path
setwd(script_path) 

kallisto_dir = "output_alt"
results_dir = 'NEW_DEGs_DESeq2'

effect_size_threshold = log2(2) # arbitrary: we will discard DEGs that hold less than abs FC = 2
count_threshold = 20 # arbitrary threshold: we will discard DEGs that hold less than 20 reads difference
tpm_threshold = 1

# 
# 1. get todays working data: kallisto output from two conditions
#

# List all files in output directory
list.files('output_alt')
list.files('output_alt/oligo01/')  # Example for checking oligo01

# Read the first abundance file as an example
df = read.csv('output_alt/oligo01/abundance.tsv', sep='\t')
View(df)

# Define WT and Adgrl3.1-/- groups
wt_folders <- c("oligo01", "oligo02", "oligo03", "oligo04")
adgrl3_folders <- c("oligo05", "oligo06", "oligo07", "oligo08")

# 2. Generate gene-to-transcript mapping for **zebrafish**
mart = biomaRt::useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "drerio_gene_ensembl",  # Change from human to zebrafish
  host = "https://www.ensembl.org",
  verbose = TRUE
)

# Define the attributes for transcript-to-gene mapping
working_attributes = c(
  'ensembl_transcript_id',
  'ensembl_gene_id',
  'external_gene_name',
  'gene_biotype',
  'description'
)

# Retrieve gene-transcript mapping
t2g = biomaRt::getBM(
  attributes = working_attributes,
  mart = mart,
  verbose = TRUE
)

# Check dimensions of the retrieved table
dim(t2g)

# Get directories but filter out oligo09 to oligo12
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames = dirnames[grepl('oligo0[1-8]$', dirnames)]  # Keep only oligo01 to oligo08

# Generate paths to the abundance.h5 files
paths = file.path(dirnames, 'abundance.tsv')

# Extract labels (folder names)
labels = sapply(strsplit(paths, split='/', fixed=TRUE), function(x) x[2])

# Define replicates and treatments
# Assuming you have 4 samples for each condition (WT and Adgrl3.1-/-)
replicates = rep(c('A', 'B', 'C', 'D'), 2)  # 4 replicates for each condition
treatments = rep(c('WT', 'WT', 'WT', 'WT', 'MUTANT', 'MUTANT', 'MUTANT', 'MUTANT'))

# Create metadata dataframe
metadata = data.frame(labels)
metadata$replicate = replicates
metadata$treatment = treatments
metadata$path = paths

# View the resulting metadata
View(metadata)

#
# 4. determine significance of change on treated vs non-treated
#

# Let's ask DESeq2 to read the files
files = metadata$path  # Using paths from metadata
print(files)

# Import the data using tximport, specify zebrafish tx2gene mapping
txi = tximport(files, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

# DESeqDataSet creation
dds = DESeqDataSetFromTximport(txi, colData=metadata, design=~treatment)

# Relevel the treatment factor so that "WT" is the reference level
dds$treatment = relevel(dds$treatment, ref="WT")

# Keep features with at least 20 counts median difference
cat(blue(paste('Size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
a = counts(dds)[, 1:4]  # Adjusted for 4 replicates in each treatment group
b = counts(dds)[, 5:8]
c = rowMedians(a) - rowMedians(b)
keep = abs(c) >= count_threshold
dds = dds[keep, ]
cat(blue(paste('Size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

# Keep features with at least a max median expression of 1 TPM
cat(blue(paste('Size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
subset = txi$abundance[names(dds), ]
a = rowMedians(subset[, 1:4])  # Adjusted for 4 replicates per treatment
b = rowMedians(subset[, 5:8])
c = pmax(a, b)
keep = c >= tpm_threshold
dds = dds[keep, ]
cat(blue(paste('Size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

# Run the test. Here you can choose between LRT or Wald test
dds = DESeq(dds, test="LRT", reduced=~1)

# Retrieve the results
res = results(dds, parallel=TRUE, alpha=0.05)

# Filter results with a threshold of adjusted p-value < 0.05 and log2FoldChange > effect_size_threshold
filtered_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtered_results = filtered_results[order(filtered_results[["padj"]]),]

# Results with no significant change
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]

cat(blue(paste('Contrast WT vs MUTANT:', dim(filtered_results)[1], sep=' ')), fill=TRUE)

# Store results in the designated directory
dir.create(results_dir)
write.table(sorted_filtered_results, file=paste(results_dir, '/NEW_WT_vs_MUTANT_filtered.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, file=paste(results_dir, '/NEW_WT_vs_MUTANT_anti.tsv', sep=''), quote=FALSE, sep='\t')


#
# 5. visualization
#

# a simple PCA
plotPCA(rlog(dds), intgroup=c('treatment')) + ggtitle('effect WT vs MUTANT')

plotPCA(rlog(dds), intgroup=c('treatment')) + 
  ggtitle('Effect WT vs MUTANT') +
  theme_minimal() +  # Makes the background white
  scale_color_manual(values = c("#33ff99", "#cc33ff")) +  # Custom colors for groups
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centered title
    legend.position = "right",  # Moves the legend to the right
    axis.title = element_text(size = 12),  # Adjusts axis label size
    axis.text = element_text(size = 10)  # Adjusts axis text size
  )

# Ensure your DESeq2 results are in a dataframe
res_df <- as.data.frame(res)  # Convert DESeq2 results to dataframe
res_df$ensembl_gene_id <- rownames(res_df)  # Make row names a column
  
NOC_df <- as.data.frame(resNOC)  # Convert DESeq2 results to dataframe
NOC_df$ensembl_gene_id <- rownames(NOC_df)  # Make row names a column

# Merge with t2g to add gene name, gene biotype, and description

# get rid of duplicates in t2g before merging

t2g_unique <- t2g[!duplicated(t2g$ensembl_gene_id), ]

NOC_annotated <- merge(NOC_df, 
                       t2g_unique[, c("ensembl_gene_id", "external_gene_name")], 
                       by = "ensembl_gene_id", all.x = TRUE)

# View the updated dataframe
head(res_annotated)

# Ensure no zero padj values (to avoid log10 issues)
res_annotated$padj[res_annotated$padj == 0] <- min(res_annotated$padj[res_annotated$padj > 0], na.rm = TRUE)

# Define plotting variables
res_annotated$logP <- -log10(res_annotated$padj)

# Set threshold for significance
padj_threshold <- 0.05
logFC_threshold <- 1

# Identify significant genes
res_annotated$significance <- "Not Significant"
res_annotated$significance[res_annotated$padj < padj_threshold & res_annotated$log2FoldChange > logFC_threshold] <- "Upregulated"
res_annotated$significance[res_annotated$padj < padj_threshold & res_annotated$log2FoldChange < -logFC_threshold] <- "Downregulated"

# Identify the top 20 most significant genes
top20 <- res_annotated[order(-res_annotated$logP), ][1:20, ]

# Save the new annotaded as a file
dir.create(results_dir)
write.table(res_annotated, file=paste(results_dir, '/NEW_WT_vs_MUTANT_annotated_update.tsv', sep=''), quote=FALSE, sep='\t')



# Create volcano plot VOLCANO PLOT VOLCANO PLOT VOLCANO PLOT VOLCANO PLOT

library(ggplot2)
install.packages("ggplot2", dependencies = TRUE)

install.packages("ggrepel")  # if not already installed
library(ggrepel)



ggplot(res_annotated, aes(x = log2FoldChange, y = logP, color = significance)) +
  geom_point(alpha = 0.8) +  # Make all points visible but faded
  scale_color_manual(values = c("Upregulated" = "#cc33ff", "Downregulated" = "#33ff99", "Not Significant" = "grey")) + 
  geom_text_repel(data = top20, aes(label = external_gene_name), 
                  size = 3, color = "black", # Set text color to black
                  box.padding = 1, point.padding = 1, # Increase padding to reduce overlap
                  max.overlaps = Inf, force = 4) + # Increase force to push labels apart
  labs(x = "log2 Fold Change", y = "-log10 Adjusted P-Value", title = "Differential Gene Expression: Mutant vs. Wildtype") +
  theme_classic()




# Thresholds
padj_threshold <- 0.05
logFC_threshold <- 1

# Categorize genes
res_annotated$Significance <- "Not Significant"
res_annotated$Significance[res_annotated$padj < padj_threshold & res_annotated$log2FoldChange > logFC_threshold] <- "Upregulated"
res_annotated$Significance[res_annotated$padj < padj_threshold & res_annotated$log2FoldChange < -logFC_threshold] <- "Downregulated"



library(dplyr)

# Define the genes of interest
genes_of_inflammatory_interest <- c("il6st", "il6r", "aif1l", "tnfrsf1a", "tnfrsf11b", "il10rb", "il1rapl2", "il15l", "ifngr1l", "cxcl12b", "cxcl12a", "icam3", "il1b", "il1r1", "il1rap", "tnfa", "tnfrsf1a", "tnfrsf1b", "il6", "il6r", "il6st", "cxcl8a", "cxcr1", "cxcr2", "ccL34a.4", "ccr1", "ccL35.1", "ccr5", "tgfb1a", "tgfbr1a", "tgfbr2a", "ifng1-1", "ifngr1", "ifngr2", "bdnf", "ntrk2a", "ntrk2b", "cd86", "cd28", "il10", "il10ra", "il10rb")

# Filter dataset
CNS_inflammation <- res_annotated %>% 
  filter(external_gene_name %in% genes_of_inflammatory_interest)

# Categorize genes
CNS_inflammation$Significance <- "Not Significant"
CNS_inflammation$Significance[CNS_inflammation$padj < 0.05 & CNS_inflammation$log2FoldChange > 0.5] <- "Upregulated"
CNS_inflammation$Significance[CNS_inflammation$padj < 0.05 & CNS_inflammation$log2FoldChange < -0.5] <- "Downregulated"

library(ggplot2)
library(ggrepel)
ggplot(CNS_inflammation, aes(x = log2FoldChange, y = logP, color = Significance)) +
  geom_point(size = 3, alpha = 1) +  # Make all points visible but faded
  scale_color_manual(values = c("Upregulated" = "#cc33ff", "Downregulated" = "#33ff99", "Not Significant" = "grey")) + 
  geom_text_repel(aes(label = ifelse(grepl("^ENSDARG", external_gene_name), NA, external_gene_name)), 
                  size = 4, color = "black", box.padding = 0.5, max.overlaps = 10) +
  theme_classic() +
  labs(title = "Volcano Plot of CNS Inflammation Genes",
       x = "Log2 Fold Change",
       y = "-Log10 p-value")

ggplot(CNS_inflammation, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(color = "blue", alpha = 0.7) +
  

CNS_myelination_GO %>% 
  group_by(external_gene_name) %>% 
  summarise(count = n()) %>% 
  filter(count > 1)
scale_color_manual(values = c("Upregulated" = "#cc33ff", "Downregulated" = "#33ff99", "Not Significant" = "grey"))
res_annotated %>% 
  group_by(external_gene_name) %>% 
  summarise(count = n()) %>% 
  filter(count > 1)

negative_regulation_of_dephosphorylation_(GO:0035305) 143 19 5.03688878993798e-05
genes_of_deph <- c("hephl1a", "tnfrsf1a", "ENSDARG00000099511", "fosl2", "jak1", "fosl1a", "ENSDARG00000103324", "fosb", "zgc:113314", "ENSDARG00000110878", "egfra", "b3gnt3.4", "ftr52p", "ptgis", "ENSDARG00000101992", "ENSDARG00000099384", "crema", "ifngr1l", "cflara")

CNS_deph <- res_annotated %>% 
  filter(external_gene_name %in% genes_of_deph)

CNS_deph$Significance <- "Not Significant"
CNS_deph$Significance[CNS_deph$padj < 0.05 & CNS_deph$log2FoldChange > 0.5] <- "Upregulated"
CNS_deph$Significance[CNS_deph$padj < 0.05 & CNS_deph$log2FoldChange < -0.5] <- "Downregulated"

ggplot(CNS_deph, aes(x = log2FoldChange, y = logP, color = Significance)) +
  geom_point(size = 3, alpha = 1) +  # Make all points visible but faded
  scale_color_manual(values = c("Upregulated" = "#cc33ff", "Downregulated" = "#33ff99", "Not Significant" = "grey")) + 
  geom_text_repel(aes(label = ifelse(grepl("^ENSDARG", external_gene_name), NA, external_gene_name)), 
                  size = 4, color = "black", box.padding = 0.5, max.overlaps = 10) +
  theme_classic() +
  labs(title = "Volcano Plot of negative regulation of dephosphorylation",
       x = "Log2 Fold Change",
       y = "-Log10 p-value")

# Same but for protein deph
genes_of_dephP <- c("hephl1a", "tnfrsf1a", "ENSDARG00000099511", "fosl2", "jak1", "fosl1a", "ENSDARG00000103324", "fosb", "zgc:113314", "ENSDARG00000110878", "egfra", "b3gnt3.4", "ftr52p", "ptgis", "ENSDARG00000101992", "ENSDARG00000099384", "crema", "ifngr1l", "cflara")

CNS_dephP <- res_annotated %>% 
  filter(external_gene_name %in% genes_of_dephP)

CNS_dephP$Significance <- "Not Significant"
CNS_dephP$Significance[CNS_dephP$padj < 0.05 & CNS_dephP$log2FoldChange > 0.5] <- "Upregulated"
CNS_dephP$Significance[CNS_dephP$padj < 0.05 & CNS_dephP$log2FoldChange < -0.5] <- "Downregulated"

ggplot(CNS_dephP, aes(x = log2FoldChange, y = logP, color = Significance)) +
  geom_point(size = 3, alpha = 1) +  # Make all points visible but faded
  scale_color_manual(values = c("Upregulated" = "#cc33ff", "Downregulated" = "#33ff99", "Not Significant" = "grey")) + 
  geom_text_repel(aes(label = ifelse(grepl("^ENSDARG", external_gene_name), NA, external_gene_name)), 
                  size = 4, color = "black", box.padding = 0.5, max.overlaps = 10) +
  theme_classic() +
  labs(title = "Volcano Plot of negative regulation of protein dephosphorylation",
       x = "Log2 Fold Change",
       y = "-Log10 p-value")
