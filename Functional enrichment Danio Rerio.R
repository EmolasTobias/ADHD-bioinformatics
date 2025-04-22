if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Dr.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("tictoc")
BiocManager::install("wesanderson")
install.packages("patchwork")
BiocManager::install("clusterProfiler", force = TRUE)

#
# 0. load libraries
#
library(crayon)
library(clusterProfiler)
library(enrichplot)
library(tictoc)
library(viridis)
library(ggplot2)
library(wesanderson)
library(org.Dr.eg.db)
library(this.path)
#
# 2. read files and generate lists of genes
#

# set your working directory. This is an option, but you can change as you prefer.

setwd("C:/Users/thorb/Documents/OfvirkiFiskurinn/RNASequencinganalysis/020_vik/NEW_DEGs_DESeq2") 

# Create working files for downregulated and upregulated

filename = 'NEW_WT_vs_MUTANT_filtered.tsv'
df = read.table(filename, sep='\t', header=TRUE)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ]
list(filename)

# Testing new gene annotations

library(biomaRt)

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "drerio_gene_ensembl")

# Convert ENSEMBL IDs to ENTREZ IDs
convertedIDsBIO <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
                      filters = "ensembl_gene_id", 
                      values = row.names(df_up), 
                      mart = ensembl)

# Filter out NA values
convertedIDsBIO <- na.omit(convertedIDsBIO)

# Get ENTREZ IDs
list_one_up <- convertedIDsBIO$entrezgene_id
length(list_one_up)


# Convert ENSEMBL IDs to ENTREZ IDs
convertedIDsBIOdown <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
                      filters = "ensembl_gene_id", 
                      values = row.names(df_down), 
                      mart = ensembl)

# Filter out NA values
convertedIDsBIOdown <- na.omit(convertedIDsBIOdown)

# Get ENTREZ IDs
list_one_down <- convertedIDsBIOdown$entrezgene_id
length(list_one_down)

# END of test



ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Dr.eg.db')
list_one_up = convertedIDs$ENTREZID
length(list_one_up)

ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Dr.eg.db')
list_one_down = convertedIDs$ENTREZID
length(list_one_down)

geneLists = list('Upregulated'=list_one_up, 
                'Downregulated'=list_one_down
                )

#
# 3. run the analysis on different Ontologies
#


library(ReactomePA)  # Only if you want to use Reactome

tic()
ck = compareCluster(geneLists, 
                    fun = "enrichKEGG",  # Use KEGG for better zebrafish support
                    pvalueCutoff = 0.05, 
                    organism = "dre")  # "dre" = Danio rerio (zebrafish)
toc()


# Convert the compareCluster results to a data frame
kegg_result_z <- as.data.frame(ck)

# Create a more readable table with pathway information and gene lists
pathway_gene_table_z <- data.frame(
  Cluster = kegg_result_z$Cluster,
  Pathway_ID = kegg_result_z$ID,
  Description = kegg_result_z$Description,
  P_value = kegg_result_z$pvalue,
  Adj_P_value = kegg_result_z$p.adjust,
  Gene_Count = kegg_result_z$Count,
  Genes = kegg_result_z$geneID,
  stringsAsFactors = FALSE
)

# Display the table
View(pathway_gene_table_z)

# Save to a CSV file if needed
write.csv(pathway_gene_table_z, "kegg_pathway_genes_z.csv", row.names = FALSE)

# Load the required library for mapping gene IDs
library(org.Hs.eg.db)

# Function to convert Entrez IDs to gene symbols
entrez_to_symbol_z <- function(entrez_string) {
  entrez_ids <- strsplit(entrez_string, "/")[[1]]
  symbols <- mapIds(org.Dr.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  paste(symbols, collapse = "/")
}

# Add a column with gene symbols
pathway_gene_table_z$Gene_Symbols <- sapply(pathway_gene_table$Genes, entrez_to_symbol_z)

# View updated table
View(pathway_gene_table)

# Save the enhanced table
write.csv(pathway_gene_table, "kegg_pathway_genes_with_symbols.csv", row.names = FALSE)





p1 = dotplot(ck, size = 'count', showCategory = 5, font.size = 8) 
print(p1)

dotplot(ck, 
        showCategory = 15,
        size = "Count",
        color = "p.adjust",
        font.size = 18) +
  scale_color_continuous("blue") +
  scale_size(range = c(3, 8), name = "Gene Count") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(title = "KEGG Pathway Enrichment Analysis",
       x = "Cluster", 
       y = "Pathway")

# please use viridis and a more intuitive direction for P values, where lower is better 
p3 = p1 +  scale_fill_viridis(direction=-1)
print(p3)

# but log scale is more appropriate
my_log_breaks = round(log10(0.05)):round(log10(min(ck@compareClusterResult$p.adjust)))
my_breaks = 10**my_log_breaks
p4 = p1 +  scale_fill_viridis(direction=-1, trans="log", breaks = my_breaks)
print(p4)

# and I have a preference for cividis, but this is just personal preference
p5 = p1 +  scale_fill_viridis(direction=-1, trans="log", breaks = my_breaks, option='cividis')
print(p5)
# arguably this plot communicates best data patterns, IMHO
library(org.Dr.eg.db)
ck2 = compareCluster(geneLists, 
                    fun = "enrichGO",
                    pvalueCutoff = 0.05,
                    OrgDb = org.Dr.eg.db, 
                    ont = "ALL")  # Choose "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)

GO1 = dotplot(ck2, size = 'count', showCategory = 5, font.size = 8) 
print(GO1)

# importantly, store your fuctional enrichment in a form of table which will be a supplementary file of your paper

storage_file = 'clusterProfiler_enrichments_RP.tsv'
write.table(ck@compareClusterResult, storage_file, quote=FALSE, sep='\t')

write.table(ck2@compareClusterResult, 'funnyname.tsv', quote=FALSE, sep='\t')

write.table(ck@compareClusterResult, 'kegg_results.tsv', quote=FALSE, sep='\t')

# Extract genes from KEGG results
kegg_genes = ck@compareClusterResult %>% dplyr::select(Description, geneID)

# Extract genes from GO results
go_genes = ck2@compareClusterResult %>% dplyr::select(Description, geneID)


intersect(go_genes, kegg_genes)  # Find common genes



for (x in convertedIDs) {
  cat(paste(x, "\n"))
}

#

ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb='org.Dr.eg.db')
list_one_up = convertedIDs$GENEID
length(list_one_up_gene)

ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb='org.Dr.eg.db')
list_one_down_gene = convertedIDs$GENEID
length(list_one_down_gene)

for (x in convertedIDs) {
  cat(paste(x, "\n"))
}

# Call fort ck to look at the 12 genes that correlate to Phagosome pathway

kegg_tsv = 'kegg_results.tsv'
kegg_phagosome = read.table(kegg_tsv, sep='\t', header=TRUE)

converted_kegg_phagosome = bitr(kegg_phagosome$geneID, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb='org.Dr.eg.db')
list(kegg_phagosome)

# Load necessary library
library(clusterProfiler)
library(org.Dr.eg.db)  # Ensure you have this installed for zebrafish

# Split the geneID column into individual ENTREZ IDs
phagosome_gene_list <- strsplit(kegg_phagosome$geneID, "/")[[1]]

# Convert ENTREZ IDs to SYMBOLs
converted_genes <- bitr(phagosome_gene_list, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Dr.eg.db)

# Create a mapping of ENTREZID to SYMBOL
gene_map <- setNames(converted_genes$SYMBOL, converted_genes$ENTREZID)

# Replace ENTREZ IDs with their corresponding SYMBOLs in the dataset
kegg_phagosome$geneID <- sapply(strsplit(kegg_phagosome$geneID, "/"), function(ids) {
  paste(gene_map[ids], collapse = "/")
})

# View the updated dataset
print(kegg_phagosome)


# CLAUDE

# Convert the GO compareCluster results to a data frame
go_result_ck2 <- as.data.frame(ck2)

# Create a table with GO term information and gene lists
go_gene_table_ck2 <- data.frame(
  Cluster = go_result_ck2$Cluster,
  GO_ID = go_result_ck2$ID,
  Description = go_result_ck2$Description,
  P_value = go_result_ck2$pvalue,
  Adj_P_value = go_result_ck2$p.adjust,
  Gene_Count = go_result_ck2$Count,
  Genes = go_result_ck2$geneID,
  stringsAsFactors = FALSE
)

# Display the table
View(go_gene_table_ck2)

# Save to a CSV file if needed
write.csv(go_gene_table_ck2, "go_term_genes_ck2.csv", row.names = FALSE)

# Already have org.Hs.eg.db loaded from your GO analysis

# Function to convert Entrez IDs to gene symbols (if not defined earlier)
entrez_to_symbol <- function(entrez_string) {
  entrez_ids <- strsplit(entrez_string, "/")[[1]]
  symbols <- mapIds(org.Dr.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  paste(symbols, collapse = "/")
}

# Add a column with gene symbols
go_gene_table_ck2$Gene_Symbols <- sapply(go_gene_table$Genes, entrez_to_symbol)

# View updated table
View(go_gene_table)

# Save the enhanced table
write.csv(go_gene_table, "go_term_genes_with_symbols.csv", row.names = FALSE)

# Create a more customized GO dotplot
dotplot(ck2, 
        showCategory = 5,
        size = "Count",
        color = "p.adjust",
        font.size = 12) +
  scale_color_continuous(low = "red", high = "blue", name = "Adjusted p-value") +
  scale_size(range = c(3, 8), name = "Gene Count") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(title = "GO Biological Process Enrichment Analysis",
       x = "Cluster", 
       y = "GO Term")
