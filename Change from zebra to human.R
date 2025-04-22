# Install Bioconductor if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

# Load library
library(biomaRt)
script_path <- "C:/Users/thorb/Documents/OfvirkiFiskurinn/RNASequencinganalysis/020_vik/NEW_DEGs_DESeq2"
script_path
setwd(script_path) 
filename = 'NEW_WT_vs_MUTANT_annotated_update.tsv'
res_annotated = read.table(filename, sep='\t', header=TRUE)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ]
# Connect to Ensembl
mart_drerio <- useMart("ensembl", dataset = "drerio_gene_ensembl")

# Define gene list
zebrafish_genes <- c("hephl1a", "tnfrsf1a", "ENSDARG00000099511", "fosl2", 
                     "jak1", "fosl1a", "ENSDARG00000103324", "fosb", "zgc:113314",
                     "ENSDARG00000110878", "egfra", "b3gnt3.4", "ftr52p", "ptgis", 
                     "ENSDARG00000101992", "ENSDARG00000099384", "crema", "ifngr1l", "cflara")

# Get human orthologs
orthologs <- getBM(attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
                   filters = "external_gene_name", 
                   values = res_annotated$external_gene_name, 
                   mart = mart_drerio)

# Print results
print(orthologs)

res_annotated_ortho <- merge(res_annotated, orthologs, by.x = "external_gene_name", 
                       by.y = "external_gene_name", all.x = TRUE)

# Rename the new column for clarity
colnames(res_annotated_ortho)[ncol(res_annotated_ortho)] <- "human_ortholog"

# Print the updated dataframe
head(res_annotated_ortho)

# Total number of genes in res_annotated
total_genes <- nrow(res_annotated_ortho)

# Count missing orthologs (both NA and empty strings)
missing_orthologs <- sum(is.na(res_annotated_ortho$human_ortholog) | res_annotated_ortho$human_ortholog == "")

# Count genes that have a human ortholog
found_orthologs <- total_genes - missing_orthologs

# Print results
cat("Total genes:", total_genes, "\n")
cat("Genes with human orthologs:", found_orthologs, "\n")
cat("Genes without human orthologs:", missing_orthologs, "\n")

# Count occurrences of each human ortholog
ortholog_counts <- table(res_annotated_ortho$human_ortholog)

# Convert to a data frame for easier viewing
ortholog_df <- as.data.frame(ortholog_counts)

# Rename columns
colnames(ortholog_df) <- c("human_ortholog", "count")

# Sort by count in descending order
ortholog_df <- ortholog_df[order(-ortholog_df$count), ]

# Print the first few rows
head(ortholog_df)
View(ortholog_df)

# Create the filtered dataset
res_annotated_ortho_filtered <- res_annotated_ortho[
  res_annotated_ortho$padj < 0.05 & abs(res_annotated_ortho$log2FoldChange) > 0.5, 
]

# Print the number of rows in the filtered dataset
cat("Number of significant genes:", nrow(res_annotated_ortho_filtered), "\n")

# Show the first few rows
head(res_annotated_ortho_filtered)

# Remove rows where human_ortholog is NA or an empty string
res_annotated_ortho_filtered <- res_annotated_ortho_filtered[
  !is.na(res_annotated_ortho_filtered$human_ortholog) & 
    res_annotated_ortho_filtered$human_ortholog != "", 
]

# Remove duplicates, keeping only the first occurrence of each human_ortholog
res_annotated_ortho_filtered <- res_annotated_ortho_filtered[
  !duplicated(res_annotated_ortho_filtered$human_ortholog), 
]

# Remove duplicates of zebrafish genes_ortholog
res_annotated_ortho_filtered <- res_annotated_ortho_filtered[
  !duplicated(res_annotated_ortho_filtered$external_gene_name), 
]

# Print the number of remaining genes
cat("Number of unique human orthologs:", nrow(res_annotated_ortho_filtered), "\n")

# Show the first few rows
head(res_annotated_ortho_filtered)

# Install required packages if not installed
if (!requireNamespace("clusterProfiler", quietly = TRUE)) install.packages("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("enrichplot", quietly = TRUE)) install.packages("enrichplot")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("DOSE", quietly = TRUE)) BiocManager::install("DOSE")

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # Human gene annotation database
library(enrichplot)
library(ggplot2)
library(DOSE)

# Convert human gene symbols to Entrez IDs
gene_list <- res_annotated_ortho_filtered$human_ortholog

# Map to Entrez IDs using org.Hs.eg.db
entrez_ids <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# Merge with original dataset to ensure correct mapping
res_annotated_ortho_filtered <- merge(res_annotated_ortho_filtered, entrez_ids, by.x="human_ortholog", by.y="SYMBOL", all.x=TRUE)

# Remove rows where Entrez ID is NA (genes that couldn't be mapped)
res_annotated_ortho_filtered <- na.omit(res_annotated_ortho_filtered)

go_results <- enrichGO(
  gene          = res_annotated_ortho_filtered$ENTREZID, 
  OrgDb         = org.Hs.eg.db, 
  keyType       = "ENTREZID", 
  ont           = "BP",  # Biological Process
  pAdjustMethod = "BH", 
  pvalueCutoff  = 0.05, 
  qvalueCutoff  = 0.05
)

# View top GO terms
head(go_results)

# Bar plot of top enriched GO terms
barplot(go_results, showCategory=15, title="GO Enrichment Analysis")

# Dot plot (alternative visualization)
dotplot(go_results, showCategory=15, title="GO Enrichment Analysis")

if (!requireNamespace("ReactomePA", quietly = TRUE)) BiocManager::install("ReactomePA")
library(ReactomePA)

reactome_results <- enrichPathway(
  gene         = res_annotated_ortho_filtered$ENTREZID,
  organism     = "human",
  pvalueCutoff = 0.05
)

# Plot Reactome enrichment results
dotplot(reactome_results, showCategory=15, title="Reactome Pathway Enrichment")

kegg_results <- enrichKEGG(
  gene         = res_annotated_ortho_filtered$ENTREZID,
  organism     = "hsa",  # Homo sapiens
  pvalueCutoff = 0.05
)

# View KEGG results
head(kegg_results)

# Plot KEGG enrichment results
dotplot(kegg_results, showCategory=15, title="KEGG Pathway Enrichment")

res_annotated_ortho_filtered$Significance <- "Not Significant"
res_annotated_ortho_filtered$Significance[CNS_dephP$padj < 0.05 & res_annotated_ortho_filtered$log2FoldChange > 0.5] <- "Upregulated"
res_annotated_ortho_filtered$Significance[CNS_dephP$padj < 0.05 & res_annotated_ortho_filtered$log2FoldChange < -0.5] <- "Downregulated"

#compare clusters
# Define thresholds
logFC_cutoff <- 1  # Change this if you want a stricter threshold

# Upregulated genes
upregulated_genes <- res_annotated_ortho_filtered[
  res_annotated_ortho_filtered$log2FoldChange > logFC_cutoff, 
]

# Downregulated genes
downregulated_genes <- res_annotated_ortho_filtered[
  res_annotated_ortho_filtered$log2FoldChange < -logFC_cutoff, 
]

# Print gene counts
cat("Upregulated genes:", nrow(upregulated_genes), "\n")
cat("Downregulated genes:", nrow(downregulated_genes), "\n")

# Convert upregulated genes to Entrez IDs
upregulated_entrez <- bitr(
  upregulated_genes$human_ortholog, 
  fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db
)

# Convert downregulated genes to Entrez IDs
downregulated_entrez <- bitr(
  downregulated_genes$human_ortholog, 
  fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db
)

library(clusterProfiler)
library(org.Hs.eg.db)

# Define threshold for log2FC
logFC_cutoff <- 1  # Change if needed

# Create a named list of gene groups
gene_list <- list(
  Upregulated = res_annotated_ortho_filtered$human_ortholog[
    res_annotated_ortho_filtered$log2FoldChange > logFC_cutoff
  ],
  Downregulated = res_annotated_ortho_filtered$human_ortholog[
    res_annotated_ortho_filtered$log2FoldChange < -logFC_cutoff
  ]
)

# Convert gene symbols to Entrez IDs
gene_list_entrez <- lapply(gene_list, function(genes) {
  bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
})

# Run compareCluster for GO enrichment
go_clusters <- compareCluster(
  geneCluster = gene_list_entrez,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Plot the results
dotplot(go_clusters, showCategory=15, title="GO Enrichment: Upregulated vs. Downregulated Genes")

kegg_clusters <- compareCluster(
  geneCluster = gene_list_entrez,
  fun = "enrichKEGG",
  organism = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

dotplot(kegg_clusters, showCategory=15, title="KEGG Pathway Enrichment: Upregulated vs. Downregulated Genes") +
  scale_size_continuous(name="Gene Count") + geom_text_repel(aes(label = Count), size = 6)
dotplot(go_clusters, showCategory=10, title="GO Enrichment") +
  scale_size_continuous(name="Gene Count") + geom_text_repel(aes(label = Count), size = 4, nudge_x = 0.1, direction = "y") # Changes legend title to "Gene Count"

go_clusters@compareClusterResult$GeneRatio <- go_clusters@compareClusterResult$Count

dotplot(go_clusters, showCategory=15, title="GO Enrichment")

kegg_clusters@compareClusterResult$GeneRatio <- kegg_clusters@compareClusterResult$Count

dotplot(kegg_clusters, showCategory = 15, title = "Kegg Enrichment") + 
  geom_point(aes(size = Count)) +  # Use 'Count' for dot size
  geom_text_repel(aes(label = Count), size = 6) +  # Add labels with gene count
  scale_size_continuous(name = "Gene Count") 

### Kóði frá CLAUDE hér fyrir neðan, bingchilling

# Convert the compareCluster results to a data frame
kegg_result_df <- as.data.frame(kegg_clusters)

# Create a more readable table with pathway information and gene lists
pathway_gene_table <- data.frame(
  Cluster = kegg_result_df$Cluster,
  Pathway_ID = kegg_result_df$ID,
  Description = kegg_result_df$Description,
  P_value = kegg_result_df$pvalue,
  Adj_P_value = kegg_result_df$p.adjust,
  Gene_Count = kegg_result_df$Count,
  Genes = kegg_result_df$geneID,
  stringsAsFactors = FALSE
)

# Display the table
View(pathway_gene_table)

# Save to a CSV file if needed
write.csv(pathway_gene_table, "kegg_pathway_genes.csv", row.names = FALSE)

# Load the required library for mapping gene IDs
library(org.Hs.eg.db)

# Function to convert Entrez IDs to gene symbols
entrez_to_symbol <- function(entrez_string) {
  entrez_ids <- strsplit(entrez_string, "/")[[1]]
  symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  paste(symbols, collapse = "/")
}

# Add a column with gene symbols
pathway_gene_table$Gene_Symbols <- sapply(pathway_gene_table$Genes, entrez_to_symbol)

# View updated table
View(pathway_gene_table)

# Save the enhanced table
write.csv(pathway_gene_table, "kegg_pathway_genes_with_symbols.csv", row.names = FALSE)

# Create a more customized dotplot
dotplot(kegg_clusters, 
        showCategory = 15,
        size = "Count",
        color = "p.adjust",
        font.size = 18) +
  scale_color_continuous(low = "red", high = "blue", name = "Adjusted p-value") +
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

png("kegg_dotplot.png", width = 1000, height = 800, res = 120)

dotplot(kegg_clusters, 
        showCategory = 15,
        size = "Count",
        color = "p.adjust")
dev.off()

###### GO NOW

# Convert the GO compareCluster results to a data frame
go_result_df <- as.data.frame(go_clusters)

# Create a table with GO term information and gene lists
go_gene_table <- data.frame(
  Cluster = go_result_df$Cluster,
  GO_ID = go_result_df$ID,
  Description = go_result_df$Description,
  P_value = go_result_df$pvalue,
  Adj_P_value = go_result_df$p.adjust,
  Gene_Count = go_result_df$Count,
  Genes = go_result_df$geneID,
  stringsAsFactors = FALSE
)

# Display the table
View(go_gene_table)

# Save to a CSV file if needed
write.csv(go_gene_table, "go_term_genes.csv", row.names = FALSE)

# Already have org.Hs.eg.db loaded from your GO analysis

# Function to convert Entrez IDs to gene symbols (if not defined earlier)
entrez_to_symbol <- function(entrez_string) {
  entrez_ids <- strsplit(entrez_string, "/")[[1]]
  symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  paste(symbols, collapse = "/")
}

# Add a column with gene symbols
go_gene_table$Gene_Symbols <- sapply(go_gene_table$Genes, entrez_to_symbol)

# View updated table
View(go_gene_table)

# Save the enhanced table
write.csv(go_gene_table, "go_term_genes_with_symbols.csv", row.names = FALSE)

# Create a more customized GO dotplot
dotplot(go_clusters, 
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
png("go_dotplot.png", width = 1000, height = 800, res = 120)
dotplot(go_clusters, 
        showCategory = 15,
        size = "Count",
        color = "p.adjust")
dev.off()