# ======================================================================= #
# Gene Ontology                                                           #
# ======================================================================= #
enrichGO_results <- enrichGO(gene = unique(Nodes_Genes_Grouped$entrez_id),
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID",
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             readable      = TRUE)
# Disgenet
dgn_enrich_results <- enrichDGN(gene = unique(Nodes_Genes_Grouped$entrez_id),
                                pAdjustMethod = "fdr", # FDR -> Benjamini–Hochberg? (BH)?
                                pvalueCutoff  = 0.05)


enricher_results <- enricher(gene = unique(Nodes_Genes_Grouped$entrez_id),
)

glimpse(Nodes_Genes_Grouped)
# ======================================================================= #
# Making plots
# Make sure this script is in the same folder as your downloaded file
DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
dir.create(Folder_of_your_choice,"Folder_images")
setwd("Folder_images")


# ======================================================================= #
# Treeplot -> Gene Ontology                                               #
# ======================================================================= #
go_sim <- enrichplot::pairwise_termsim(enrichGO_results)
filename <- paste0("GO_Treeplot.png")
png(filename , width = 3000, height = 4000, res = 150)
plot(enrichplot::treeplot(go_sim, nCluster = 5, showCategory = 30))
dev.off()
# ======================================================================= #
# Treeplot -> Disgenet                                                    #
# ======================================================================= #
dgn_sim <- enrichplot::pairwise_termsim(dgn_enrich_results)
filename <- paste0("DGN_Treeplot.png")
png(filename , width = 3000, height = 4000, res = 150)
plot(enrichplot::treeplot(dgn_sim, nCluster = 5, showCategory = 20))
dev.off()
# ======================================================================= #
# Other Plots                                                             #
# ======================================================================= #
# dotplot
dotplot(enrichGO_results, showCategory = 20) + 
  ggtitle("Top Biological Processes")

dotplot(dgn_enrich_results, showCategory = 20) + 
  ggtitle("Top Biological Disease Processes")
# ======================================================================= #
# Function to count the "counts"
# ======================================================================= #
gene_count_calculation = function(enrichment_results) {
  if(is.null(enrichment_results)) return(NULL)
  
  else
    enrichment_dataframes <- as.data.frame(enrichment_results)  %>% 
      dplyr::filter(p.adjust < 0.05 ) %>%
      tidyr::separate_rows(geneID, sep ="/") %>%
      dplyr::group_by(geneID) %>%
      dplyr::summarize(gene_count = n())
  return(as.data.frame(enrichment_dataframes))
}
# ======================================================================= #
GO_count <- gene_count_calculation(enrichGO_results)
DGN_count <- gene_count_calculation(dgn_enrich_results)
