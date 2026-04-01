# ============================================ #
#                   Step 0:                    #
#         Installing & Loading packages        #
# ============================================ #
Packages_needed <- c("tidyverse",
                     "msigdbr",
                     "biomaRt",
                     "STRINGdb",
                     "clusterProfiler",
                     "org.Hs.eg.db",
                     "DOSE",
                     "enrichplot")

install_and_load <- function(package_names) {
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  Available_on_BiocManager <- BiocManager::available()
  
  for (package_name in package_names) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      if (package_name %in% Available_on_BiocManager) {
        BiocManager::install(package_name, ask = FALSE, update = FALSE)
      } else {
        install.packages(package_name)
      }
    }
  }
  invisible(lapply(package_names, library, character.only = TRUE))
  message("Packages are installed & Loaded!")
}

install_and_load(Packages_needed)

# ======================================================================= #
# STEP 3: Enrichment Analysis                                             #
# ======================================================================= #
# Part 1:  Over-representation analysis                                   #
#           Gene and Disease ontology                                     #
# Part 2: Counting the genes that are enriched                            #
#           GO over-representation analysis                               #
# ======================================================================= #
#                         Last Part is Optional                           #
# ======================================================================= #
# Part 3: Pairwise Analysis and making plots                              #
# ======================================================================= #
# ======================================================================= #
# Part 1: Gene Ontology                                                   #
#   GO over-representation analysis                                       #
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
# ======================================================================= #
# Disease Ontology                                                        #
# ======================================================================= #
dgn_enrich_results <- enrichDGN(gene = unique(Nodes_Genes_Grouped$entrez_id),
                                pAdjustMethod = "BH", 
                                pvalueCutoff  = 0.05)


# ======================================================================= #
# Part 2:                                                                 #
# ======================================================================= #
# Function to count the "counts"                                          #
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
# ======================================================================= #



# ======================================================================= #
# Part 3:                                                                 #
# ======================================================================= #
# If making plots or gene enrichment regarding the entire network is      #
# needed                                                                  #
# ======================================================================= #

Folder_of_your_choice <- rstudioapi::selectDirectory()
Folder_images_export_data <- file.path(Folder_of_your_choice, "Folder_export_data")
dir.create(Folder_export_data, recursive = TRUE)
setwd(Folder_export_data)
# ======================================================================= #
# Making plots
# Make sure this script is in the same folder as your downloaded file
Folder_of_your_choice <- rstudioapi::selectDirectory()
Folder_images_GeneOntology <- file.path(Folder_of_your_choice, "Folder_images_GeneOntology")
dir.create(Folder_images_GeneOntology, recursive = TRUE)
setwd(Folder_images_GeneOntology)
# ======================================================================= #
# Treeplot -> Gene Ontology                                               #
# ======================================================================= #
go_sim <- enrichplot::pairwise_termsim(enrichGO_results)
GO_Treeplot.png <- paste0("GO_Treeplot.png")
png(GO_Treeplot.png , width = 3000, height = 4000, res = 150)
plot(enrichplot::treeplot(go_sim, nCluster = 5, showCategory = 30))
dev.off()
# ======================================================================= #
# Treeplot -> Disgenet                                                    #
# ======================================================================= #
dgn_sim <- enrichplot::pairwise_termsim(dgn_enrich_results)
DGN_Treeplot.png <- paste0("DGN_Treeplot.png")
png(DGN_Treeplot.png , width = 3000, height = 4000, res = 150)
plot(enrichplot::treeplot(dgn_sim, nCluster = 5, showCategory = 20))
dev.off()
# ======================================================================= #
# Other Plots                                                             #
# ======================================================================= #
# dotplot
# ======================================================================= #
png("GO_dotplot.png", width = 3000, height = 4000, res = 150)
dotplot_enrichGO_results <- dotplot(enrichGO_results, showCategory = 20) +
  ggtitle("Top Biological Processes")
print(dotplot)   
dev.off()

# Create plot
p <- dotplot(enrichGO_results, showCategory = 15) +
  ggtitle("Top Enriched Disease-Related Processes") +
  scale_color_gradient(low = "#4575B4", high = "#D73027", name = "Adjusted p-value") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
ggsave("GO_dotplot_publication.png", plot = p,
       width = 8, height = 10, dpi = 300)
# ======================================================================= #
png("DNG_dotplot.png", width = 3000, height = 4000, res = 150)
dotplot_dgn_results <- dotplot(dgn_enrich_results, showCategory = 20) + 
  ggtitle("Top Biological Disease Processes")
print(dotplot_dgn_results)
dev.off()



# ======================================================================= #
enrichGO_results_DF <- as.data.frame(enrichGO_results) %>%
  tidyr::separate_rows(geneID, sep ="/") %>%
  dplyr::select(geneID, Description, p.adjust, qvalue)%>%
  dplyr::group_by(Description, p.adjust, qvalue) %>%
  dplyr::summarize()

dgn_enrich_results_DF <- as.data.frame(dgn_enrich_results) %>%
  tidyr::separate_rows(geneID, sep ="/") %>%
  dplyr::select(geneID, Description, p.adjust, qvalue)%>%
  dplyr::group_by(Description, p.adjust, qvalue) %>%
  dplyr::summarize()

# ======================================================================= #
write.table(enrichGO_results_DF, "enrichGO_results.csv", sep = ";", row.names = FALSE)
write.table(dgn_enrich_results_DF, "dgn_enrich_results.csv", sep = ";", row.names = FALSE)