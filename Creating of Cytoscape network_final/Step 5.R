# ============================================ #
#                   Step 0:                    #
#         Installing & Loading packages        #
# ============================================ #
Packages_needed <- c("tidyverse", 
                     "RCy3", 
                     "rstudioapi", 
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
# STEP 5: Post Cytoscape Network 2 Analysis                               #
# ======================================================================= #
# ======================================================================= #
# STEP 4: Post Cytoscape Network-Analysis                                 #
#   Part 0: Making a folder for exporting data and images                 #
#   Part 1: Extracting Data from the network                              #
#   Part 2: Enrichment for clusters                                       #
#       I. Gene Ontology & exporting as .csv-file                         #
#       II. Disease Ontology & exporting as .csv-file                     #
#   Part 3: Enriching and plotting per cluster                            #
#       I. Separating Clusters to two lists                               #
#           (Cluster name and NCBI/Entrez-ID)                             #  
#           EnrichGO & EnrichDGN                                          #
#     (I)I. optional Perform pairwise_termsim on clusters)                #
#     (I)II. Plotting & exporting the plots                               #
#             and exporting data file (optional)                          #
# ======================================================================= #
# Part 0: Making a folder for exporting data and images                   #
# ======================================================================= #
# For the raw-data/.csv-files                                             #
# ======================================================================= #
Folder_of_your_choice <- rstudioapi::selectDirectory()
Folder_export_data <- file.path(Folder_of_your_choice, "Folder_export_data")
dir.create(Folder_export_data, recursive = TRUE)
# ======================================================================= #
# For the images/plots                                                    #
# ======================================================================= #
Folder_of_your_choice <- rstudioapi::selectDirectory()
Folder_images_clusters <- file.path(Folder_of_your_choice, "Folder_images_clusters")
dir.create(Folder_images_clusters, recursive = TRUE)

# ======================================================================= #
# Part 1: Extracting Data from the network                                #
#                                                                         #
# ======================================================================= #
# Before you extract the columns from the network; make sure you          #
# make column for AutoAnnotate                                            #
# ======================================================================= #
# CREATE COLUMN -> AUTOANNOTARE -> 3 Stripes -> CREATE COLUMN OF          #
#                                               CLUSTERLABLES             #
# ======================================================================= #



node_table_clustered <- RCy3::getTableColumns(table = "node")



node_table_with_ClusterNames <- node_table_clustered %>%
  dplyr::mutate(cluster_name = AutoAnnotate_GLay_Clusters, cluster_id = `__glayCluster`) %>%
  dplyr::mutate(cluster_name = stringr::str_to_title(cluster_name)) %>% # Capitalize
  dplyr::mutate(cluster_name = stringr::str_remove_all(cluster_name,"Wp|Kegg|Reactome")) %>%
  dplyr::mutate(cluster_name = stringr::str_replace_all(cluster_name, "_", " ")) %>%
  dplyr::mutate(cluster_name = stringr::str_trim(cluster_name)) %>%  
  dplyr::select(-`__glayCluster`, -AutoAnnotate_GLay_Clusters) 


# ======================================================================= #
# Dataframe with Cluster_ID and ClusterNames exporting it                 #
# ======================================================================= #
Cluster_ID_ClusterNames <- node_table_with_ClusterNames %>%
  dplyr::select(cluster_id, cluster_name)
# ======================================================================= #
# Setting folder & exporting data                                         #
# ======================================================================= #
setwd(Folder_export_data)
# ======================================================================= #
write.table(Cluster_ID_ClusterNames, "Cluster_ID_ClusterNames.csv", sep = ";", row.names = FALSE)
# ======================================================================= #
# cluster_pathway_summary and exporting them                              #
# ======================================================================= #
cluster_pathway_summary <- node_table_clustered %>%
  dplyr::filter(type == "Pathway") %>%
  dplyr::select(cluster = `__glayCluster`, label, interaction) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarize(
    pathway_count = n(),
    # Create a comma-separated list of the pathways in this cluster
    pathway_list = paste(unique(label), collapse = "; "),
    # Identify which databases are contributing to this cluster
    sources = paste(unique(interaction), collapse = ", "),
  ) %>%
  dplyr::arrange(desc(pathway_count))

cluster_pathway_summary_with_names <- cluster_pathway_summary %>%
  dplyr::left_join(names_for_clusters, by = c("cluster" = "cluster_id"))  %>%
  dplyr::mutate(sources = stringr::str_remove_all(sources, "CP:")) %>%
  dplyr::mutate(sources = stringr::str_replace_all(sources, "_", " ")) %>%
  dplyr::mutate(pathway_list = stringr::str_remove_all(pathway_list, "KEGG_MEDICUS_|WP_|REACTOME_")) %>%
  dplyr::mutate(pathway_list = stringr::str_replace_all(pathway_list, "_", " "))

# ======================================================================= #
# Setting folder & exporting data                                         #
# ======================================================================= #
setwd(Folder_images_export_data)
write.table(cluster_pathway_summary_with_names, 
            "cluster_pathway_summary_with_names.csv", 
            sep = ";", row.names = FALSE)

# ======================================================================= #
# Step 2: Enrichment for clusters                                         #
#     Part 1:                                                             #
#        Gene Ontology & exporting as .csv-file                           #
#     Part 2:                                                             #
#         Disease Ontology & exporting as .csv-file                       #
# ======================================================================= #
# Part 1:                                                                 #
#   Gene Ontology                                                         #
# ======================================================================= #
cluster_go_results <- node_table_with_ClusterNames %>%
  dplyr::filter(type != "Pathway") %>% # Only genes
  dplyr::group_by(cluster_name) %>%
  dplyr::do(as.data.frame(
    clusterProfiler::enrichGO(
      gene = .$entrez_id, 
      OrgDb = org.Hs.eg.db, 
      ont = "BP", 
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    )
  ))

# Selecting top 2 terms for each cluster
cluster_go_results_df <- cluster_go_results %>%
  dplyr::group_by(cluster_name) %>%
  dplyr::slice_head(n = 2) %>%
  dplyr::select(cluster_name , Description, p.adjust, qvalue)
# ======================================================================= #
# Setting folder & exporting data                                         #
# ======================================================================= #
setwd(Folder_export_data)
write.table(cluster_go_results_df, "cluster_GO_results.csv", sep = ";", row.names = FALSE)


# ======================================================================= #
# Part 2:                                                                 #
#   Disease Ontology                                                      #
# ======================================================================= #
cluster_dgn_results <- node_table_with_ClusterNames %>%
  dplyr::filter(type != "Pathway") %>% # Only genes
  dplyr::group_by(cluster_name) %>%
  dplyr::do(as.data.frame(
    DOSE::enrichDGN(.$entrez_id,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    readable      = TRUE)
  ))
# Selecting top 2 terms for each cluster
cluster_dgn_results_df <- cluster_dgn_results %>%
  left_join(names_for_clusters, by = "cluster_id") %>%
  dplyr::group_by(cluster_name) %>%
  dplyr::arrange(desc(p.adjust))  %>%
  dplyr::slice_head(n = 2) %>%
  dplyr::select(cluster_id, cluster_name , Description, p.adjust, qvalue)
# ======================================================================= #
# Setting folder & exporting data                                         #
# ======================================================================= #
setwd(Folder_export_data)
write.table(cluster_dgn_results_df, "cluster_DGN_results.csv", sep = ";", row.names = FALSE)

# ======================================================================= #
#                                                                         #
# ======================================================================= #
# Step 2: Enrichment for clusters                                         #
#     Part 1:                                                             #
#        Gene Ontology & exporting as .csv-file                           #
#     Part 2:                                                             #
#         Disease Ontology & exporting as .csv-file                       #
# ======================================================================= #
# Part 1:                                                                 #
#   Gene Ontology                                                         #
# ======================================================================= #
# Step 3: Enriching and plotting per cluster                              #
#     Part 1:                                                             #
#           Separating Clusters                                           #
#           to two lists (Cluster name and NCBI/Entrez-ID)                #  
#           EnrichGO & EnrichDGN                                          #
#     Part 2: (optional)                                                  #
#           Perform pairwise_termsim on clusters                          #
#     Part 3:                                                             #
#         Plotting & exporting the plots                                  #
#     Part 4:                                                             #
#         Exporting data file (optional)                                  #
# ======================================================================= #
# ======================================================================= #
#  Step 0:                                                                #
# ======================================================================= #

# ======================================================================= #
setwd(Folder_images_clusters)
# ======================================================================= #
# Step 1:                                                                 #
#   Separating Clusters to two lists (Cluster name and NCBI/Entrez-ID)    #
# ======================================================================= #
gene_list_by_cluster <- node_table_with_ClusterNames %>% 
  dplyr::filter(type != "Pathway", !is.na(entrez_id)) %>% 
  dplyr::select(entrez_id, cluster_name) %>% # Splits into a list where each element is a cluster's genes { split(.$entrez_id, .$cluster_name) }
  { split(.$entrez_id, .$cluster_name) }


# If there are NAs
gene_list_by_cluster <- lapply(gene_list_by_cluster, function(x) x[!is.na(x)])
# ======================================================================= #
# Preforming EnrichGO on the clusters                                     #
# ======================================================================= #
Cluster_comparsion_go <- clusterProfiler::compareCluster(geneCluster = gene_list_by_cluster, 
                                                         fun = "enrichGO", 
                                                         OrgDb = org.Hs.eg.db,
                                                         ont="BP",
                                                         pvalueCutoff = 0.05)
# ======================================================================= #
# Preforming EnrichDGN on the clusters                                    #
# ======================================================================= #
Cluster_comparsion_dgn <- clusterProfiler::compareCluster(
  geneCluster = gene_list_by_cluster, 
  fun = "enrichDGN", 
  pvalueCutoff = 0.05
)
# ======================================================================= #
# Step 2: (optional?)                                                     #
#   Perform pairwise_termsim on clusters                                  #
# ======================================================================= #
# EnrichGO TERMSIM
Cluster_comparsion_go_sim <- enrichplot::pairwise_termsim(Cluster_comparsion_go)
Cluster_comparsion_go_sim@compareClusterResult$Description <- str_wrap(Cluster_comparsion_go_sim@compareClusterResult$Description, width = 45)


# EnrichDGN TERMSIM
Cluster_comparsion_dgn_sim <- enrichplot::pairwise_termsim(Cluster_comparsion_dgn)
Cluster_comparsion_dgn_sim@result$Description <- str_wrap(Cluster_comparsion_dgn_sim@result$Description, width = 40)
# ======================================================================= #
# Step 3:                                                                 #
#   Making the Dot Plot                                                   #
# ======================================================================= #
# EnrichGO                                                                #
# ======================================================================= #
# Setting folder                                                          #
# ======================================================================= #
p_dotplot_GO <- dotplot(Cluster_comparsion_go, showCategory = 2) +
  ggtitle("Top Biological Process per Cluster") +
  scale_color_gradient(low = "#4575B4", high = "#D73027") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Rotate cluster names
    axis.text.y = element_text(size = 9),                        
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) 

# Saving dotplot
ggsave("Cluster_GO.png", plot = p_dotplot_GO, width = 14, height = 20, dpi = 300)
# ======================================================================= #
# EnrichDGN                                                               #
# ======================================================================= #
p_dotplot_dng <- dotplot(Cluster_comparsion_dgn, showCategory = 2) +
  ggtitle("Top Disease Associations per Cluster") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Rotate cluster names
    axis.text.y = element_text(size = 9),                        # Make disease names readable
    panel.grid.major = element_line(color = "grey90")
  )

# Saving dotplot
ggsave("Cluster_DGN.png", plot = p_dotplot_dng, width = 12, height = 10, dpi = 300)
