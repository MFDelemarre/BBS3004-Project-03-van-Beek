# ============================================ #
#                   Step 0:                    #
#         Installing & Loading packages        #
# ============================================ #
Packages_needed <- c("tidyverse", 
                     "RCy3", 
                     "rstudioapi", 
                     "igraph")

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
# STEP 4: Preparation for And launching                                   #
#         Cytoscape:                                                      #
#   1. Finishing the Edges                                                #
#   2. Finishing the Nodes                                                #
#   3. Creating the network                                               #
# ======================================================================= #
# Part 1: Edges                                                           #
#   Pathway-Edges                                                         #
# ======================================================================= #
pathway_edges_clean <- pathway_edges_preference %>%
  dplyr::select(source = source,
                target = target,       # This is the Pathway Name
                interaction)  %>%      # Pathway name acts as its own symbol
  dplyr::filter(!is.na(source) & !is.na(target))
# ======================================================================= #
# Combining the two edges                                                 #
# ======================================================================= #
all_edges <- bind_rows(Edges_PPI_final, pathway_edges_clean) %>%
  distinct(source, target, .keep_all = TRUE)
# ======================================================================= #
# STEP 4: Launching Cytoscape:                                            #
#   2. Finishing the Nodes                                                #
# ======================================================================= #
# Part 2: Nodes
#   Pathway-Nodes
# ======================================================================= #
node_pathways <- pathway_edges_clean %>%
  dplyr::select(id = target, interaction) %>%
  unique() %>%
  dplyr::mutate(label = id, type = "Pathway", log_bio_impact = 0.2, visual_group = interaction)
# ======================================================================= #
# Part 2: Nodes
#   Gene-Nodes
# ======================================================================= #
node_genes <- Nodes_Genes_Grouped %>%
  dplyr::left_join(GO_count, by = c("label" = "geneID"))  %>%
  dplyr::mutate(GO_count = tidyr::replace_na(gene_count, 0))  %>% #?na_if???
  dplyr::select(-gene_count) %>%
  dplyr::left_join(DGN_count, by = c("entrez_id" = "geneID"))  %>%
  dplyr::mutate(DGN_count = tidyr::replace_na(gene_count, 0)) %>%
  dplyr::select(-gene_count) %>%
  dplyr::mutate(bio_impact = GO_count + DGN_count) %>%
  dplyr::mutate(log_bio_impact = log1p(bio_impact), interaction = "Gene", visual_group = type) %>%
  as.data.frame()
# ======================================================================= #
# Combining the two nodes
# ======================================================================= #
all_nodes <- bind_rows(node_genes, node_pathways) %>% 
  dplyr::distinct(id, .keep_all = TRUE) 
# ======================================================================= #
# STEP 4: Launching Cytoscape: 
# ======================================================================= #
RCy3::cytoscapePing()

RCy3::createNetworkFromDataFrames(nodes = all_nodes, 
                                  edges = all_edges, 
                                  title = "FUMA_Enriched_Network_test_fuck")


style_name <- "Cytoscape_Network_GWAS_T2DM"

if (!(style_name %in% RCy3::getVisualStyleNames())) {
  RCy3::createVisualStyle(style_name)
}

RCy3::setEdgeLineWidthDefault(5, style.name = style_name)



RCy3::setNodeSizeMapping('log_bio_impact', 
                         table.column.values = range(all_nodes$log_bio_impact), 
                         sizes = c(30, 100), 
                         mapping.type = "c", # continous
                         style.name = style_name) # 'c' for continuous mapping


# c("PPI", "CP:KEGG_MEDICUS", "CP:REACTOME", "CP:WIKIPATHWAYS"),
# c("#999999", "#E69F00", "#56B4E9", "#009E73"

# 3. Apply mappings directly (these are more robust than mapVisualProperty for basic styles)
RCy3::setNodeShapeMapping('visual_group', 
                          c("SNP", "eQTL", "SNP & eQTL", "eQTL & SNP", "CP:KEGG_MEDICUS", "CP:REACTOME", "CP:WIKIPATHWAYS"),
                          c("ELLIPSE", "ELLIPSE", "ELLIPSE", "ELLIPSE", "RECTANGLE", "RECTANGLE", "RECTANGLE"), 
                          style.name = style_name)


RCy3::setNodeColorMapping('visual_group',
                          c("SNP", "eQTL", "SNP & eQTL", "eQTL & SNP", 
                            "CP:KEGG_MEDICUS", "CP:REACTOME", "CP:WIKIPATHWAYS"),
                          c("#FF7F00", "#377EB8", "#984EA3", "#984EA3", # Gene Colors (Orange, Blue, Purple)
                            "#E69F00", "#56B4E9", "#009E73"), 
                          mapping.type = "d",
                          style.name = style_name)


RCy3::setNodeColorMapping('type', 
                          c("SNP",    "eQTL",     "SNP & eQTL", "eQTL & SNP", "CP:KEGG_MEDICUS", "CP:REACTOME", "CP:WIKIPATHWAYS"), 
                          c("#FF7F00", "#377EB8", "#984EA3", "#984EA3", "#999999", "#E69F00", "#56B4E9", "#009E73"), 
                          mapping.type = "d",
                          style.name = style_name)

# Use the node-'label' as label
RCy3::setNodeLabelMapping('label', style.name = style_name)

# Set Edge style based on interaction type
RCy3::setEdgeLineStyleMapping('interaction', c("CP:WIKIPATHWAYS", "CP:REACTOME", "CP:KEGG_MEDICUS", "PPI"), 
                              c("DOT", "LONG_DASH", "EQUAL_DASH", "SOLID"), 
                              style.name = style_name)

RCy3::setEdgeColorMapping('interaction',
                          c("PPI", "CP:KEGG_MEDICUS", "CP:REACTOME", "CP:WIKIPATHWAYS"),
                          c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                          mapping.type = "d",
                          style.name = style_name)



# ======================================================================= #
# SAVING THE NETWORK IMAGE                                                #
# ======================================================================= #
Folder_of_your_choice <- rstudioapi::selectDirectory()
Folder_images_Cytoscape <- file.path(Folder_of_your_choice, "Folder_images_Cytoscape")
dir.create(Folder_images_Cytoscape, recursive = TRUE)
setwd(Folder_images_Cytoscape)

RCy3::exportImage()