# ======================================================================= #
# STEP 4: Preparation for And launching 
#         Cytoscape: 
#   1. Finishing the Edges
#   2. Finishing the Nodes
#   3. Creating the network
# ======================================================================= #
# Part 1: Edges
#   Pathway-Edges
# ======================================================================= #
pathway_edges_clean <- pathway_edges_preference %>%
  dplyr::select(source = source,
                target = target,       # This is the Pathway Name
                interaction)  %>%      # Pathway name acts as its own symbol
  dplyr::filter(!is.na(source) & !is.na(target))


# ======================================================================= #
# Combining the two edges
# ======================================================================= #
all_edges <- bind_rows(Edges_PPI_final, pathway_edges_clean) %>%
  distinct(source, target, .keep_all = TRUE)


glimpse(all_edges)
# ======================================================================= #
# STEP 4: Launching Cytoscape:
#   2. Finishing the Nodes
# ======================================================================= #
# Part 2: Nodes
#   Pathway-Nodes
# ======================================================================= #
node_pathways <- data.frame(id = unique(pathway_edges_clean$target)) %>%
  dplyr::mutate(label = id, type = "Pathway", bio_impact = 0)

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
  dplyr::mutate(bio_impact = GO_count + DGN_count)  %>%
  as.data.frame()

# ======================================================================= #
# Combining the two nodes
# ======================================================================= #
all_nodes <- bind_rows(node_genes, node_pathways) %>% 
  dplyr::distinct(id, .keep_all = TRUE)


glimpse(all_nodes)
# ======================================================================= #
# STEP 4: Launching Cytoscape: 
# ======================================================================= #
RCy3::cytoscapePing()


style_name <- "FUMA_Scientific_Style_test"

RCy3::createVisualStyle(style_name)
style_name %in% getVisualStyleNames()


RCy3::createNetworkFromDataFrames(nodes = all_nodes, 
                                  edges = all_edges, 
                                  title = "FUMA_Enriched_Network_test")


style_name <- "FUMA_Scientific_Style_test"
RCy3::createVisualStyle(style_name, defaults)
style_name %in% RCy3::getVisualStyleNames()
setVisualStyle(style_name, network = "FUMA_Enriched_Network_test")


RCy3::setNodeSizeMapping('bio_impact', 
                         table.column.values = c(0, max(node_pathways$bio_impact)), 
                         sizes = c(30, 100), 
                         mapping.type = "c",
                         style.name = style_name) # 'c' for continuous mapping




# 3. Apply mappings directly (these are more robust than mapVisualProperty for basic styles)
RCy3::setNodeShapeMapping('type', c("SNP", "eQTL", "Pathway"), c("ELLIPSE", "ELLIPSE", "RECTANGLE"), style.name = "FUMA_Scientific_Style_test")
RCy3::setNodeColorMapping('type', 
                          c("SNP",    "eQTL",     "SNP & eQTL", "Pathway"), 
                          c("#FF7F00", "#377EB8", "#984EA3",    "#4DAF4A"), 
                          mapping.type = "d",
                          style.name = style_name)

# Use the node-'label' as label
RCy3::setNodeLabelMapping('label', style.name = "FUMA_Scientific_Style_test")

# Set Edge style based on interaction type
RCy3::setEdgeLineStyleMapping('interaction', c("CP:WIKIPATHWAYS", "CP:REACTOME", "CP:KEGG_MEDICUS", "PPI"), 
                              c("DOT", "LONG_DASH", "EQUAL_DASH", "SOLID"), 
                              style.name = style_name)
