# ============================================ #
#                   Step 0:                    #
#         Installing & Loading packages        #
# ============================================ #
Packages_needed <- c("tidyverse", 
                     "msigdbr", 
                     "biomaRt", 
                     "RCy3", 
                     "rWikiPathways", 
                     "rstudioapi", 
                     "STRINGdb", 
                     "Rgraphviz", 
                     "CePa", 
                     "clusterProfiler", 
                     "org.Hs.eg.db",
                     "DOSE") # "igraph")

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
# STEP 1: READ FUMA OUTPUT                                                #
#         1. Selecting Folder/data-directory                              #
#            And reading/opening files                                    #
#              For ANNOVAR and                                            #
#              For eQTL                                                   #
#         2. Filtering of Data                                            #
#              ANNOVAR and eQTL                                           #
#              Combining the two                                          #
# ======================================================================= #
# Part 1: Selecting Folder/data-directory                                 #
#   MAP TO ENTREZ GENE IDENTIFIERS                                        #
# ======================================================================= #
# You will need: annov.txt & eqtl.txt
# Make sure this script is in the same folder as your downloaded file


SNP2GENE_folder  <- rstudioapi::selectDirectory()
Folder_of_your_choice <- rstudioapi::selectDirectory()

Choosing_file_path <- function(file_path) {
  if(!dir.exists(file_path))
   return(NULL) 
      setwd(file_path)
    }



Choosing_file_path(Folder_of_your_choice)
Choosing_file_path(SNP2GENE_folder)

# Opening of the files:
annov <- read.table(paste0("annov.txt"), header = T, sep = "\t")            # Opens the file
eqtl <- read.table(paste0("eqtl.txt"), header = T, sep = "\t")              # Opens the file

# Without dyplr 
# ANNOVAR -> File 
annov.filt <- annov[annov$annot %in% c("intronic","exonic","UTR3","UTR5"),] # Filter based on "intronic";"exonic";"UTR3"; "UTR5"
annov.filt <- annov.filt[,c("gene", "symbol")]                              # Only selecting "gene" & "symbol"
annov.filt$Type <- "SNP"                                                    # Making a new column-> "Type" with only "SNP" in it

# eQTL-mapping -> File
eqtl.filt <- eqtl[,c("gene", "symbol")]                                     # Only selecting "gene" & "symbol"
eqtl.filt$Type <- "eQTL"                                                    # Making a new column-> Type" with only  "eQTL" in it


# With dyplr
# ANNOVAR -> File
annov.filt <- annov %>%
  filter(annot %in% c("intronic","exonic","UTR3","UTR5")) %>%               # Filter based on "intronic";"exonic";"UTR3"; "UTR5
  dplyr::select(gene, symbol) %>%                                           # Only selecting "gene" & "symbol
  mutate(Type = "SNP")                                       # Making a new columns -> "Memory" with only "PW & "Type" with only "SNP"
# eQTL-mapping -> File
eqtl.filt <- eqtl %>% 
  dplyr::select(gene, symbol) %>%                                           # Only selecting "gene" & "symbol"
  mutate(Type = "eQTL")                                      # Making a new columns -> "Memory" with only "PW &  Type" with only "eQTL" in it

# ============================================
# COMBINES MOST DELETRIOUS GENE (CONSEQUENCE) & EFFECT IN TISSUE
# ============================================
genes <- rbind(annov.filt, eqtl.filt, stringsAsFactors = FALSE)
# wanted to keep information if SNP or eQTL but code didn't work - need to fix
# genes_unique <- unique(genes[,c(1:4)]) # SEE THIS LATER



# ======================================================================= #
# STEP 2: Intergrating different databases:                               #
#         1. Biomart                                                      #
#              MAP TO ENTREZ GENE IDENTIFIERS                             #
#         2. STRING DB                                                    #
#              Mapping/making networks from STRING DB                     #
#         3. MSigDB (Gene Set Enrichtment Analysis?)                      #
#              Wikipathways; Reactome; KEGG                               #
#                                                                         #
# ======================================================================= #
# Part 1: Biomart                                                         #
#   MAP TO ENTREZ GENE IDENTIFIERS                                        #
# ======================================================================= #
# Making sure the "biomart" is correct
# ("ENSEMBL_MART_ENSEMBL" is (always) correct)
ensemble_mart_HS <- biomaRt::useMart(biomart = 
                                       "ENSEMBL_MART_ENSEMBL", 
                                     dataset = 
                                       "hsapiens_gene_ensembl")

# Mapping of the Ensemble ID to genes 
mapped_genes_ENTREZ <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                                     "entrezgene_id"),
                                      filters = "ensembl_gene_id",
                                      values = unique(genes$gene), 
                                      mart = ensemble_mart_HS)

# Combining dataframes 
# X -> Ensembl IDs/genes (from FUMA) to | y -> bioMart databases
genes_combined <- merge(genes, 
                        mapped_genes_ENTREZ, 
                        by.x = "gene", 
                        by.y = "ensembl_gene_id")

rm(mapped_genes_ENTREZ)
# ======================================================================= #
# STEP 2: Intergrating different databases:                               #
# ======================================================================= #
# Part 2: STRING DB                                                       #
#   Mapping/making networks from STRING DB                                #
# ======================================================================= #
# 1. Setting up string
string_db <- STRINGdb::STRINGdb$new(version="12",                # (Probably) the latest version
                                    species=9606,                # Is human
                                    score_threshold=900,         # Normal threshold is 400
                                    network_type="full",        
                                    link_data='combined_only',  
                                    input_directory="",         
                                    protocol="http")             # Otherwise it would not work 

mapped_genes_STRING <- string_db$map(genes_combined,
                                     "entrezgene_id",
                                     removeUnmappedRows = TRUE) 

# Seeing mapping success
nrow(mapped_genes_STRING) / nrow(genes_combined)

# 2. Map gene ids to protein ids
  # extract protein ids from the human network
mapped_genes_STRING$ensembl_peptide_id <- gsub("9606.", "", mapped_genes_STRING$STRING_id)

  # Mapping of the Ensemble ID to genes 
Biomart_hit_results <- biomaRt::getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id", "hgnc_symbol", "external_gene_name"), # maybe use hgnc_symbol???
                                      filters = "ensembl_peptide_id",
                                      values = unique(mapped_genes_STRING$ensembl_peptide_id), 
                                      mart = ensemble_mart_HS)

# 3. Combining dataframes 
# X -> Ensembl IDs/genes (from FUMA) to | y -> bioMart databases
final_mapped_data <- merge(mapped_genes_STRING, 
                           Biomart_hit_results, 
                           by = "ensembl_peptide_id", 
                           all.x = TRUE)

# 4. Plotting subnetwork for Edges-protein-protein-interactions
hits <- final_mapped_data$STRING_id
network <- string_db$get_subnetwork(hits)

Edges <- igraph::as_data_frame(network, what = "edges")
Genes_PPI <- Edges %>% # What-default is "edges"
  dplyr::rename(source = from, target = to) %>%
  dplyr::mutate(interaction = "PPI")

# ======================================================================= #
# For the edges                                                           #
# ======================================================================= #
# 0. Making-String-network ready for Biomart
Genes_PPI <- Genes_PPI %>%
  dplyr::mutate(ensembl_peptide_id_source = gsub("9606.", "", source), # ("9606\\.", "", 
                ensembl_peptide_id_target = gsub("9606.", "", target))


# 1. Cleaning peptide IDs for BioMart
ppi_peptides_ids <- unique(c(Genes_PPI$ensembl_peptide_id_source,
                             Genes_PPI$ensembl_peptide_id_target))


# Making mapping table from Biomart
mapping_table <- biomaRt::getBM(attributes = c("ensembl_peptide_id", 
                                               "ensembl_gene_id", 
                                               "entrezgene_id", 
                                               "external_gene_name",
                                               "hgnc_symbol"), # maybe use hgnc_symbol???
                                filters = "ensembl_peptide_id",
                                values = ppi_peptides_ids, 
                                mart = ensemble_mart_HS)

# 4. Filtering gene symbols
mapping_table <- mapping_table %>%
  dplyr::mutate(hgnc_symbol = na_if(hgnc_symbol, ""),                            # If hgnc_symbol is empty -> NA
                external_gene_name = na_if(external_gene_name, ""),              # If external_gene_name is empty -> NA
                # The actual fallback: pick the first non-NA value
                final_label = coalesce(hgnc_symbol, external_gene_name, ensembl_gene_id))    # First non-NA -> is chosen


glimpse(mapping_table)
# Putting back the edges-ppi-table with gene and peptide-id
# Ensembl peptide id to gene id 
ensembl_gene_peptide_id <- mapping_table$ensembl_gene_id
names(ensembl_gene_peptide_id) <- mapping_table$ensembl_peptide_id

# Ensembl peptide id to gene symbol 
peptide_id_to_symbol <- mapping_table$final_label
names(peptide_id_to_symbol) <- mapping_table$ensembl_peptide_id


Edges_PPI_final <- Genes_PPI %>%
  dplyr::mutate(source_gene_id = unname(ensembl_gene_peptide_id[ensembl_peptide_id_source]),
                target_gene_id = unname(ensembl_gene_peptide_id[ensembl_peptide_id_target]),
                source_peptide_id = ensembl_peptide_id_source,
                target_peptide_id = ensembl_peptide_id_target,
                source_symbol = unname(peptide_id_to_symbol[ensembl_peptide_id_source]),
                target_symbol = unname(peptide_id_to_symbol[ensembl_peptide_id_target]),
                source = source_gene_id,
                target = target_gene_id)

glimpse(Edges_PPI_final_test)
glimpse(Edges_PPI_final)

glimpse(all_nodes)
# ======================================================================= #
# For the Nodes                                                           #
# ======================================================================= #
Nodes_Genes_Grouped <- final_mapped_data %>% 
  dplyr::group_by(ensembl_gene_id) %>% # groups by 
  dplyr::summarise(                    # Makes every group a row
    label = dplyr::first(                     # Choose a label/symbol
      coalesce(                        # Picks 1st non-NA
        na_if(hgnc_symbol, ""),        # 1. If empty -> NA (& not picked)
        na_if(external_gene_name, ""), # 2. If empty -> NA (& not picked)
        ensembl_gene_id)),             # 3. Last option (if both are missing)
    
    entrez_id = dplyr::first(entrezgene_id),  # Because it needs to know what it needs 'summarize'
    string_id = dplyr::first(STRING_id),      # Because it needs to know what it needs 'summarize'
    hgnc_name = dplyr::first(hgnc_symbol),    # Because it needs to know what it needs 'summarize'
    ensembl_peptide = dplyr::first(ensembl_peptide_id),
    type = paste(unique(Type), collapse = " & ")) %>% # If it's both SNP and eQTL -> combined
  dplyr::ungroup() %>%
  dplyr::rename(id = ensembl_gene_id)


glimpse(Nodes_Genes_Grouped)
# ======================================================================= #
# STEP 2: Integrating different databases:                                #
# ======================================================================= #
# Part 3: MSigDB                                                          #
#   Wikipathways; Reactome; KEGG                                          #
# ======================================================================= #
pathways_raw <- msigdbr(species = "Homo sapiens", collection = "C2") %>%
  dplyr::filter(gs_subcollection %in% c("CP:KEGG_MEDICUS", "CP:REACTOME", "CP:WIKIPATHWAYS")) %>%
  dplyr::select(gs_name, gs_subcollection, ncbi_gene, ensembl_gene, gene_symbol)


pathway_edges_preference <- pathways_raw_priority %>%
  dplyr::filter(ensembl_gene %in% final_mapped_data$ensembl_gene_id) %>%
  dplyr::mutate(
    database_priority = dplyr::case_when(
      gs_subcollection == "CP:WIKIPATHWAYS" ~ 1,
      gs_subcollection == "CP:REACTOME" ~ 2,
      gs_subcollection == "CP:KEGG_MEDICUS" ~ 2, # 3?
      TRUE ~ 3
    )
  ) %>%
  # For every gene, keep only the pathways from the highest priority DB available
  dplyr::group_by(ensembl_gene) %>%
  dplyr::filter(database_priority == min(database_priority)) %>% 
  
  dplyr::slice_head(n = 3) %>%
  
  dplyr::ungroup() %>%
  dplyr::select(source = ensembl_gene, target = gs_name, interaction = "gs_subcollection")


# ========================================================================================== # 
# TO BE LOOKED AT!
pathway_edges_preference_cleaned <- pathway_edges_preference %>%
  dplyr::mutate(target = stringr::str_remove_all(target, "KEGG_MEDICUS_|WP_|REACTOME_")) %>%
  dplyr::mutate(target = stringr::str_replace_all(target, "_", " "))
# ========================================================================================== #  

# ======================================================================= #
# STEP 3: Enrichment Analysis                                             #
# ======================================================================= #
# Part 1: Gene Ontology (GO)                                              #
#           GO over-representation analysis                               #
# Part 2: DisGeNET                                                        #
#   Ma                                                                    #
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
# Disgenet
dgn_enrich_results <- enrichDGN(gene = unique(Nodes_Genes_Grouped$entrez_id),
                                pAdjustMethod = "fdr", # FDR -> Benjamini–Hochberg? (BH)?
                                pvalueCutoff  = 0.05)

glimpse(Nodes_Genes_Grouped)
# ======================================================================= #
# Making plots
# Make sure this script is in the same folder as your downloaded file
DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
DATA.DIR - dirname(rstudioapi::getActiveDocumentContext()$path)
Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
dir.create("Folder_images")
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
  dplyr::mutate(label = id, type = "Pathway", bio_weight = 0)

# ======================================================================= #
# Part 2: Nodes
#   Gene-Nodes
# ======================================================================= #
node_genes <- Nodes_Genes_Grouped %>%
  dplyr::left_join(GO_count, by = c("label" = "geneID"))  %>%
  dplyr::mutate(GO_count = tidyr::replace_na(gene_count, 0))  %>%
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



# ======================================================================= #
# STEP 4: Launching Cytoscape: 
# ======================================================================= #
RCy3::cytoscapePing()


style_name <- "FUMA_Scientific_Style_test"

RCy3::createVisualStyle(style_name)
style_name %in% getVisualStyleNames()


RCy3::cytoscapeVersionInfo()
RCy3::commandsGET("cytoscape version")

RCy3::createNetworkFromDataFrames(nodes = all_nodes, 
                                  edges = all_edges, 
                                  title = "FUMA_Enriched_Network_test")

str(all_nodes)
str(all_edges)


style_name <- "FUMA_Scientific_Style_test"
RCy3::createVisualStyle(style_name, defaults)
style_name %in% RCy3::getVisualStyleNames(style.name, defaults)

RCy3::setNodeSizeMapping('bio_impact', 
                         table.column.values = c(0, max(node_pathways$bio_weight)), 
                         sizes = c(30, 100), 
                         mapping.type = "c",
                         style.name = "FUMA_Scientific_Style_test") # 'c' for continuous mapping



setVisualStyle(style_name, network = "FUMA_Enriched_Network")
# 3. Apply mappings directly (these are more robust than mapVisualProperty for basic styles)
RCy3::setNodeShapeMapping('type', c("SNP", "eQTL", "Pathway"), c("ELLIPSE", "ELLIPSE", "RECTANGLE"), style.name = "FUMA_Scientific_Style_test")
RCy3::setNodeColorMapping('type', 
                          c("SNP", "eQTL", "SNP & eQTL", "Pathway"), 
                          c("#FF7F00", "#377EB8", "#984EA3", "#4DAF4A"), 
                          mapping.type = "d",
                          style.name = style_name)

# Use the node-'label' as label
RCy3::setNodeLabelMapping('label', style.name = "FUMA_Scientific_Style_test")

# Set Edge style based on interaction type
RCy3::setEdgeLineStyleMapping('interaction', c("Pathway-Gene", "PPI"), c("DOT", "SOLID"), style.name = "FUMA_Scientific_Style_test")

