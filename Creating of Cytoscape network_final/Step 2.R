# ============================================ #
#                   Step 0:                    #
#         Installing & Loading packages        #
# ============================================ #

Packages_needed <- c("tidyverse",
                     "biomaRt",
                     "STRINGdb",
                     "igraph",
                     "msigdbr",
                     "rstudioapi")

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

glimpse(mapped_genes_ENTREZ)

# Combining dataframes 
# X -> Ensembl IDs/genes (from FUMA) to | y -> bioMart databases
genes_combined <- merge(genes, 
                        mapped_genes_ENTREZ, 
                        by.x = "gene", 
                        by.y = "ensembl_gene_id")



genes_combined_current_biomart <- genes_combined
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
Biomart_hit_results <- biomaRt::getBM(attributes = c("ensembl_peptide_id", 
                                                     "ensembl_gene_id", 
                                                     "hgnc_symbol", 
                                                     "external_gene_name"), # maybe use hgnc_symbol???
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

Edges <- igraph::as_data_frame(network, what = "edges") # What-default is "edges"
Genes_PPI <- Edges %>% 
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


pathway_edges_preference <- pathways_raw %>%
  dplyr::filter(ensembl_gene %in% final_mapped_data$ensembl_gene_id) %>% # Keeps pathways in relation to genes
  dplyr::mutate(
    database_priority = dplyr::case_when(        # If 1 exists; if 2 exists else -> 3  
      gs_subcollection == "CP:WIKIPATHWAYS" ~ 1, # 1st preference
      gs_subcollection == "CP:REACTOME" ~ 2,     # 2nd/3rd preference
      gs_subcollection == "CP:KEGG_MEDICUS" ~ 2, # 2nd/3rd preference
      TRUE ~ 3                                   # 
    )
  ) %>%
  # For each gene, only the pathways from highest preference database available
  dplyr::group_by(ensembl_gene) %>% # group by gene
  dplyr::filter(database_priority == min(database_priority)) %>% # lowest value/(or highest preference is kept)
  
  dplyr::ungroup() %>%
  
  dplyr::select(source = ensembl_gene, target = gs_name, interaction = "gs_subcollection") # interaction -> is it from Wikipathways; or reactome or kegg?


tibble(pathway_edges_preference)
# ========================================================================================== # 
# TO BE LOOKED AT!
pathway_edges_preference_cleaned <- pathway_edges_preference %>%
  dplyr::mutate(target = stringr::str_remove_all(target, "KEGG_MEDICUS_|WP_|REACTOME_")) %>%
  dplyr::mutate(target = stringr::str_replace_all(target, "_", " "))
# ========================================================================================== #  

