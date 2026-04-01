# ======================================================================= #
#                           STEP 0: Setup                                 #
#         Installing & Loading packages, Setting Directories              #
# ======================================================================= #
Packages_needed <- c("tidyverse", "msigdbr", "biomaRt", "RCy3", 
                     "rWikiPathways", "rstudioapi", "STRINGdb", 
                     "Rgraphviz", "CePa", "clusterProfiler", 
                     "org.Hs.eg.db", "DOSE", "devtools", "pak",
                     "enrichplot") # <--- whops 

install_and_load <- function(package_names) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
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
}

install_and_load(Packages_needed)

# Install disgenet2r if missing
if (!requireNamespace("disgenet2r", quietly = TRUE)) {
  devtools::install_gitlab("medbio/disgenet2r")
}
library(disgenet2r)

# Setup Directories
DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(DATA.DIR)

out_dir  <- file.path(DATA.DIR, "output")
img_dir  <- file.path(DATA.DIR, "Folder_images")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(img_dir)) dir.create(img_dir, recursive = TRUE)

message("Step 0 Complete: Packages and Directories Set!")

# ======================================================================= #
# STEP 1: READ FUMA OUTPUT                                                #
# ======================================================================= #
# Ensure annov.txt and eqtl.txt are in your DATA.DIR!
annov <- read.table("annov.txt", header = TRUE, sep = "\t")
eqtl <- read.table("eqtl.txt", header = TRUE, sep = "\t")

# Filter ANNOVAR
annov.filt <- annov %>%
  dplyr::filter(annot %in% c("intronic","exonic","UTR3","UTR5")) %>%
  dplyr::select(gene, symbol) %>%
  dplyr::mutate(Type = "SNP")

# Filter eQTL
eqtl.filt <- eqtl %>% 
  dplyr::select(gene, symbol) %>%
  dplyr::mutate(Type = "eQTL")

# Combine them
genes <- rbind(annov.filt, eqtl.filt, stringsAsFactors = FALSE)
key_genes_list <- unique(genes$gene) # Saving this for Step 5!

message("Step 1 Complete: FUMA Data Loaded!")

# ======================================================================= #
# STEP 2: Integrating Databases (Biomart, STRING, MSigDB)                 #
# ======================================================================= #
# 1. Biomart Mapping
ensemble_mart_HS <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
mapped_genes_ENTREZ <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                                      filters = "ensembl_gene_id", values = unique(genes$gene), 
                                      mart = ensemble_mart_HS)
genes_combined <- merge(genes, mapped_genes_ENTREZ, by.x = "gene", by.y = "ensembl_gene_id")

# 2. STRING DB Mapping
string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=900, network_type="full", input_directory="")
mapped_genes_STRING <- string_db$map(genes_combined, "entrezgene_id", removeUnmappedRows = TRUE) 
mapped_genes_STRING$ensembl_peptide_id <- gsub("9606\\.", "", mapped_genes_STRING$STRING_id)

Biomart_hit_results <- biomaRt::getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id", "hgnc_symbol", "external_gene_name"),
                                      filters = "ensembl_peptide_id", values = unique(mapped_genes_STRING$ensembl_peptide_id), 
                                      mart = ensemble_mart_HS)

final_mapped_data <- merge(mapped_genes_STRING, Biomart_hit_results, by = "ensembl_peptide_id", all.x = TRUE)

# Base FUMA PPI Edges
hits <- final_mapped_data$STRING_id
network <- string_db$get_subnetwork(hits)
Genes_PPI <- igraph::as_data_frame(network, what = "edges") %>% 
  dplyr::rename(source = from, target = to) %>%
  dplyr::mutate(interaction = "PPI")

Genes_PPI <- Genes_PPI %>%
  dplyr::mutate(ensembl_peptide_id_source = gsub("9606\\.", "", source),
                ensembl_peptide_id_target = gsub("9606\\.", "", target))

ppi_peptides_ids <- unique(c(Genes_PPI$ensembl_peptide_id_source, Genes_PPI$ensembl_peptide_id_target))

mapping_table <- biomaRt::getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id", "entrezgene_id", "external_gene_name", "hgnc_symbol"),
                                filters = "ensembl_peptide_id", values = ppi_peptides_ids, mart = ensemble_mart_HS) %>%
  dplyr::mutate(hgnc_symbol = na_if(hgnc_symbol, ""), external_gene_name = na_if(external_gene_name, ""),
                final_label = coalesce(hgnc_symbol, external_gene_name, ensembl_gene_id))

ensembl_gene_peptide_id <- setNames(mapping_table$ensembl_gene_id, mapping_table$ensembl_peptide_id)
peptide_id_to_symbol <- setNames(mapping_table$final_label, mapping_table$ensembl_peptide_id)

Edges_PPI_final <- Genes_PPI %>%
  dplyr::mutate(source_gene_id = unname(ensembl_gene_peptide_id[ensembl_peptide_id_source]),
                target_gene_id = unname(ensembl_gene_peptide_id[ensembl_peptide_id_target]),
                source_symbol = unname(peptide_id_to_symbol[ensembl_peptide_id_source]),
                target_symbol = unname(peptide_id_to_symbol[ensembl_peptide_id_target])) 

# FUMA Nodes
Nodes_Genes_Grouped <- final_mapped_data %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::summarise(
    label = dplyr::first(coalesce(na_if(hgnc_symbol, ""), na_if(external_gene_name, ""), ensembl_gene_id)),
    entrez_id = dplyr::first(entrezgene_id), 
    string_id = dplyr::first(STRING_id), 
    hgnc_name = dplyr::first(hgnc_symbol),
    ensembl_peptide = dplyr::first(ensembl_peptide_id), 
    type = paste(unique(Type), collapse = " & ")) %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(id = ensembl_gene_id)

# 3. MSigDB Pathways
pathways_raw <- msigdbr(species = "Homo sapiens", collection = "C2") %>%
  dplyr::filter(gs_subcollection %in% c("CP:KEGG_MEDICUS", "CP:REACTOME", "CP:WIKIPATHWAYS"))

pathway_edges_preference <- pathways_raw %>%
  dplyr::filter(ensembl_gene %in% final_mapped_data$ensembl_gene_id) %>%
  dplyr::mutate(db_preference = dplyr::case_when(gs_subcollection == "CP:WIKIPATHWAYS" ~ 1, gs_subcollection == "CP:REACTOME" ~ 2, TRUE ~ 3)) %>%
  dplyr::group_by(ensembl_gene) %>% dplyr::filter(db_preference == min(db_preference)) %>% dplyr::ungroup() %>%
  dplyr::select(source = ensembl_gene, target = gs_name, interaction = gs_subcollection)

message("Step 2 Complete: Database Mapping Done!")
# ======================================================================= #
# STEP 3: Enrichment Analysis (GO & General DisGeNET)                     #
# ======================================================================= #
enrichGO_results <- enrichGO(gene = unique(Nodes_Genes_Grouped$entrez_id), 
                             OrgDb = org.Hs.eg.db, keyType = "ENTREZID", 
                             ont = "BP", pAdjustMethod = "BH", 
                             pvalueCutoff = 0.05, readable = TRUE)

dgn_enrich_results <- enrichDGN(gene = unique(Nodes_Genes_Grouped$entrez_id), 
                                pAdjustMethod = "fdr", pvalueCutoff = 0.05, 
                                readable = TRUE)

# Make plots
setwd(img_dir)
png("GO_Treeplot.png", width = 3000, height = 4000, res = 150)
plot(treeplot(pairwise_termsim(enrichGO_results), nCluster = 5, showCategory = 30))
dev.off()

png("DGN_Treeplot.png", width = 3000, height = 4000, res = 150)
plot(treeplot(pairwise_termsim(dgn_enrich_results), nCluster = 5, showCategory = 20))
dev.off()
setwd(DATA.DIR)

# --- PERMANENT FIX 1: Using tidyr::separate_rows ---
gene_count_calculation = function(res) {
  if(is.null(res)) return(NULL)
  
  as.data.frame(res) %>% 
    dplyr::filter(p.adjust < 0.05) %>% 
    tidyr::separate_rows(geneID, sep ="/") %>% 
    dplyr::group_by(geneID) %>% 
    dplyr::summarize(gene_count = n())
}

GO_score <- gene_count_calculation(enrichGO_results)
DGN_score <- gene_count_calculation(dgn_enrich_results)

# ======================================================================= #
# STEP 4: Cytoscape Network 1 (General FUMA Enrichment)                   #
# ======================================================================= #
RCy3::cytoscapePing()

# --- PERMANENT FIX 2: Matching Edge IDs and converting to Base R Dataframes ---

# 1. Prepare and fix the Edges (Swap Peptide IDs for Gene IDs)
edges_ppi_fixed <- Edges_PPI_final %>%
  dplyr::select(source = source_gene_id, target = target_gene_id, interaction)

pathway_edges_clean <- pathway_edges_preference %>%
  dplyr::select(source, target, interaction) 

all_edges_fixed <- dplyr::bind_rows(edges_ppi_fixed, pathway_edges_clean) %>% 
  dplyr::filter(!is.na(source) & !is.na(target)) %>%
  dplyr::distinct(source, target, .keep_all = TRUE)

# 2. Prepare the Nodes
node_pathways <- data.frame(id = unique(pathway_edges_clean$target)) %>% 
  dplyr::mutate(label = id, type = "Pathway", bio_weight = 0)

node_genes <- Nodes_Genes_Grouped %>%
  dplyr::left_join(GO_score, by = c("label" = "geneID")) %>% 
  dplyr::mutate(GO_count = tidyr::replace_na(gene_count, 0)) %>% 
  dplyr::select(-gene_count) %>%
  dplyr::left_join(DGN_score, by = c("hgnc_name" = "geneID")) %>% 
  dplyr::mutate(DGN_count = tidyr::replace_na(gene_count, 0)) %>% 
  dplyr::select(-gene_count) %>%
  dplyr::mutate(bio_impact = GO_count + DGN_count)

all_nodes_fixed <- dplyr::bind_rows(node_genes, node_pathways) %>% 
  dplyr::filter(!is.na(id)) %>%
  dplyr::distinct(id, .keep_all = TRUE)

# 3. Force conversion to plain dataframes and characters for RCy3
all_nodes_fixed <- as.data.frame(all_nodes_fixed)
all_edges_fixed <- as.data.frame(all_edges_fixed)
all_nodes_fixed$id <- as.character(all_nodes_fixed$id)
all_edges_fixed$source <- as.character(all_edges_fixed$source)
all_edges_fixed$target <- as.character(all_edges_fixed$target)

# 4. Build Network 1 in Cytoscape
RCy3::createNetworkFromDataFrames(nodes = all_nodes_fixed, 
                                  edges = all_edges_fixed, 
                                  title = "FUMA_Enriched_Network")

# 5. Apply Visual Styling
style_name <- "FUMA_Scientific_Style"
if(!(style_name %in% RCy3::getVisualStyleNames())) {
  RCy3::createVisualStyle(style_name, defaults = list(NODE_SHAPE="ELLIPSE"))
}

RCy3::setVisualStyle(style_name, network = "FUMA_Enriched_Network")
RCy3::setNodeSizeMapping('bio_impact', c(0, max(all_nodes_fixed$bio_impact, na.rm = TRUE)), c(30, 100), mapping.type = "c", style.name = style_name)
RCy3::setNodeColorMapping('type', c("SNP", "eQTL", "SNP & eQTL", "Pathway"), c("#FF7F00", "#377EB8", "#984EA3", "#4DAF4A"), mapping.type = "d", style.name = style_name)
RCy3::setNodeLabelMapping('label', style.name = style_name)
RCy3::layoutNetwork("force-directed")

message("Step 4 Complete: Cleaned FUMA Network generated & styled in Cytoscape!")

# ======================================================================= #
# STEP 5: Cytoscape Network 2 (Targeted Diabetes PPI)                     #
# ======================================================================= #

# 1. Fetch Diabetes Disease Genes via DisGeNET API
api_key <- "c036b920-6851-4f1a-8957-992bb8e50e35"
Sys.setenv(DISGENET_API_KEY = api_key)

message("Fetching Diabetes genes from DisGeNET...")
diabetes_disgenet <- disease2gene(disease = "UMLS_C0011860", database = "CURATED", score = c(0, 1))
diabetes_genes <- data.frame(gene = unique(diabetes_disgenet@qresult[["gene_symbol"]]), 
                             score = diabetes_disgenet@qresult[["score"]]) 

# 2. Convert HGNC to Ensembl using biomaRt (GRCh37)
message("Mapping DisGeNET genes to Ensembl...")
grch37.gene <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                       dataset = "hsapiens_gene_ensembl", 
                       host = "https://useast.ensembl.org") # Using the US East mirror
disgenet_ensembl <- getBM(mart = grch37.gene, values = diabetes_genes$gene, filters = "hgnc_symbol", attributes = c("ensembl_gene_id", "hgnc_symbol"))

# 3. Combine HD Genes with FUMA Key Genes (from Step 1)
combined_gene_dataframe <- data.frame(
  gene = c(disgenet_ensembl$ensembl_gene_id, key_genes_list),
  groups = c(rep("Diabetes genes", nrow(disgenet_ensembl)), rep("FUMA Key genes", length(key_genes_list)))
) %>% dplyr::distinct(gene, .keep_all = TRUE)

# 4. Map STRING to Ensembl (WITH CHUNKING TO PREVENT 504 TIMEOUT)
message("Extracting STRINGdb network...")
g_all <- string_db$get_graph()
ppi_all <- igraph::as_data_frame(g_all, what = "edges")

# Extract protein IDs from BOTH 'from' and 'to' columns
protein_ids <- unique(c(sapply(strsplit(ppi_all$from, '\\.'), function(x) x[2]), 
                        sapply(strsplit(ppi_all$to, '\\.'), function(x) x[2])))

message("Querying Ensembl in chunks to prevent server timeout... This may take a minute.")
chunk_size <- 1000
protein_chunks <- split(protein_ids, ceiling(seq_along(protein_ids)/chunk_size))

# Loop through chunks to safely download mapping
mart_results_list <- lapply(seq_along(protein_chunks), function(i) {
  cat("Processing chunk", i, "of", length(protein_chunks), "...\n")
  biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
                 filters = "ensembl_peptide_id", 
                 values = protein_chunks[[i]], 
                 mart = grch37.gene)
})

# Combine mapping results
mart_results <- dplyr::bind_rows(mart_results_list) %>% 
  dplyr::distinct(ensembl_peptide_id, .keep_all = TRUE)

# 5. Apply mapping to the PPI network
message("Filtering the PPI Network...")
ppi_all$from_clean <- sapply(strsplit(ppi_all$from, '\\.'), function(x) x[2])
ppi_all$to_clean <- sapply(strsplit(ppi_all$to, '\\.'), function(x) x[2])

ppi_all$from_gene_id <- mart_results$ensembl_gene_id[match(ppi_all$from_clean, mart_results$ensembl_peptide_id)]
ppi_all$to_gene_id <- mart_results$ensembl_gene_id[match(ppi_all$to_clean, mart_results$ensembl_peptide_id)]

# Remove NAs mapped to genes
ppi_all <- ppi_all[complete.cases(ppi_all$from_gene_id, ppi_all$to_gene_id), ]

# 6. Filter network for Diabetes + Key Genes
target_ids <- combined_gene_dataframe$gene
ppi_diabetes_filtered <- ppi_all %>%
  dplyr::filter(from_gene_id %in% target_ids | to_gene_id %in% target_ids) %>%
  dplyr::filter(combined_score >= 900)

# 7. Prepare Nodes and Edges for Cytoscape
message("Preparing data with Association Scores...")

# Create the base node list from the PPI network
nodes_df <- data.frame(id = unique(c(ppi_diabetes_filtered$from_gene_id, ppi_diabetes_filtered$to_gene_id)), 
                       stringsAsFactors = FALSE)

# Join with your group definitions (FUMA vs Diabetes)
nodes_df <- merge(nodes_df, combined_gene_dataframe, by.x = "id", by.y = "gene", all.x = TRUE)

# --- NEW CALCULATION: Mapping DisGeNET Scores & Evidence ---
# We merge the HGNC symbol mapping with the original score data from Part 1
score_lookup <- merge(disgenet_ensembl, diabetes_genes, by.x = "hgnc_symbol", by.y = "gene")

# Join the scores and evidence into our main node table
nodes_df <- merge(nodes_df, score_lookup[,c("ensembl_gene_id", "score")], 
                  by.x = "id", by.y = "ensembl_gene_id", all.x = TRUE)

# Final cleanup for the table
nodes_df$groups[is.na(nodes_df$groups)] <- "Other Interactors"
nodes_df$score[is.na(nodes_df$score)] <- 0  # Genes not in DisGeNET get a 0 score

edges_df <- data.frame(source = ppi_diabetes_filtered$from_gene_id, 
                       target = ppi_diabetes_filtered$to_gene_id,
                       interaction = "PPI", 
                       weight = ppi_diabetes_filtered$combined_score, 
                       stringsAsFactors = FALSE)

# --- STRICT RCy3 FORMATTING FIXES ---
nodes_df <- as.data.frame(nodes_df[!is.na(nodes_df$id), ])
edges_df <- as.data.frame(edges_df[!is.na(edges_df$source) & !is.na(edges_df$target), ])

nodes_df$id <- as.character(nodes_df$id)
edges_df$source <- as.character(edges_df$source)
edges_df$target <- as.character(edges_df$target)

# 8. Send to Cytoscape and Apply Visual Scoring
message("Building Network in Cytoscape...")
RCy3::cytoscapePing()

# Create the Network
RCy3::createNetworkFromDataFrames(nodes = nodes_df, edges = edges_df, 
                                  title = "Diabetes-PPI-Network", collection = "DisGeNET")

# --- VISUAL CALCULATION: Sizing nodes by Association Score ---
# This makes the "High Evidence" genes physically larger in your network
RCy3::setNodeSizeMapping('score', 
                         table.column.values = c(0, 1), 
                         sizes = c(30, 100), 
                         mapping.type = "c")

RCy3::layoutNetwork("force-directed")

# Load the full table (including scores) into Cytoscape's Node Table
RCy3::loadTableData(nodes_df, data.key.column = "id", table.key.column = "name", table = "node")

message("Step 5 Complete: Targeted Diabetes Network with Scores generated!")