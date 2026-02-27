#------------------------------------------------------------------------------#
#                                                                              #
#                                                                              #
#                                                                              #
#                SCRIPT TO COLLECT AND FILTER FUMA RESULTS                     #
#               SECOND PART IS FOR CYTOSCAPE NETWORK ANALYSIS                  #
#                                                                              #
#                                                                              #
#------------------------------------------------------------------------------#

# ============================================
# SETUP
# ============================================


# I tried making it into one go
# You will need: 
# "tidyverse", "msigdbr", "biomaRt", "RCy3", "rWikiPathways", "rstudioapi
Packages_needed <- c("tidyverse", "msigdbr", "biomaRt", "RCy3", "rWikiPathways", "rstudioapi")
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

# Data-prep SNP2GENE
# You will need: annov.txt & eqtl.txt
# Make sure this script is in the same folder as your downloaded file
DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(DATA.DIR)

# Choosing DATA.DIR SNP2GENE 
DATA.DIR_SNP2GENE  <- rstudioapi::selectDirectory()
setwd(DATA.DIR_SNP2GENE)
# Choosing DATA.DIR GENE2FUNC
DATA.DIR_GENE2FUNC <- rstudioapi::selectDirectory()

# Choosing a miscellaneous folder  
DATA.DIR_MISC <- rstudioapi::selectDirectory()
setwd(DATA.DIR_SNP2GENE)



# ============================================
# READ FUMA OUTPUT
# ============================================


file_path <- file.path(DATA.DIR_SNP2GENE)
read.table(DATA.DIR_SNP2GENE, "FUMA_iwrd", "annov.txt", header = TRUE)

# Full results of ANNOVAR -> "annov"
# ANNOVAR annotates "consequence" function by prioritising the most deleterious annotation for SNPs which are locating a genomic region where multiple genes are overlapped.

# eQTL mapping -> "eqtl" 
# Contains unique pair of SNP-gene-tissue

# Penn word 
annov <- read.table(paste0("FUMA_iwrd/annov.txt"), header = T, sep = "\t")  # Opens the file
annov.filt <- annov[annov$annot %in% c("intronic","exonic","UTR3","UTR5"),] # Filter based on "intronic";"exonic";"UTR3"; "UTR5"
annov.filt <- annov.filt[,c("gene", "symbol")]                              # Only selecting "gene" & "symbol"
annov.filt$Memory <- "PW"                                                   # Making a new column-> "Memory" with only "PW" in it
annov.filt$Type <- "SNP"                                                    # Making a new column-> "Type" with only "SNP" in it


eqtl <- read.table(paste0("FUMA_iwrd/eqtl.txt"), header = T, sep = "\t")    # Opens the file
eqtl.filt <- eqtl[,c("gene", "symbol")]                                     # Only selecting "gene" & "symbol"
eqtl.filt$Memory <- "PW"                                                    # Making a new column-> "Memory" with only "PW" in it
eqtl.filt$Type <- "eQTL"                                                    # Making a new column-> Type" with only  "eQTL" in it

# Penn word (Max)
annov <- read.table(paste0("annov.txt"), header = T, sep = "\t")            # Opening the file

annov.filt <- annov %>%
  filter(annot %in% c("intronic","exonic","UTR3","UTR5")) %>%               # Filter based on "intronic";"exonic";"UTR3"; "UTR5
  dplyr::select(gene, symbol) %>%                                           # Only selecting "gene" & "symbol
  mutate(Memory = "PW", Type = "SNP")                                       # Making a new columns -> "Memory" with only "PW & "Type" with only "SNP"

eqtl <- read.table(paste0("eqtl.txt"), header = T, sep = "\t")              # Opening the file

eqtl.filt <- eqtl %>% 
  dplyr::select(gene, symbol) %>%                                           # Only selecting "gene" & "symbol"
  mutate(Memory = "PW", Type = "eQTL")                                      # Making a new columns -> "Memory" with only "PW &  Type" with only "eQTL" in it

# ============================================
# COMBINE DIFFERENT MEMORY TESTS
# ============================================

genes <- rbind(annov.filt, eqtl.filt, stringsAsFactors = FALSE)
# wanted to keep information if SNP or eQTL but code didn't work - need to fix
genes_test <- unique(genes[,c(1:3)]) # SEE THIS LATER

# ============================================
# MAP TO ENTREZ GENE IDENTIFIERS
# ============================================


# Making sure the "biomart" is correct
  # ("ENSEMBL_MART_ENSEMBL" is correct)
ensemble_mart_HS <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# Mapping of the Ensemble ID to genes 
mapping <- biomaRt::getBM(attributes = c("ensembl_gene_id","entrezgene_id"), filters = "ensembl_gene_id",
                          values = unique(genes$gene), 
                          mart = ensemble_mart_HS)
glimpse(mapping)

genes <- merge(genes, mapping, by.x = "gene", by.y = "ensembl_gene_id")


# ============================================
# LOAD PATHWAY INFO
# ============================================

#Wikipathways in Cyto?
installApp('WikiPathways')

# For making the a new directory called: "databases"
dir.create("databases")

# Install the Files
# Link for all 
  # https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/

# Kegg -> 
  # https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c2.cp.kegg_medicus.v2023.2.Hs.entrez.gmt
kegg <- rWikiPathways::readGMT(file.path(getwd(), "databases", "c2.cp.kegg_medicus.v2023.2.Hs.entrez.gmt"))
# Reactome -> 
  # https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c2.cp.reactome.v2023.2.Hs.entrez.gmt
reactome <- rWikiPathways::readGMT(file.path(getwd(), "databases", "c2.cp.reactome.v2023.2.Hs.entrez.gmt"))
# WikiPathways ->
  # https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c2.cp.wikipathways.v2023.2.Hs.entrez.gmt
wp <- rWikiPathways::readGMT(file.path(getwd(), "databases", "c2.cp.wikipathways.v2023.2.Hs.entrez.gmt"))

# Combinging results in 1 factor. & removing
pathways <- rbind(kegg, reactome, wp)
rm(kegg, reactome, wp)

# OR
# Use "msigdbr"-package
Canonical_Pathways_humans <- msigdbr(species = "Homo sapiens", collection = "C2")

kegg_v2 <- Canonical_Pathways_humans %>%
  dplyr::filter(gs_subcollection == "CP:KEGG_MEDICUS") %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  dplyr::rename(gene = ncbi_gene, term = gs_name)

reactome_v2 <- Canonical_Pathways_humans %>%
  dplyr::filter(gs_subcollection == "CP:REACTOME") %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  dplyr::rename(gene = ncbi_gene, term = gs_name)

wp_v2 <- Canonical_Pathways_humans %>%
  dplyr::filter(gs_subcollection == "CP:WIKIPATHWAYS") %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  dplyr::rename(gene = ncbi_gene, term = gs_name)

# Combinging results in 1 factor. & removing the rest
pathways <- rbind(kegg_v2, reactome_v2, wp_v2)
rm(kegg_v2, reactome_v2, wp_v2)

# Continuation
pathways.filt <- pathways[pathways$gene %in% genes$entrezgene_id,]
colnames(pathways.filt) <- c("source", "target")

# Use "msigdbr"-package, even smarter
pathways_v2 <- Canonical_Pathways_humans  %>%
  dplyr::filter(gs_subcollection %in% c("CP:KEGG_MEDICUS", "CP:REACTOME", "CP:WIKIPATHWAYS")) %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  dplyr::rename(gene = ncbi_gene, term = gs_name) %>%
  dplyr::mutate(gene = as.integer(gene))

# Continuation
pathways.filt <- pathways[pathways$gene %in% genes$entrezgene_id,]
colnames(pathways.filt) <- c("source", "target")

# Being lazy! (For pathways)

if(!exists("pathways")){
  Databases_to_merge <- list("kegg_v2", "reactome_v2", "wp_v2")
  pathways <- do.call(rbind(Databases_to_merge))
  rm(kegg_v2, reactome_v2, wp_v2, Databases_to_merge)
  message("Combined Databases into 'pathways'.")
  # Filtration
  pathways.filt <- pathways[pathways$gene %in% genes$entrezgene_id, ]
  colnames(pathways.filt) <- c("source", "target")
  
  message("Filtered pathways into 'pathways.filt'.")
} else {
  message("'pathways' already exists.")
}

# Prep for Cytoscape

nodes.g <- genes[,c(5,2)]                 # This will select: "entrezgene_id & symbol"
colnames(nodes.g) <- c("id", "label")     # This will change: "entrezgene_id & symbol" into "id" & label"
nodes.p <- unique(pathways.filt[,c(1,1)]) # This will select: "source" & "target"
colnames(nodes.p) <- c("id", "label")     # This will change: "source" & "target" into "id" & label"
nodes <- rbind(nodes.g, nodes.p)


# Or
nodes_g_v2 <- genes %>%
  dplyr::select(entrezgene_id, symbol) %>%
  dplyr::rename(id = entrezgene_id, label = symbol)

nodes.p_v2 <- pathways.filt %>%
  dplyr::distinct(source, target) %>%
  dplyr::rename(id = source, label = target) %>% # to be fixed
  distinct()

nodes_v2 <- rbind(nodes_g_v2, nodes.p_v2)


# Or

nodes <- bind_rows(
  # Gene nodes
  genes %>%
    dplyr::select(entrezgene_id, symbol) %>%
    dplyr::rename(id = entrezgene_id, label = symbol),
  # Pathway nodes
  pathways.filt %>%
    dplyr::distinct(source, target) %>%
    dplyr::rename(id = source, label = target)
) %>%
  dplyr::distinct()

# Making the network
RCy3::createNetworkFromDataFrames(nodes=nodes,edges=pathways.filt, title="FUMA-all-memory-tests-merged", collection="HPA-FUMA")
RCy3::loadTableData(genes, data.key.column = "entrezgene_id", table.key.column = "id")

RCy3::createVisualStyle("pathway-gene")
setNodeColorMapping(table.column = "Memory",
                    table.column.values = c("PW", "LS", "PS"),
                    mapping.type = "d",
                    colors = c("#ff9933", "#6699cc", "#99cc99"), style.name = "pathway-gene")


fumajobs <- list.files(pattern = "^FUMA_")

# START SETTING UP THE SNP-GENE NETWORK
for (i in 1:length(fumajobs)) {
  snpgenenet <- total[[i]]
  
  # FILTER ON SNPs RELATED TO GENES PRESENT IN THE CANDIDATE LIST
  snpgenenet <- snpgenenet %>% filter(snpgenenet$gene %in% snpgenes[[i]]$ensembl_gene_id)
  
  # CREATE NODE TABLE FOR THE SNP-GENE NETWORK
  genenodes <- as.data.frame(unique(snpgenenet$gene))
  colnames(genenodes)[1] <- "nodes"
  snpnodes <- as.data.frame(unique(snpgenenet$rsID))
  colnames(snpnodes)[1] <- "nodes"
  node_table <- as.data.frame(rbind(genenodes, snpnodes))
  node_table <- na.omit(node_table)
  
  # CREATE EDGE TABLE FOR THE SNP-GENE NETWORK
  edge_table <- snpgenenet[,c("rsID", "gene")]
  edge_table <- as.data.frame(unique(edge_table))
  edge_table <- na.omit(edge_table)
  
  # CREATE THE SNP-GENE NETWORK
  nodes <- data.frame(id = node_table$nodes)
  edges <- data.frame(source = edge_table$rsID, target = edge_table$gene)
  createNetworkFromDataFrames(nodes, edges, title = paste0(fumajobs[i],"_snp_gene"), collection = "SNP-gene")
  rm(nodes, edges, node_table, edge_table)
  
  mapTableColumn(column = "id", species = "Human", 
                 map.from = "Ensembl", map.to = "Ensembl", force.single = T, network = paste0(fumajobs[i],"_snp_gene"))
  
  mapTableColumn(column = "id", species = "Human", 
                 map.from = "Ensembl", map.to = "Entrez Gene", force.single = T, network = paste0(fumajobs[i],"_snp_gene"))
  
  mapTableColumn(column = "id", species = "Human", 
                 map.from = "Ensembl", map.to = "HGNC", force.single = T, network = paste0(fumajobs[i],"_snp_gene"))
  
  # ADD TYPE DATA TO SNPs
  snpdata <- snpgenenet[,c("rsID", "CADD", "gwasP")]
  snpdata <- unique(snpdata)
  colnames(snpdata)[1] <- "name"
  snpdata$type <- "SNP in LD"
  snpdata$type[!is.na(snpdata$gwasP)] <- "SNP from GWAS"
  
  loadTableData(snpdata, data.key.column = "name", network = paste0(fumajobs[i],"_snp_gene"))
  rm(snpdata)
}

#
#
#
#
#
# MERGE BASED ON HGNC SYMBOL (MAP COLUMNS FIRST)
#
#
#
#
#


# ADD TYPE DATA TO GENES IN PATHWAYS
gene_pathnodes <- getTableColumns(table = "node", column = c("Ensembl","Entrez Gene","name"), network = paste0(fumajobs[i],"_merged"))
gene_pathnodes <- gene_pathnodes[!is.na(gene_pathnodes$'Ensembl'),]
gene_pathnodes$type[gene_pathnodes$'Entrez Gene' %in% databases$entrezgene] <- "Gene in Pathway"
gene_pathnodes$type[!gene_pathnodes$'Entrez Gene' %in% databases$entrezgene] <- "Gene"
gene_pathnodes <- gene_pathnodes[,c("name","type")]
loadTableData(gene_pathnodes, data.key.column = "name", network = paste0(fumajobs[i],"_merged"))
rm(gene_pathnodes)

# CREATE A UNIFIED COLUMN FOR NODE NAMES
snp <- getTableColumns(table = "node", column = c("name", "type"), network = paste0(fumajobs[i],"_merged"))
path <- getTableColumns(table = "node", column = c("name", "type"), network = paste0(fumajobs[i],"_merged"))

snp <- snp %>% filter(snp$type=="SNP in LD" | snp$type=="SNP from GWAS")
path <- path %>% filter(path$type=="Pathway")

gene <- getTableColumns(table = "node", column = c("name", "HGNC"), network = paste0(fumajobs[i],"_merged"))
gene <- gene %>% filter(!is.na(gene$'HGNC'))

colnames(gene)[2] <- "nodename"

snp$nodename <- snp$name
path$nodename <- path$name

# MAP THE NODE NAME DATA TO A "NODENAME" COLUMN AND APPLY NODE LABEL MAPPING
loadTableData(snp, data.key.column = "name", network = paste0(fumajobs[i],"_merged"))
loadTableData(path, data.key.column = "name", network = paste0(fumajobs[i],"_merged"))
loadTableData(gene, data.key.column = "name", network = paste0(fumajobs[i],"_merged"))
setNodeLabelMapping(table.column = "nodename", network = paste0(fumajobs[i],"_merged"))

rm(snp, path, gene)

gene_iwrd <- getTableColumns(table = "node", column = c("nodename","name","type"), network = paste0("FUMA_iwrd_merged"))
gene_list <- getTableColumns(table = "node", column = c("nodename","name","type"), network = paste0("FUMA_list_merged"))
gene_pics <- getTableColumns(table = "node", column = c("nodename","name","type"), network = paste0("FUMA_pics_merged"))

gene_iwrd$test <- "iwrd"
gene_list$test <- "list"
gene_pics$test <- "pics"

gene_iwrd <- gene_iwrd[gene_iwrd$type=="Gene" | gene_iwrd$type=="Gene in Pathway", ]
gene_list <- gene_list[gene_list$type=="Gene" | gene_list$type=="Gene in Pathway", ]
gene_pics <- gene_pics[gene_pics$type=="Gene" | gene_pics$type=="Gene in Pathway", ]

gene_test <- rbind(gene_iwrd,gene_list)
gene_test <- rbind(gene_test,gene_pics)

gene_test <- gene_test[,c(2,4)]

loadTableData(gene_test, data.key.column = "name", network = "Our_network")

# CREATE COLOR AND SHAPE VISUALIZATION BASED ON TYPE
# PATHWAYS: BLUE RECTANGLES
# GENES: ELLIPTICAL, BRIGHT GREEN IN PATHWAYS AND LIGHT GREEN FOR GENERAL
# SNPs: RED FOR ORIGINAL ARRAY-SNPs AND GREY FOR IMPUTED SNPs
setVisualStyle('default', network = paste0(fumajobs[i],"_merged"))
setNodeShapeMapping(table.column = "type",
                    table.column.values = c("SNP in LD", "SNP from GWAS", "Pathway", "Gene", "Gene in Pathway"),
                    shapes = c("Diamond", "Diamond", "Rectangle", "Ellipse", "Ellipse"),
                    network = paste0(fumajobs[i],"_merged"))

setNodeColorMapping(table.column = "type",
                    table.column.values = c("SNP in LD", "SNP from GWAS", "Pathway", "Gene", "Gene in Pathway"),
                    mapping.type = "d",
                    colors = c("#eeeeee", "#eeeeee", "#4daf4a", "#fdc086", "#ff7f00"),
                    network = paste0(fumajobs[i],"_merged"))

# CREATE NETWORK WITHOUT SNPs (FIGURE 2 IN PAPER)
# Introduce edges between genes with shared SNPs
gene_gene <- getTableColumns(table = "edge", column = c("source","target"), network = "TRY")
gene_gene1 <- gene_gene[grep("ENSG", gene_gene$target),]

library(igraph)
gene_gene1 <- graph_from_edgelist(as.matrix(gene_gene1))
gene_gene1 <- connect.neighborhood(gene_gene1, order = 2)
gene_gene1 <- as.data.frame(get.edgelist(gene_gene1))
gene_gene1 <- gene_gene1[grep("ENSG", gene_gene1$V1),]
gene_gene1 <- gene_gene1[grep("ENSG", gene_gene1$V2),]
colnames(gene_gene1) <- c("source", "target")

gene_names <- getTableColumns(table = "node", column = c("name","Ensembl", "Ensembl (1)", "Ensembl (2)", "type"), network = "TRY")
gene_names <- gene_names[!is.na(gene_names$Ensembl) | !is.na(gene_names$`Ensembl (1)`) | !is.na(gene_names$`Ensembl (2)`),]
gene_names$Ensembl[is.na(gene_names$Ensembl)] <- gene_names[is.na(gene_names$Ensembl),"Ensembl (1)"]
gene_names$Ensembl[is.na(gene_names$Ensembl)] <- gene_names[is.na(gene_names$Ensembl),"Ensembl (2)"]
gene_names <- gene_names[,c(1,2)]

for (i in 1:nrow(gene_gene1)) {
  name <- gene_gene1[i,"source"]
  shared <- gene_names[gene_names$Ensembl==name,"name"]
  gene_gene1[i,"source"] <- shared
}

for (i in 1:nrow(gene_gene1)) {
  name <- gene_gene1[i,"target"]
  shared <- gene_names[gene_names$Ensembl==name,"name"]
  gene_gene1[i,"target"] <- shared
}

gene_gene1 <- as.list(as.data.frame(t(gene_gene1)))
for (i in 1:length(gene_gene1)) {
  addCyEdges(gene_gene1[[i]], network = "TRY")
}

createColumnFilter("snp_filter", column = "type", criterion = "SNP", predicate = "CONTAINS", network = "TRY")
deleteSelectedNodes(network = "TRY")

setNodeSizeMapping(table.column = "type",
                   table.column.values = c("Pathway", "Gene", "Gene in Pathway"),
                   mapping.type = "d",
                   sizes = c(25, 45, 45),
                   network = "TRY")

