# ============================================ #
#                   Step 0:                    #
#         Installing & Loading packages        #
# ============================================ #
Packages_needed <- c("tidyverse", 
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

# Reading the files
annov <- read.table(paste0("annov.txt"), header = T, sep = "\t")            # Opens the file
eqtl <- read.table(paste0("eqtl.txt"), header = T, sep = "\t")              # Opens the file


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