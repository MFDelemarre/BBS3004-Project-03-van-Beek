# Template Script for Data Preparation FUMA 

# This script has a few different steps (and for some steps; 2 different option):
  # Step 0: Checking if all the packages are loaded (or install them)
  # Step 1: Loading the GWAS-dataset
          # 2 options
  # Step 2: Ensure the column-names are correct
  # Step 3: Sorting and filtering the data
          # 2 options
  # Step 4: Exporting it into a .gz


# Step 0: 
  # Loading and/or installing packages

# First check if some necessary packages are installed and loaded
Needed_Packages <- c("rstudioapi", "tidyverse")
install_Packs <- Needed_Packages %in% rownames(installed.packages())
if (any(!install_Packs)) {
  install.packages(Needed_Packages[!install_Packs])
}
lapply(Needed_Packages, library, character.only = TRUE)

# In this script I used the dataset: 
# https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010118/
# You can modify to what ever you want:
# https://fuma.ctglab.nl/tutorial#prepare-input-files
# You will need to change the headers; or at least remember the names; 
# Headers (mention in FUMA)
# Column names are automatically detected based on the following headers (case insensitive).
# SNP | snpid | markername | rsID: rsID
# CHR | chromosome | chrom: chromosome
# BP | pos | position: genomic position (hg19)
# A1 | effect_allele | allele1 | alleleB: affected allele
# A2 | non_effect_allele | allele2 | alleleA: another allele
# P | pvalue | p-value | p_value | pval: P-value (Mandatory)
# OR: Odds Ratio
# Beta | be: Beta
# SE: Standard error

# Step 1: 
# Loading the GWAS-dataset | (Option 2 - recommended)


# Option 1: Loading data from the same folder 
  # Make sure this script is in the same folder as your downloaded file
if(!exists("DATA.DIR")) {
  DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)} # This sets the "WORKING DIRECTORY TO CURRENT SCRIPT DIRECTORY"


# Change the file name (it's now -> "TEMPLATE_DATABASE_NAME.txt")
# Or change it to what
if(!exists("gwas_data") && file.exists("TEMPLATE_DATABASE_NAME.txt")) {
  message("Loading the GWAS Data... please wait (It will take a while.")
  gwas_data <- read_delim("TEMPLATE_DATABASE_NAME.txt")
  message("Finished loading")
  }

# Option 2: 
    # Choosing the file yourself
file_path <- file.choose()      # The chooses the file(and remembers the file-path)

if(!exists("gwas_data") && file.exists(file_path)) {
  message("Loading the GWAS Data... please wait (It will take a while.")
  gwas_data <- read_delim(file_path)
  message("Finished loading!")
}

glimpse(gwas_data)
# Step 2: Ensure the column-names are correct

# Now that you have your dataset read; check for the (column)names
# Check the column names; and look the name names for FUMA/SNP2GENE
oldcolnames <- colnames(gwas_data) 
print(oldcolnames)

# Check the column names, and change them to some shorter versions
# In my case the names were changed from left column to right.
    # Look the name names for FUMA/SNP2GENE (See line 31:39 for that)
# Names were changed:

# From these (old)              To these (new)
# MarkerName                  -> "SNP" [/ "rsID"] <- this was apparently also correct
# Chrpos: Chr:pos             -> "CHR_POS"
# Chr: Chromosome             -> "CHR"
# Pos: Position (hg19)        -> "BP"
# EA: Effect allele           -> "A1"
# NEA: Other allele           -> "A2"
# EAF: Effect allele frequency-> "EAF
# Neff: Effective sample size -> "N"
# Beta: Effect estimate       -> "BETA"
# SE: Standard error 
#     of effect estimate      -> "SE"
# P: P-value                  -> "P"

new_column_names <- c("SNP", "CHR_POS" ,"CHR", "BP", "A1", "A2", "EAF", "N", "BETA", "SE", "P")

# Check the new names of the columns first
colnames(gwas_data) <- new_column_names
print(colnames(gwas_data))

# Now you can filter; as you don't need all of the columns
if(!exists("gwas_filtered")) {
  gwas_filtered <- gwas_data %>% 
    select(SNP, CHR, BP, A1, A2, EAF, N, BETA, SE, P) # Filters by: ID, chromosome, location, effect allele and other allele, beta, effect frequency, p-value
}
# Now check the new columns and make sure they are complete
  # Or did not delete some of them.
print(colnames(gwas_filtered))
glimpse(gwas_filtered)


# Step 3
  # Sorting and filtering the data | (Option 2-recommended)
  # We are going for the top & bottom 5% of the Effect allele frequency (EAF)

# Option 1 (2-steps)
  # Filtering and checking

# GWAS-data arranged by P-value and top 5 percent of EAF is here ->
if(!exists("gwas_EAF_0.05")) {
  gwas_EAF_0.05 <- gwas_filtered %>% 
    arrange(P) %>% # Arrange P-value
    filter(EAF <= 0.95 & EAF >= 0.05) # Calculates MAF; and filters to  cut-off value '0.05'
}                                     


# Here you can check it is done correctly
min(gwas_EAF_0.05$EAF)
max(gwas_EAF_0.05$EAF)

if(!exists("gwas_EAF_0.05_step_2")) {
  gwas_EAF_0.05_step_2 <- gwas_EAF_0.05 %>% 
    select(-EAF)                       # This removes the "EAF"-column
}

# Option 2 (Doing it all at once)
  # Filters to cut-off value 0.05; &  Removes "Effect allele frequency"-column
  # In other words; it's basically the previous two things in one step
if(!exists("gwas_EAF_0.05_noEAF")) {
  gwas_MAF0.05_noEAF <- gwas_filtered %>% 
    arrange(P) %>%
    filter(EAF <= 0.95 & EAF >= 0.05) %>%   # Calculates MAF & uses cut-off value '0.05'
    select(-EAF)                            # removes "Effect allele frequency"
}                           

# Here you can see it's correct (for both steps)
glimpse(gwas_EAF_0.05_step_2)   # Option 1
glimpse(gwas_MAF0.05_noEAF)     # Option 2

# Making sure BP (position) is an integer (not a fraction)
# Otherwise it will NOT work
gwas_EAF_0.05_step_2$BP <- as.integer(gwas_EAF_0.05_step_2$BP)    # Option 1
gwas_MAF0.05_noEAF$BP <- as.integer(gwas_EAF_0.05_step_2$BP)      # Option 2


# if you want to double-double check; you can see here if the dataframes are the same
# IF YOU RENDERED BOTH!!
all.equal(gwas_MAF0.05_noEAF, gwas_EAF_0.05_step_2)


# Step 4
  # Exporting it into a .gz (so you can send it to FUMA)
  # Choose the function based on what file you have chosen

# If option 1 was used
write.table(
  gwas_EAF_0.05_step_2,
  file = gzfile("GWAS_for_FUMA.txt.gz", "w"),   # Compresses to .gz
  sep = "\t",        
  row.names = FALSE,
  quote = FALSE
)

# If option 2 was used 
write.table(
  gwas_MAF0.05_noEAF,
  file = gzfile("GWAS_for_FUMA.txt.gz", "w"),  # Compresses to .gz
  sep = "\t",        
  row.names = FALSE,
  quote = FALSE
)

# Steps after data preperation 
# You need to Check
  # Length of 'N' (you will need this)
length(unique(gwas_EAF_0.05_step_2N)) # Option 1
length(unique(gwas_MAF0.05_noEAF$N))  # Option 2
# You need to make sure you change the Reference panel population
