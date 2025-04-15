#Set directory 
setwd("C:/Data/HG results (SMILES)")
library(gplots)

#source wilcox comparison function
source("wilcox_comparison_function.R")

#read in peaktables
red_peaktable <- peaktable
class_info <- as.character(annot[1, ])
classes <- unique(class_info)

#create vector to store Results

subset_tables <- list()

# Subset the data into tables for each class
for (class_name in classes) {
  # Find the column indices that match the current class name
  sample_col <- which(class_info == class_name)
  
  # Subset the data frame based on the class
  subset_data <- red_peaktable[, sample_col, drop = FALSE]
  
  # Add the subset data to the list, using class_name as the key
  subset_tables[[as.character(class_name)]] <- subset_data
}

# Carry out significance test for each class vs pET
omit_group <- c("pET28a", "NA", "NorQC", "M3-5", "QC")
control_group <- "pET28a"

wilcox_results <- list()

for(class_name in names(subset_tables)){
  #skip pET28a group as comparison is against that
  if(is.na(class_name) || class_name %in% omit_group){
    next
  }
  
  #Retrive subset data
  subset_data <- as.data.frame(lapply(subset_tables[[class_name]], as.numeric))
  
  #Set comparison function
  control_data <- as.data.frame(lapply(subset_tables[[control_group]], as.numeric))
  
  # call significant test for function
  results_data <- compare_classes(subset_data, control_data, "table", "pET28a")
  
  wilcox_results[[class_name]] <- results_data
  
}


#Create list of metabolites which are varying significantly

# Initialize an empty vector to store metabolite names
all_metabolites <- c()

# Loop through wilcox_results to collect metabolite names
for (class_name in names(wilcox_results)) {
  if (!is.null(wilcox_results[[class_name]])) {
    all_metabolites <- c(all_metabolites, rownames(wilcox_results[[class_name]]))
  }
}

# Get unique metabolite names
unique_metabolites <- unique(all_metabolites)

# Print or save the unique metabolite list
print(unique_metabolites)

# Optionally, save to a file
write.table(unique_metabolites, "unique_metabolites.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

