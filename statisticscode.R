install.packages("factoextra")
library(ggplot2)
library(factoextra)

#Define function to make PCA plots for each peaktable
process_files <- function(directory, pattern, classes_to_omit){

  files <- list.files(path = directory, pattern = pattern, full.names = TRUE)
  
  #Make for loop to do this for every peaktable in directory
  
  for(file in files){
  
    peaktable1 <- read.table(file, sep = "\t", header = TRUE, quote = "", check.names = FALSE, comment.char = "")
    
    sample_info <- grep("Area:", colnames(peaktable1))

    sample_names <- as.character(unlist(peaktable1[1, sample_info]))  # Remove the first column if it is not part of the data
    class_info <- as.character(unlist(peaktable1[2, sample_info]))

    # Remove the first two rows to get the data matrix
    data_matrix <- peaktable1[-c(1, 2), sample_info]
    data_matrix <- t(data_matrix)  # Transpose the matrix
    colnames(data_matrix) <- as.character(unlist(peaktable1[-c(1, 2), 1]))  # Set column names from the first column of the original data
   
    # Convert data to numeric
    data_matrix <- apply(data_matrix, 2, as.numeric)
    data_matrix <- as.data.frame(data_matrix)

    # Perform PCA
    pca_result <- prcomp(data_matrix, scale. = TRUE)

    # Extract PCA scores for the first two principal components
    pca_scores <- as.data.frame(pca_result$x[, 1:2])
    pca_scores$Sample <- sample_names
    pca_scores$Class <- class_info

    #Omit certain classes if necessary
    classes_to_omit <- c("QC", "M3-5", "ERED")
    pca_scores_filtered <- pca_scores[!pca_scores$Class %in% classes_to_omit, ]

    # Plot PCA results
    pca_plot <- ggplot(pca_scores_filtered, aes(x = PC1, y = PC2, color = Class)) +
      geom_point(size = 3) +
      labs(title = "PCA of Metabolomic Data",
           x = "Principal Component 1",
           y = "Principal Component 2") +
      theme_minimal() +
      theme(legend.title = element_blank())
    
    #Save plot in the directory
    ggsave(filename = paste0("PCA_plot_", basename(file), ".png"), plot = pca_plot, path = "C:/Data/HG Results/PCA plots", width = 8, height = 6)
    
    cat(file, "created in directory")  
  }
  
}

directory <- "C:/Data/HG Results"  # Define directory
pattern <- "reduced" # pattern for file names
classes_to_omit <- c("ERED", "QC", "M3-5")  # Replace Class1 and Class2 with the classes you want to omit

process_files(directory, pattern, classes_to_omit)


