doPCA <- function(peaktable, annot, x){
  
    sample_info <- grep("Area:", colnames(peaktable))
  
    sample_names <- as.character(unlist(colnames(peaktable[, sample_info])))  # Set names for matrix
    class_info <- as.character(unlist(annot[1, sample_names]))
  
    # Remove the first two rows to get the data matrix
    data_matrix <- peaktable[-c(1, 2), sample_info]
    data_matrix <- t(data_matrix)  # Transpose the matrix
    colnames(data_matrix) <- as.character(unlist(peaktable[-c(1, 2), 1]))  # Set column names from the first column of the original data
  
    # Convert data to numeric
    data_matrix <- apply(data_matrix, 2, as.numeric)
    data_matrix <- as.data.frame(data_matrix)
  
    # Perform PCA
     pca_result <- prcomp(data_matrix, scale. = FALSE) #check whether scale funciton is screwing with the PCA plots
  
    # Extract PCA scores for the first two principal components
    pca_scores <- as.data.frame(pca_result$x[, 1:2])
    ci <- as.factor(class_info)
    #Omit certain classes if necessary
    #classes_to_omit <- c("QC", "M3-5", "ERED")
    #pca_scores_filtered <- pca_scores[!pca_scores$Class %in% classes_to_omit, ]
    
    # Explain variance, Plot PCA results, comparison of data as boxplot
    var.e <- round(summary(pca_result)$importance[2, 1:2]*100, 1)
    png(paste0("Step", x, ".png"), width = 1024, height = 800)
    op <- par(mfrow = c(2, 1))
    boxplot(log2(t(data_matrix)))
    plot(pca_scores, col = as.numeric(ci), ylab = paste0("PC2 ", var.e[2], "%"), xlab = paste0("PC1 ", var.e[1], "%"), pch = c(5, 15, 19)[1+(as.numeric(ci)%/%6)])
    legend("topleft", levels(ci), col = as.numeric(as.factor(levels(ci))), pch = c(5, 15, 19)[1+(as.numeric(as.factor(levels(ci)))%/%6)], ncol = 3)
    dev.off()
    
}
