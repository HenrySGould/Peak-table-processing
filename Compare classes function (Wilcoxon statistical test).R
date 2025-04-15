compare_classes <- function(class1_df, class2_df, label1, label2) {
  # Initialize a list to store results
  results_list <- list()
  
  for (i in 1:nrow(class1_df)) {
    # Extract values for each class
    values_class1 <- as.numeric(class1_df[i, ])
    values_class2 <- as.numeric(class2_df[i, ])
    
    # Perform one-sided Wilcoxon tests
    test_up <- wilcox.test(values_class1, values_class2, alternative = "greater")
    test_down <- wilcox.test(values_class1, values_class2, alternative = "less")
    
    # Calculate fold change
    fold_change <- mean(values_class1) / mean(values_class2)
    
    # Store results in the list with separate p-values for up and down regulation
    results_list[[i]] <- c(
      p_value_up = test_up$p.value,
      p_value_down = test_down$p.value,
      fold_change = fold_change
    )
  }
  
  # Convert results list to a data frame
  results_df <- as.data.frame(do.call(rbind, results_list))
  metabolitenames <- red_peaktable[-c(1:2), 1]
  rownames(results_df) <- metabolitenames
  
  # Adjust p-values for multiple testing (FDR correction)
  results_df$adjusted_p_value_up <- p.adjust(results_df$p_value_up, method = "fdr")
  results_df$adjusted_p_value_down <- p.adjust(results_df$p_value_down, method = "fdr")
  
  # Add a comparison label
  results_df$comparison <- paste(label1, "vs", label2)
  
  # Filter for significant metabolites (up or down regulated)
  significant_up <- results_df[results_df$adjusted_p_value_up < 0.05 & results_df$fold_change > 1, ]
  significant_down <- results_df[results_df$adjusted_p_value_down < 0.05 & results_df$fold_change < 1, ]
  
  # Combine significant upregulated and downregulated metabolites
  significant_metabolites <- rbind(significant_up, significant_down)
  
  # Sort by adjusted p-value for clarity
  significant_metabolites <- significant_metabolites[order(significant_metabolites$adjusted_p_value_up, 
                                                           significant_metabolites$adjusted_p_value_down), ]
  
  return(significant_metabolites)
}


