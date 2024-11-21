# Load Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(argparse)
library(reshape2)
library(corrr)
library(purrr)

# Argument parsing
parser <- ArgumentParser(description = "Read inputs")
parser$add_argument("--cov", type = "character", required = TRUE, help = "Path to the coverage.csv file")
parser$add_argument("--cpy", type = "character", required = TRUE, help = "Path to the TEcopies.csv file")
parser$add_argument("--len", type = "character", required = TRUE, help = "Path to the median_cpyLen.csv file")
parser$add_argument("--age", type = "character", required = TRUE, help = "Path to the medLen_trees.csv file")
args <- parser$parse_args()

# Load Data 
cov <- read.csv(args$cov, header = TRUE) 
cpyNum <- read.csv(args$cpy, header = TRUE)
cpyLen <- read.csv(args$len,header=TRUE)
medLen_tbranch <- read.csv(args$age,header=TRUE)

# Format inputs

cov <- cov[!grepl("Unknown", cov$TEtype), ]
cov <- cov[!grepl("Satellite", cov$TEtype), ]
cov <- cov %>% mutate(TEtype = gsub("/Merlin|/Helitron|/PiggyBac", "/Others", TEtype))
cov <- cov %>%
  mutate(cov = genic + intergenic) %>% 
  select(1,2,3,4,7)
cpyLen <- cpyLen %>% mutate(TEtype = gsub("/Merlin|/Helitron|/PiggyBac", "/Others", TEtype))
medLen_tbranch <- medLen_tbranch %>% mutate(TEtype = gsub("/Merlin|/Helitron|/PiggyBac", "/Others", TEtype))

# Merge tables
TEfams.tab <- merge(cov, cpyNum, by = c("assembly", "order","TEtype","TEname"))
TEfams.tab <- merge(TEfams.tab, cpyLen, by = c("assembly", "order","TEtype","TEname"))
TEfams.tab <- merge(TEfams.tab, medLen_tbranch, by = c("assembly", "order","TEtype","TEname"))

# Difference function
difference <- function(A, B) {
  A - B
}

# Function to process data for specified assemblies
process_assembly_data <- function(assembly_names) {
  
  # Initialize an empty list to store results for each assembly
  results_list <- list()
  
  # Loop through each assembly name
  for (assembly_name in assembly_names) {
    
    # Filter the TEfams.tab dataset by the current species (Fhep or Fgig) 
    TEfams2.tab <- TEfams.tab[grepl(assembly_name, TEfams.tab$assembly), ]
    
    # Reshape and calculate the differences in the statistics associated to each TE family between the long-read versus the short-read assemblies in the current species:
    results <- TEfams2.tab %>%
      pivot_wider(names_from = assembly, values_from = c(cov, numCopies, numFullCopies, med_cpyLen, medlen_tbranch)) %>%
      mutate(
        cov_change = difference(get(paste0("cov_", assembly_name, "_2")), get(paste0("cov_", assembly_name, "_1"))),
        numCopies_change = difference(get(paste0("numCopies_", assembly_name, "_2")), get(paste0("numCopies_", assembly_name, "_1"))),
        numFullCopies_change = difference(get(paste0("numFullCopies_", assembly_name, "_2")), get(paste0("numFullCopies_", assembly_name, "_1"))),
        med_cpyLen_change = difference(get(paste0("med_cpyLen_", assembly_name, "_2")), get(paste0("med_cpyLen_", assembly_name, "_1"))),
        medlen_tbranch_change = difference(get(paste0("medlen_tbranch_", assembly_name, "_2")), get(paste0("medlen_tbranch_", assembly_name, "_1")))
      )
    
    # Store the result in the results list with assembly name as the key
    results_list[[assembly_name]] <- results
  }
  
  # Return the list of results for all assemblies
  return(results_list)
}

# Call the function for both Fgig and Fhep
assembly_results <- process_assembly_data(c("Fgig", "Fhep"))

# Access the results for each assembly from the list
results_Fgig <- assembly_results[["Fgig"]]
results_Fhep <- assembly_results[["Fhep"]]

# Create the df_gig and df_hep tables for correlation plotting
df_gig <- results_Fgig %>% select(cov_change, numCopies_change, numFullCopies_change, med_cpyLen_change, medlen_tbranch_change)
df_hep <- results_Fhep %>% select(cov_change, numCopies_change, numFullCopies_change, med_cpyLen_change, medlen_tbranch_change)

# Function to calculate correlations and plot the correlogram
plot_correlogram <- function(df_gig, df_hep) {
  
  # Calculate correlation matrices for both Fgig and Fhep
  cor_gig <- cor(df_gig, method = "spearman")
  cor_hep <- cor(df_hep, method = "spearman")
  
  # Melt both correlation matrices to long format
  cor_melted_gig <- melt(cor_gig)
  cor_melted_hep <- melt(cor_hep)
  
  # Rename columns
  colnames(cor_melted_gig) <- c("x", "y", "r")
  colnames(cor_melted_hep) <- c("x", "y", "r")
  
  # Calculate p-values for both datasets
  cor_melted_gig <- cor_melted_gig %>%
    mutate(p_value = map2_dbl(x, y, 
                              ~ cor.test(df_gig[[.x]], df_gig[[.y]], method = "spearman", exact = FALSE)$p.value))
  
  cor_melted_hep <- cor_melted_hep %>%
    mutate(p_value = map2_dbl(x, y, 
                              ~ cor.test(df_hep[[.x]], df_hep[[.y]], method = "spearman", exact = FALSE)$p.value))
  
  # Add significance indicator
  cor_filtered_gig <- cor_melted_gig %>%
    mutate(significance = ifelse(p_value < 0.01, "*", ""))
  
  cor_filtered_hep <- cor_melted_hep %>%
    mutate(significance = ifelse(p_value < 0.01, "*", ""))
  
  # Combine the data for both Fgig and Fhep for plotting
  cor_combined <- bind_rows(
    cor_filtered_gig %>% mutate(assembly = "Fgig"),
    cor_filtered_hep %>% mutate(assembly = "Fhep")
  )
  
  # Set factor levels for axis ordering
  cor_combined$x <- factor(cor_combined$x, levels = c("cov_change", "numCopies_change", "numFullCopies_change", "med_cpyLen_change", "medlen_tbranch_change"))
  cor_combined$y <- factor(cor_combined$y, levels = c("cov_change", "numCopies_change", "numFullCopies_change", "med_cpyLen_change", "medlen_tbranch_change"))
  
  # Create the plots of the correlograms
  ggplot(data = cor_combined, aes(x = x, y = y, fill = r)) +  
    geom_tile(color = "white") +  
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Spearman Correlation") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Correlogram") +
    geom_text(aes(label = paste0(round(r, 2), significance)), color = "black", size = 3) + 
    scale_x_discrete(labels = c(
      cov_change = expression(Delta ~ " (Coverage)"),               
      numCopies_change = expression(Delta ~ " (TE copy-number)"),          
      numFullCopies_change = expression(Delta ~ " (Full-length TE copy-number)"),
      med_cpyLen_change = expression(Delta ~ " (TE copy-length (median))"),  
      medlen_tbranch_change = expression(Delta ~ " (TE age)")  
    )) +
    scale_y_discrete(labels = c(
      cov_change = "", 
      numCopies_change = "",
      numFullCopies_change = "",
      med_cpyLen_change = "",
      medlen_tbranch_change = "       "
    )) +
    facet_wrap(~assembly, scales = "free")
}


# Plot the correlograms:
gp <- plot_correlogram(df_gig, df_hep)

ggsave("correlograma.tiff",plot = gp, dpi = 300,width =10,height = 5)
