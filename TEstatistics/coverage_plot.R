# Load necessary libraries
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)

# Parse command line arguments
parse_args <- function(args) {
  input_path <- NULL
  output_path <- NULL
  format <- NULL
  
  for (i in seq_along(args)) {
    if (args[i] == "--input" && (i + 1) <= length(args)) {
      input_path <- args[i + 1]
    } else if (args[i] == "--output" && (i + 1) <= length(args)) {
      output_path <- args[i + 1]
    } else if (args[i] == "--format" && (i + 1) <= length(args)) {
      format <- args[i + 1]
    }
  }
  
  return(list(input = input_path, output = output_path, format = format))
}

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
parsed_args <- parse_args(args)

# Assign paths
input_path <- parsed_args$input
output_file <- parsed_args$output

# Check if paths are provided
if (is.null(input_path) || is.null(output_file)) {
  cat("Usage: Rscript plot_coverage.R --input <input_file> --output <output_path> [--format <pdf|tiff>]\n")
  cat("  --input:   Path to the input CSV file.\n")
  cat("  --output:  Path for saving the output plot (without file extension).\n")
  cat("  --format:  Output format for the plot (pdf or tiff). Default is pdf.\n")
  stop("Please provide valid arguments.\n")
}

# Set format
if (is.null(parsed_args$format)) {  
  format <- "pdf"  
} else {
  format <- parsed_args$format  
}

# Load summary table
TE_summ <- read.csv(file = input_path, header = TRUE)

# Rename minor TEtypes as "others"
TE_summ$TEtype <- gsub("Merlin", "Others", TE_summ$TEtype)
TE_summ$TEtype <- gsub("Helitron", "Others", TE_summ$TEtype)
TE_summ$TEtype <- gsub("PiggyBac", "Others", TE_summ$TEtype)

# Data type conversion
TE_summ$assembly <- as.factor(TE_summ$assembly)
TE_summ$TEtype <- as.factor(TE_summ$TEtype)
TE_summ$genic <- as.numeric(TE_summ$genic)
TE_summ$intergenic <- as.numeric(TE_summ$intergenic)

# Summarize total
TE_summ <- select(TE_summ, 1, 3, 5, 6)
TE_summ_2 <- TE_summ %>%
  group_by(assembly, TEtype) %>%
  summarise(genic = sum(genic), intergenic = sum(intergenic))

# Convert to larger format
TE_summ_2_lger <- TE_summ_2 %>%
  pivot_longer(!(assembly | TEtype), names_to = "Genomic_location", values_to = "size") %>%
  mutate(size = size / 1000000)  # Convert bp to Mb

# Set headers' names
names(TE_summ_2_lger) <- c("Assembly", "Repeat type", "Genomic location", "Size")

# Adjust values based on genomic location
TE_summ_2_lger <- TE_summ_2_lger %>%
  mutate(Size = ifelse(`Genomic location` == "intergenic", Size, -1 * Size))

# Set factor levels and define colors
TE_summ_2_lger$Assembly <- factor(TE_summ_2_lger$Assembly, levels = c("Fbus", "Fgig_2", "Fgig_1", "Fhep_2", "Fhep_1"))
TE_summ_2_lger$`Repeat type` <- factor(TE_summ_2_lger$`Repeat type`, levels = c("DNA/CMC", "DNA/Mutator-like", "DNA/Tc1-Mariner", "DNA/Others", "DNA", "SINE", "PLE/Poseidon", "PLE", "LTR/Ty3-Gypsy", "LTR/Bel-Pao", "LTR", "LINE/CR1", "LINE/CR1-Zenon", "LINE/Rex-Babar", "LINE/RTE", "LINE", "Satellite", "Unknown"))
scale.tmp <- scale_fill_manual(values = c("#E0F2DA", "#B0DEA2", "#5CB640", "#3D782A", "#203E16", "#FFD700", "#DB99D0", "#8B317C", "#A4CDE6", "#3182BD", "#04396d", "#FCDED4", "#FDAF9D", "#FB633F", "#A12103", "#4E1002", "#C0C0C0", "#646464"))

# Create the plot
plot <- TE_summ_2_lger %>%
  ggplot(aes(x = Assembly, y = Size, fill = `Repeat type`)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(name = "Genome coverage (Mb)", 
                     breaks = c(-200, -100, 0, 100, 200, 300, 400, 500, 600, 700, 800),
                     labels = abs(c(-200, -100, 0, 100, 200, 300, 400, 500, 600, 700, 800))) +
  theme_minimal(base_size = 20) +
  scale.tmp +
  geom_hline(yintercept = 0, linetype = 1, color = "black", linewidth = 1) + theme(legend.position = "none")

# Save the plot with specified dimensions
ggsave(paste0(output_file, ".", format), plot = plot, dpi = 300, width = 10, height = 5)
