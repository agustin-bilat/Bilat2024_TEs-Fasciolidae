library(ggplot2)
library(waffle)
library(dplyr)
library(tidyr)
library(argparse)

# Argument parsing
parser <- ArgumentParser(description = "Read inputs")
parser$add_argument("-a", "--aut", type = "character", required = TRUE,
                    help = "Path to the autonomous.tsv file")

args <- parser$parse_args()

# Load Data 
data <- read.table(args$aut, header = TRUE) 


# Transform data

data_lger <- data %>% 
  pivot_longer(!c(species, order, Repeat_type), names_to = "Functional classification", values_to = "count")

# Convert columns to factors and set levels
data_lger$species <- as.factor(data_lger$species)
data_lger$order <- factor(data_lger$order, levels = c("DNA", "LINE", "PLE", "LTR", "small_TEs", "Miscellaneous_Repeats"))
data_lger$Repeat_type <- factor(data_lger$Repeat_type, levels = c("DNA/CMC", "DNA/Mutator-like", "DNA/Tc1-Mariner", "DNA/Others", "DNA", "SINE", "PLE/Poseidon", "PLE", "LTR/Ty3-Gypsy", "LTR/Bel-Pao", "LTR", "LINE/CR1", "LINE/CR1-Zenon", "LINE/Rex-Babar", "LINE/RTE", "LINE", "Satellite", "Unknown"))
data_lger$`Functional classification` <- as.factor(data_lger$`Functional classification`)


# Palette of colors
pal_fill <- c("DNA/CMC" = "#E0F2DA", "DNA/Mutator-like" = "#B0DEA2", "DNA/Tc1-Mariner" = "#5CB640", 
              "DNA/Others" = "#3D782A", "DNA" = "#203E16", "SINE" = "#FFD700", 
              "PLE/Poseidon" = "#DB99D0", "PLE" = "#8B317C", "LTR/Ty3-Gypsy" = "#A4CDE6", 
              "LTR/Bel-Pao" = "#3182BD", "LTR" = "#04396d", "LINE/CR1" = "#FCDED4", 
              "LINE/CR1-Zenon" = "#FDAF9D", "LINE/Rex-Babar" = "#FB633F", "LINE/RTE" = "#A12103", 
              "LINE" = "#4E1002", "Satellite" = "#C0C0C0", "Unknown" = "#646464", 'z' = alpha('white', 0))

#Filter data
summary_data <- data_lger %>% 
  filter(count > 0) %>% 
  arrange(species, Repeat_type)



# Create the waffle plot
gp <- ggplot(summary_data, aes(fill = Repeat_type, values = count)) +
  waffle::geom_waffle(n_rows = 10, size = 0.2, color = "black", make_proportional = FALSE) +
  facet_grid(`Functional classification` ~ species) +
  scale_fill_manual(values = pal_fill) +  # Use my color palette
  coord_equal()+
theme(
  panel.grid = element_blank(),          # Remove grid lines
  axis.title = element_blank(),          # Remove axis titles
  axis.text = element_blank(),           # Remove axis text
  axis.ticks = element_blank(),          # Remove axis ticks
  #legend.position = "none",
  #strip.text= element_blank(),
  
  panel.spacing = unit(0.5, "lines")   # Adjust facet size
)

ggsave("waffle.tiff",plot = gp, dpi = 300,width =8,height = 6)
