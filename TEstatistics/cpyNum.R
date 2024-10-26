library(ggplot2)
library(dplyr)
library(scales)
library(argparse)

# Argument parsing
parser <- ArgumentParser()
parser$add_argument("--input", type = "character", help = "Path to the input file insertions.csv", required = TRUE)

args <- parser$parse_args()

# Load data
df <- read.csv(args$input, header = TRUE)

# Process Input Data

df2 <- df %>% 
  group_by(assembly,order,TEtype,TEname) %>% 
  summarise(numCopies = n(),  # TE Copy-number
            numFullCopies=sum(insLen > 0.85 * consLen & insLen < 1.15 * consLen), # Full-length Copy-number (85% to 115% of consensus length)
            .groups = 'drop') 

df3 <- df2 %>% mutate(TEtype = gsub("/Merlin|/Helitron|/PiggyBac", "/Others", TEtype)) # Substitute TEtype-names

# Set Factor levels:
df3$TEtype <- factor(df3$TEtype, 
                     levels = c("DNA/CMC", "DNA/Mutator-like", "DNA/Tc1-Mariner", 
                                "DNA/Others", "DNA", "SINE", "PLE/Poseidon",
                                "PLE", "LTR/Ty3-Gypsy", "LTR/Bel-Pao", "LTR", 
                                "LINE/CR1", "LINE/CR1-Zenon", "LINE/Rex-Babar", 
                                "LINE/RTE", "LINE"))

df3$order <- factor(df3$order,levels=c("DNA","LINE","PLE","LTR","SINE"))

# Define color scale
scale.color <- scale_fill_manual(values = c("#E0F2DA", "#B0DEA2", "#5CB640", "#3D782A", "#203E16", 
                                            "#FFD700", "#DB99D0", "#8B317C", "#A4CDE6", "#3182BD", 
                                            "#04396d", "#FCDED4", "#FDAF9D", "#FB633F", "#A12103", 
                                            "#4E1002"), breaks = levels(df3$TEtype))


# Functions

scientific <- function(x) {
  ifelse(x == 0, "0", parse(text = gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

make_plot <- function(data, y_col, y_label, file_name) {
  ggplot(data, aes(x = assembly, y = !!sym(y_col), fill = TEtype)) + 
    geom_bar(stat = "identity", width = 0.7) + 
    facet_wrap(~ order, ncol = 5) + 
    scale_y_continuous(name = y_label, labels = scientific) +
    scale_x_discrete(name = "Assembly name") +
    theme_linedraw(base_size = 16) + 
    theme(axis.ticks = element_line(colour = "black"),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          legend.position = "none") + 
    scale.color
  ggsave(file_name, dpi = 300, width = 10, height = 5)
}

# generate output files:

# 1/ Plots:
make_plot(df3, "numCopies", "TE copy-number", "cpyNum.tiff")

make_plot(df3, "numFullCopies", "Full-length TE copy-number", "full-CpyNum.tiff")

# 2/ Table
write.csv(df3,"TEcopies.csv",row.names = FALSE)
