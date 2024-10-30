# Load Libraries

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(argparse)

# Argument parsing
parser <- ArgumentParser(description = "Read inputs")
parser$add_argument("-b", "--blen", type = "character", required = TRUE,
                    help = "Path to the medLen_trees.csv file")
parser$add_argument("-c", "--cpy", type = "character", required = TRUE,
                    help = "Path to the TEcopies.csv file")

args <- parser$parse_args()

# Load Data 
blen <- read.csv(args$blen, header = TRUE) 
cpy <- read.csv(args$cpy, header = TRUE)


# Process Data
cpy <- cpy %>% mutate(TEtype = gsub("/Merlin|/Helitron|/PiggyBac", "/Others", TEtype))
blen <- blen %>% mutate(TEtype = gsub("/Merlin|/Helitron|/PiggyBac", "/Others", TEtype))

blen <- merge(cpy, blen, by = c("assembly", "order","TEtype","TEname")) 

# Set Format
blen$assembly <- as.factor(blen$assembly)
blen$TEtype <- as.factor(blen$TEtype)
blen$numCopies <- as.numeric(blen$numCopies)

# Set Factor levels:

blen$order <- factor(blen$order,levels=c("DNA","LINE","PLE","LTR","SINE"))
blen$assembly <- factor(blen$assembly,levels=c("Fbus","Fgig_1","Fgig_2","Fhep_1","Fhep_2"))
blen$TEtype <- factor(blen$TEtype, 
                      levels = c("DNA/CMC", "DNA/Mutator-like", "DNA/Tc1-Mariner", 
                                 "DNA/Others", "DNA", "SINE", "PLE/Poseidon",
                                 "PLE", "LTR/Ty3-Gypsy", "LTR/Bel-Pao", "LTR", 
                                 "LINE/CR1", "LINE/CR1-Zenon", "LINE/Rex-Babar", 
                                 "LINE/RTE", "LINE"))

# Define color scale
scale.color <- scale_fill_manual(values = c("#E0F2DA", "#B0DEA2", "#5CB640", "#3D782A", "#203E16", 
                                            "#FFD700", "#DB99D0", "#8B317C", "#A4CDE6", "#3182BD", 
                                            "#04396d", "#FCDED4", "#FDAF9D", "#FB633F", "#A12103", 
                                            "#4E1002"), breaks = levels(blen$TEtype))



# Plot

gp <- blen %>% ggplot(aes(x=numCopies, y=medlen_tbranch,fill=TEtype)) + geom_point(shape=21,size=3) + scale.color + 
  scale_y_continuous(name="Median terminal branch length\n(substitutions per site)") + 
  scale_x_log10(name="TE copy-number",breaks = c(10^1,10^2,10^3,10^4),labels=trans_format("log10",math_format(10^.x))) + 
  theme_bw(base_size = 15)  + geom_hline(yintercept = 0.05, linetype="dashed",color = "red", linewidth=1.0) + facet_wrap(~assembly,ncol = 5) + 
  annotation_logticks() #  + theme(legend.position = "none")

ggsave("branch_lengths.tiff",plot = gp, dpi = 300,width =10,height = 5)
