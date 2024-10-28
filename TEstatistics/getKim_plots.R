# Load Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(argparse)

# Argument parsing
parser <- ArgumentParser()
parser$add_argument("--input", type = "character", help = "Path to the input file obtained after joining with 'cat' all
                    the outputs obtained with the bash script 'getKimura_from_align.sh' from the alignment (.aln) files, 
                    required = TRUE")

args <- parser$parse_args()

# Load data
kim.tab <- read.table(args$input, header = FALSE)

# Format data
names(kim.tab) <- c("assembly","TEname","TEtype","percent_of_genome","kimura")
kim.tab <- kim.tab %>% mutate(TEtype = gsub("/Merlin|/Helitron|/PiggyBac", "/Others", TEtype)) # Substitute TEtype-names
kim.tab <- kim.tab %>% mutate(TEtype = gsub("LTR/Unknown","LTR",TEtype))

# Set factors and color palette
kim.tab$TEtype <- factor(kim.tab$TEtype, levels=c("DNA/CMC","DNA/Mutator-like","DNA/Tc1-Mariner","DNA/Others","DNA","SINE","PLE/Poseidon","PLE","LTR/Ty3-Gypsy","LTR/Bel-Pao","LTR","LINE/CR1","LINE/CR1-Zenon","LINE/Rex-Babar","LINE/RTE","LINE","Unknown"))
scale.tmp <- scale_fill_manual(values=c("#E0F2DA","#B0DEA2","#5CB640","#3D782A","#203E16","#FFD700","#DB99D0","#8B317C","#A4CDE6","#3182BD","#04396d","#FCDED4","#FDAF9D","#FB633F","#A12103","#4E1002","#646464"))

# PLOTS
# 1. Stacked Kimura

gp <- ggplot(data=kim.tab,aes(x=kimura,y=percent_of_genome,fill=TEtype)) + geom_bar(stat="identity",color="black") + 
  scale_y_continuous(name="Percent of the Genome") + scale_x_continuous(name="Kimura substitution level (CpG adjusted)",limits=c(-1,50))  + 
  labs(fill="") + ggtitle("") + theme_minimal(base_size = 20) + facet_wrap(~assembly,ncol = 5) + geom_col(width= 0.8) + theme(legend.position="right") + 
  scale.tmp

ggsave('kim_stack.tiff', plot = gp, device = "tiff", dpi = 300,height = 10,width = 20)

# 2. TEtype kimura

#kim.tab <- subset(kim.tab,!grepl("Unknown",TEtype))

gp2 <- ggplot(kim.tab,aes(x=kimura,y=percent_of_genome,fill=TEtype)) + 
  geom_bar(stat="identity",color="black",show.legend = FALSE) + 
  scale_y_continuous(n.breaks=3,name="Percent of the Genome") + 
  scale_x_continuous(name="Kimura substitution level (CpG adjusted)",limits=c(-1,50)) + 
  facet_grid(TEtype~assembly,scale='free_y') + labs(fill="") + scale.tmp + geom_col(width= 0.8) + 
  theme_bw(base_size = 14) + theme(strip.text.y = element_blank(),legend.position = "right")  

ggsave('kim_tetype.tiff', plot=gp2, dpi = 300, height = 8,width = 8)
