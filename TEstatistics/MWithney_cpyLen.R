#Mann Withney U Test

#THE INPUT FILE ("inserts.csv") SHOULD LOOKS LIKE THIS:

#     assembly  order   TEtype    TEname          insLen  consLen
#     Fhep_1    SINE    SINE      bus-Botija2#    307     272
#     Fhep_1    SINE    SINE      bus-Botija2#    304     272
#     Fhep_1    SINE    SINE      bus-Botija2#    304     272
#     Fhep_1    SINE    SINE      bus-Botija2#    303     272
#     Fhep_1    SINE    SINE      bus-Botija2#    302     272

#THE ASSEMBLIES ARE INDICATED WITH THE NAME OF THE SPECIES FOLLOWED BY "_1" AND "_2" CORRESPONDINGLY.
#One-tailed test are made comparing if insertions in Fhep_2 are greater than the ones in Fhep_1 (Idem for F. gigantica)

# Load Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(argparse)

# Argument parsing
parser <- ArgumentParser()
parser$add_argument("--input", type = "character", help = "Path to the input file insertions.csv", required = TRUE)

args <- parser$parse_args()

# Load data
inserts <- read.csv(args$input, header = TRUE)

# Format input
inserts$TEname <- as.factor(inserts$TEname)
inserts$TEtype <- as.factor(inserts$TEtype)
inserts$assembly <- as.factor(inserts$assembly)
inserts$order <- as.factor(inserts$order)

# Function

make_MWUtest <- function(data, species) {
  
  # Select species
  TE_inserts <- subset(data,grepl(species,assembly))
  species <- c(species)
  
  # Vector of Family names
  namesTab <- unique(select(TE_inserts,TEname))
  namesVec <- as.character(namesTab[['TEname']])
  
  # Create an Empty Data Frame
  insMWUt <- data.frame(nameFam = character(0), pval =numeric(0))
  
  # Mann Withney U Test (one-tailed)
  for (i in 1:length(namesVec)) { 
      x1 <- namesVec[i]
      assembly_1 <- subset(TE_inserts, grepl(sprintf("%s_1", species), assembly) & grepl(x1, TEname))
      assembly_2 <- subset(TE_inserts, grepl(sprintf("%s_2", species), assembly) & grepl(x1, TEname))
      
      if (nrow(assembly_1) > 0 && nrow(assembly_2) > 0) {
        wilc <- wilcox.test(x = assembly_2$insLen, y = assembly_1$insLen, alternative = "greater", paired = FALSE)
        pval <- as.numeric(wilc[['p.value']])
        insMWUt <- rbind(insMWUt, data.frame(nameFam = x1, pval = pval))
      } else {
        warning(paste("Not enough data for TE family:", x1))
      }
  }

  insMWUt$pval <- as.numeric(insMWUt$pval)
  pval_adj <-p.adjust(insMWUt$pval, method ="BH") # Adjust p-values for multiple testing (Benjamini-Hochberg)
  MWUtest  <- cbind(insMWUt,pval_adj)
  
  # Formatting the Output
  compared_assemblies <- sprintf("%s_2-vs-%s_1", species, species)
  MWUtest <- cbind(assemblies = rep(compared_assemblies, nrow(MWUtest)), MWUtest)
  classification <- select(inserts,2,3,4) %>% unique() 
  # Merge and rearrange columns
  MWUTest <- MWUtest  %>%
    inner_join(classification, by = c("nameFam" = "TEname")) %>%
    select(assemblies,order, TEtype, nameFam, pval,pval_adj)
  
  return(MWUtest)
}


##############
#    MAIN    #
##############


# 1. Calculate Median Length of insertions (copies)
medianLen <- inserts %>% 
    group_by(assembly,order,TEtype,TEname) %>% 
    summarise(med_cpyLen = median(insLen),
              .groups = 'drop') 

write.csv(medianLen,"median_cpyLen.csv",row.names = FALSE)

# 2. Make MWUTest and plot the fraction of TE families with significantly increased copy-length in the improved assembly for each species.

# MWUTest
hep <- make_MWUtest(inserts,"Fhep")
gig <- make_MWUtest(inserts,"Fgig")
df <- rbind(hep, gig)
rm(hep,gig)

write.csv(df,"MannWithney_cpyLen.csv",row.names = FALSE)


# 3. PIE CHARTS

# Create a new column for significance
df <- df %>%
  mutate(significant = ifelse(pval_adj < 0.01, "Significant", "Not Significant"))  # Adjusted p-value: 0.01

# Calculate counts
count_df <- df %>%
  group_by(assemblies,significant) %>%
  summarise(count = n(),
           .groups = 'drop')

# Create a pie chart

count_df$assemblies <- gsub("-", " ", count_df$assemblies)
gp <- ggplot(count_df, aes(x = "", y = count, fill = significant)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(fill = "Result") +
  scale_fill_manual(values = c("Significant" = "#4D4D4D", "Not Significant" = "#BFBFBF")) +  
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  facet_wrap(~assemblies) + 
  geom_text(aes(label = count), 
            position = position_stack(vjus = 0.5), 
            color = "white", 
            size = 8)  

ggsave("pie_chart.tiff", plot = gp, dpi = 300, width = 10, height = 5)
