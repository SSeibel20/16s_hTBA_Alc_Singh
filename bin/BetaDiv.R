# Beta Diversity
# Samantha Seibel
# Singh TBA Alcohol Mouse Study
# November 27th, 2024

# ---- beta diversity  ----
library(vegan)
library(pairwiseAdonis)
library(microViz)
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# check wd
getwd()

# set wd if needed
setwd("/Users/sls6550/work/hTBA_Alc_Singh")

# check wd to confirm change
getwd()

# SET PATHS
RDATA = "R"

#for reproducibility
set.seed(123245)

#load ps-filt
ps_decon <- readRDS("R/ps_decontam.Rds")

ntaxa(ps_decon) #788

#exploratory
ps_decon %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  dist_calc("aitchison") %>%
  dist_permanova(
    variables = c("Intervention", "Background"),
    perm = 9999
  ) #Intervention p = 0.001, Background p = 0.017

# clr transform phyloseq objects at Genus level
beta <- ps_decon %>%
  # clr transform to normalize
  tax_fix() %>%
  tax_transform(trans = "clr", rank = "Genus")

ps_get(beta) # 103 taxa

# generate distance matrix
psdist <- phyloseq::distance(beta, method = "euclidean")

#ADONIS test
# Intervention
adonis2(psdist ~ sample_data(beta)$Intervention, permutations = 10000)
#R^2 = 0.07533, p < 0.00001

# Background
adonis2(psdist ~ sample_data(beta)$Background, permutations = 10000)
#R^2 = 0.07391, p = 0.0381 (this is up for sig debate, depends on perm)

#pairwise for Sample_type
statdf <- pairwise.adonis(psdist, phyloseq::sample_data(beta)$Sample_Type)
statdf

#save table
write.table(statdf, file = file.path("out/pairwise-adonis.txt"), sep = "\t")

earth_pal <- c("#582F0E", "#A68A64", "#A4AC86", "#333D29")

#beta diversity
beta %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>%
  ord_plot(color = "Sample_Type", shape = "Intervention", size = 10, alpha=0.5) +
  scale_color_manual(values = earth_pal) +
  stat_ellipse(aes(colour = Sample_Type))+
  theme(axis.title= element_text(size = 20, color = "black"),
        axis.text = element_text(size = 20, color = "black"),
        legend.position = "right")

#save file of plot
ggsave("viz/betadiv_bytype_genus.pdf", plot = last_plot())
