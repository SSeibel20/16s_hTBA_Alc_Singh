library(ggplot2)
library(tidyverse)
library(viridis)
library(phyloseq)

# check wd
getwd()

# set wd if needed
setwd("/Users/sls6550/work/hTBA_Alc_Singh")

# check wd to confirm change
getwd()

# SET PATHS
RDATA = "R"

# load psdecontam  from ps script
ps_decon_rel <- readRDS("R/ps_decon_rel.rds")

ntaxa(ps_decon_rel) # 616

# filter out negative control in metadata (sample_dat)
ps_relab <- ps_decon_rel %>%
  ps_filter(!str_detect(sample_names(ps_decon_rel), "Buff"))

ntaxa(ps_relab) # 614 taxa

saveRDS(ps_relab, file = file.path(RDATA, "ps_rel_no_contr.rds"))

# fix and transform tax
ps_relab <- ps_relab %>%
  tax_fix()

# check if we have NAs
anyNA(tax_table(ps_relab)[,"Genus"])
anyNA(tax_table(ps_relab)[,"Family"])
anyNA(tax_table(ps_relab)[,"Phylum"])

# melt ps obj to make it useable in ggplot
ps_relab_melt <- ps_relab %>% psmelt()

# filter based on a threshold
abundance_threshold <- 0.05  
# 5% threshold (all taxa with 5% or greater abundance included)

ps_relab_filt_gen <- ps_relab_melt %>%
  filter(Abundance >= abundance_threshold) # 84 ASVs

fire_pal <- c("#03071E","#370617", "#6A040F", "#9D0208", "#D00000", "#DC2F02",
         "#E85D04", "#F48C06", "#FAA307", "#FFBA08", "#ffd97d")

earth_pal <- c("#03071E","#582F0E", "#7F4F24", "#936639", "#A68A64", "#B6AD90",
              "#C2C5AA", "#A4AC86", "#656D4A", "#414833", "#333D29")
         

# relative abundance with ggplot2
ggplot(ps_relab_filt_gen, aes(Sample, Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample", y = "Abundance") +
  ylim(0, 1) +
  scale_fill_manual(values = earth_pal) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = guide_legend(ncol = 4)) + # Set number of columns in legend
  theme(axis.text= element_text(size = 20),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

ggsave("viz/relab_genus.pdf", plot = last_plot(), height = 8, width = 26, units = "in")

# relative abundance by Family

# filter based on a threshold
abundance_threshold <- 0.01 
# 1% threshold (all taxa with 1% or greater abundance included)

ps_relab_filt_fam <- ps_relab_melt %>%
  filter(Abundance >= abundance_threshold)

# palettes
earth <- c("#03071E","#582F0E", "#7F4F24", "#936639", "#A68A64", "#B6AD90",
               "#C2C5AA", "#A4AC86", "#656D4A", "#414833", "#333D29", "#31572c")

#relative abundance with ggplot2
ggplot(ps_relab_filt_fam, aes(Sample, Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample", y = "Abundance") +
  ylim(0, 1) +
  scale_fill_manual(values = earth) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = guide_legend(ncol = 4)) + # Set number of columns in legend
  theme(axis.text= element_text(size = 20),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

ggsave("viz/relab_family.pdf", plot = last_plot(), height = 8, width = 26, units = "in")

#relative abundance by Phylum

# Aggregate data by taxa (e.g. Genus) and filter based on a threshold
abundance_threshold <- 0.01
# 1% threshold (all taxa with 1% or greater abundance included)

ps_relab_filt_fam <- ps_relab_melt %>%
  filter(Abundance >= abundance_threshold)

# palettes
earth <- c("#582F0E", "#936639", "#A68A64",
           "#C2C5AA", "#656D4A", "#333D29")

# relative abundance with ggplot2
ggplot(ps_relab_filt_fam, aes(Sample, Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample", y = "Abundance") +
  ylim(0, 1) +
  scale_fill_manual(values = earth) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = guide_legend(ncol = 4)) + # Set number of columns in legend
  theme(axis.text= element_text(size = 20),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

ggsave("viz/relab_phylum.pdf", plot = last_plot(), height = 8, width = 26, units = "in")



