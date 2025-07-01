library(ggplot2)
library(tidyverse)
library(microViz)
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
ps_decon <- readRDS("R/ps_decontam.rds")

ntaxa(ps_decon) # 788

# fix and transform tax
ps_relab <- ps_decon %>%
  tax_fix() %>%
  tax_transform("compositional")

# check if we have NAs
anyNA(tax_table(ps_relab)[,"Genus"])
anyNA(tax_table(ps_relab)[,"Family"])
anyNA(tax_table(ps_relab)[,"Phylum"])
  
# melt ps obj to make it useable in ggplot
ps_relab_melt <- ps_relab %>% psmelt()

# reorder samples by groups
reorder <- c("CT-1", "CT-2", "CT-3", "CT-4", "CT-13", "CT-14", "CT-5", "CT-6", "CT-7", "CT-15", "CT-16", "CT-17", "HT-16", "HT-18", "HT-28", "HT-32", "HT-40", "HT-41", "HT-19", "HT-20", "HT-21", "HT-8", "HT-30", "HT-37")

ps_relab_melt <- ps_relab_melt %>%
  mutate(
    Sample = factor(Sample, levels = reorder)
  )

# 5% threshold (all taxa with 1% or greater abundance included for descriptive)
abundance_threshold = 0.01
ps_relab_high <- ps_relab_melt %>%
  filter(Abundance >= abundance_threshold) # 496 ASVs

ps_relab_gen_low <- ps_relab_melt %>%
  filter(Abundance < abundance_threshold) # 18,912 ASVs

# Replace these low-abundance Genera with "Other"
ps_relab_melt_grouped <- ps_relab_melt %>%
  mutate(
    Genus_grouped = ifelse(Abundance < 0.01, "Other", Genus)
  )

# relative abundance with ggplot2 Genus
ggplot(ps_relab_melt_grouped, aes(Sample, Abundance, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(xlab = "Sample", ylab = "Abundance (Proportion)", title = "Relative Abundance by Genus", fill = "Genus", caption = "ASVs less than 1% abundance listed as Other") +
  ylim(0, 1) +
  scale_fill_viridis_d(option = "magma") +
  guides(fill = guide_legend(ncol = 8)) + # Set number of columns in legend
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 16),
        axis.text.y = element_text(size = 20), 
        legend.text = element_text(size = 10),
        legend.title = element_text(angle = 90, size = 20),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_text(angle = 360,  vjust = 0.5, hjust = 0.5, size = 24),
        plot.title = element_text(size = 24, hjust = 0.5),
       legend.position = "bottom")

ggsave("viz/relab_genus1.pdf", plot = last_plot(), height = 8, width = 16, units = "in")

# relative abundance by Family

# 1% threshold (all taxa with 1% or greater abundance included for descriptive)

# Replace these low-abundance Family with "Other"
ps_relab_melt_grouped_fam <- ps_relab_melt %>%
  mutate(
    Family_grouped = ifelse(Abundance < 0.01, "Other", Family)
  )

# relative abundance with ggplot2 Family
ggplot(ps_relab_melt_grouped_fam, aes(Sample, Abundance, fill = Family_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(xlab = "Sample", ylab = "Abundance (Proportion)", title = "Relative Abundance by Family", fill = "Family", caption = "ASVs less than 1% abundance listed as Other") +
  ylim(0, 1) +
  scale_fill_viridis_d(option = "magma") +
  guides(fill = guide_legend(ncol = 7)) + # Set number of columns in legend
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 16),
        axis.text.y = element_text(size = 20), 
        legend.text = element_text(size = 10),
        legend.title = element_text(angle = 90, size = 20),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_text(angle = 360,  vjust = 0.5, hjust = 0.5, size = 24),
        plot.title = element_text(size = 24, hjust = 0.5),
        legend.position = "bottom")

ggsave("viz/relab_family.pdf", plot = last_plot(), height = 8, width = 16, units = "in")

#relative abundance by Phylum

# Aggregate data by taxa (e.g. Genus) and filter based on a threshold
abundance_threshold <- 0.01
# 1% threshold (all taxa with 1% or greater abundance included)

ps_relab_filt_phy <- ps_relab_melt %>%
  filter(Abundance >= abundance_threshold)

# Replace these low-abundance Phylum with "Other"
ps_relab_melt_grouped_phy <- ps_relab_melt %>%
  mutate(
    Phylum_grouped = ifelse(Abundance < 0.01, "Other", Phylum)
  )

# relative abundance with ggplot2 Phylum
ggplot(ps_relab_melt_grouped_phy, aes(Sample, Abundance, fill = Phylum_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(xlab = "Sample", ylab = "Abundance (Proportion)", title = "Relative Abundance by Phylum", fill = "Phylum", caption = "ASVs less than 1% abundance listed as Other") +
  ylim(0, 1) +
  scale_fill_viridis_d(option = "magma") +
  guides(fill = guide_legend(ncol = 6)) + # Set number of columns in legend
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 16),
        axis.text.y = element_text(size = 20), 
        legend.text = element_text(size = 10),
        legend.title = element_text(angle = 90, size = 20),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_text(angle = 360,  vjust = 0.5, hjust = 0.5, size = 24),
        plot.title = element_text(size = 24, hjust = 0.5),
        legend.position = "bottom")

ggsave("viz/relab_phylum.pdf", plot = last_plot(), height = 8, width = 16, units = "in")
