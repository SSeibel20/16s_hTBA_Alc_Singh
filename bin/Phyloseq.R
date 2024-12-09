# Building a Phyloseq Object 
# Samantha Seibel
# Singh TBA Alcohol Mouse Study
# November 26th, 2024

# ---- building ps object ----
library(phyloseq)
library(tidyverse)
library(tidylog)
library(dada2)
library(BiocManager)

#check wd
getwd()

#set wd if needed
setwd("/Users/sls6550/work/hTBA_Alc_Singh")

#check wd to confirm change
getwd()

# SET PATHS
RDATA = "R"

#need ASV table, Taxonomy, and Sample Metadata
#metadata first
meta <- read.csv("Metadata_16S_hTBA.csv")
head(meta)

#ensure no NAs
meta[is.na(meta)] <- "None"

#make SampleID names the row names to match with ASV
meta <- column_to_rownames(meta, var = "SampleID")

# LOAD ASV table from Dada2 script
ASV <- readRDS("R/asv_table.rds") 

# load tax table from Dada2 script # called taxa
taxa <-readRDS("R/tax_table.rds") 

# matrices in R
dim(ASV)
dim(taxa)

#construct into useable ps value
asv <- otu_table(ASV, taxa_are_rows=FALSE)
taxa <- tax_table(taxa)
meta <- sample_data(meta)


#create phyloseq object
ps <- phyloseq(asv, meta, taxa)

# get number of taxa
ntaxa(ps)

#get taxa ranks
rank_names(ps)

# access the data "slots" with @
head(ps@tax_table)
head(ps@otu_table)

#fix ASV names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
DNAps <- merge_phyloseq(ps, dna)
taxa_names(DNAps) <- paste0("ASV", seq(ntaxa(DNAps)))
DNAps

# access the data "slots" with @
head(DNAps@tax_table)
head(DNAps@otu_table)

# get taxa names
head(taxa_names(DNAps))
head(sample_names(DNAps))

# save as Rimage
saveRDS(DNAps, file = file.path(RDATA, "ps_raw.rds"))

#view ps
DNAps

# save image
save.image(file = file.path(RDATA, "building_ps.RData"))
