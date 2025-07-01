# ---- Decontam ----
#based on negative controls
library(decontam)
library (microViz)
library(phyloseq)

#get working directory
getwd()

#set working directory
setwd("/Users/sls6550/work/hTBA_Alc_Singh")

#load ps-filt
psraw <- readRDS("R/ps_raw.rds")

ntaxa(psraw) # 794

get_taxa_unique(psraw, "Kingdom") # Bacteria only no need to filter others
get_taxa_unique(psraw, "Phylum") # 7 different phyla including NA

# get only bacteria
#psf <- subset_taxa(ps, Kingdom == "Bacteria")

# remove NA phylum
psraw <- subset_taxa(psraw, !is.na(Phylum) & !Phylum %in% c("", "NA"))

ntaxa(psraw) # 790

# Validate the phyloseq object and remove undetected taxa
psf <- psraw %>%
  # Assuming tax_fix is fixing taxonomy issues
  tax_fix() %>%
  # Validate the object and remove undetected taxa
  phyloseq_validate(remove_undetected = TRUE) %>%
  # Create a new column 'SampleBinary' based on sample names
  ps_mutate(SampleBinary = if_else(str_detect(sample_names(psraw), "Neg"),
                                true = "Control", false = "Sample"))

#check
sample_data(psf)

ntaxa(psf) # 790

# what do raw controls look like
neg <- psf %>%
  # get negatives
  ps_filter(str_detect(SampleBinary, "Control")) %>%
  tax_fix()

factor(sample_data(neg)$SampleBinary)

ntaxa(neg) # 2 taxa
view(tax_table(neg)) # Ruminococcus_1 and Roseburia

## inspect library sizes
sampdf <- data.frame(sample_data(psf))
sampdf$LibrarySize <- sample_sums(otu_table(psf))
sampdf <- sampdf[order(sampdf$LibrarySize), ]
sampdf$Index <- seq(nrow(sampdf))
ggplot(data = sampdf, aes(x = Index, y = LibrarySize, color = SampleBinary)) +
  geom_point() 

# the prevalence method of decontam is used
# we cannot use frequency due to lack of DNA quantification from experiment

sample_data(psf)$is.neg <- sample_data(psf)$SampleBinary == "Control"
contamdf.prev <- isContaminant(psf, method="prevalence", neg="is.neg", threshold=0.5) 
table(contamdf.prev$contaminant) # 790 taxa, no contaminants

#filter out negative control in metadata (sample_dat)
psdecon <- psf %>%
  ps_filter(!str_detect(sample_names(psf), "Buff")) 

ntaxa(psdecon) # 788

#save
saveRDS(psdecon, file = file.path("R/ps_decontam.Rds"))

# above RDS will be used for Aldex2 (differential relative abundance)

# Filter by relative abundance (less than 10e-5 not relevant in literature)
# using psf to have negative controls still present for decontam

#transform to relative abundance
#pst <- psdecon %>%
  #tax_transform("compositional") 

# below does the same thing as above, use whichever you prefer
#pstest <- transform_sample_counts(psdecon, function(x) x / sum(x) ) #788

# Check if the OTU tables are identical in values
# Assuming you have two phyloseq objects: ps1 and ps2
#otu1 <- otu_table(pst)
#otu2 <- otu_table(pstest)
#identical(otu1, otu2)

#ntaxa(pst) #788 taxa

#save
#saveRDS(pst, file = file.path("R/ps_decontam_rel.Rds"))

#remove taxa with total relative abundance less than 10e-5
#psr <- filter_taxa(pst, function(x) mean(x) > 1e-5, TRUE) #614 taxa (173 removed)

#ntaxa(psr) # 616

#save decontam ps as RImage
#saveRDS(psr, file = file.path("R/ps_decon_rel.rds"))

