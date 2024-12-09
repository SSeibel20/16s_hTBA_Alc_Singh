library(ALDEx2)
library(tidyverse)
library(microViz)
library(phyloseq)

# check wd
getwd()

# set wd if needed
setwd("/Users/sls6550/work/hTBA_Alc_Singh")

# check wd to confirm change
getwd()

# load ps from decontam
ps_decon <- readRDS("R/ps_decontam.Rds")

# Trying multiple variable statistics in Aldex2
otu_table_decontam <- as.data.frame(phyloseq::otu_table(ps_decon))
meta_decontam <- data.frame(phyloseq::sample_data(ps_decon))

# confirm class type
class(otu_table_decontam)
class(meta_decontam)


# check assumptions of ttest
# normality is assumed to be false due to count-based data, zero-inflation, skew, and non-independence 

# equal variance
# Convert OTU table into long format: one row per OTU sample pair
otu_lev <- otu_table_decontam %>%
  as.data.frame() %>%
  mutate(Sample_Type = meta_decontam$Sample_Type) %>%
  gather(key = "OTU", value = "Count", -Sample_Type)

# Run Levene's Test (Check variance across different groups)
library(car)
leveneTest(Count ~ Sample_Type, data = otu_lev) # equal variance

# need to make sure formating of tables is good for aldex.clr
all(colnames(otu_table_decontam) == rownames(meta_decontam))  
# Should return FALSE, check visually

dim(otu_table_decontam)  # Returns features, samples (24, 785)
dim(meta_decontam) # Returns samples, features (24, 7) 
# NEED samples to be the same

# Transpose the rows to columns
otu_table_decontam <- t(otu_table_decontam)

all(colnames(otu_table_decontam) == rownames(meta_decontam))  
# Should return TRUE, check visually

# test clr on Intervention to check
clr_int <- aldex.clr(otu_table_decontam, conds = meta_decontam$Intervention)

#---- Run sensitivity analysis function (Derived from SK) ----

##Plotting TP and FP over different levels of gamma
gamma.to.test <- c(0, 0.1, 0.3, 0.5, 1, 3, 5)
sen_res_int <- aldex.senAnalysis(clr_int, gamma = gamma.to.test)
head(sen_res_int[[1]])

##Plotting the sensitivity results.
plotGamma(sen_res_int, thresh = 0.5)
#not much deviation or clear patterns, implies results are not sens to gamma

#---- Diff Rel Ab Intervention ----
#aldex.clr
clr_int <- aldex.clr(otu_table_decontam, conds = meta_decontam$Intervention, 
                     gamma = 0.5, mc.samples = 1000, denom = "none")

clr_int@scaleSamps

#plot scale
graph.df <- data.frame("scale" = clr_int@scaleSamps[,1], 
                       group = meta_decontam$Intervention)

ggplot(graph.df, aes(x = group, y = scale))+
  geom_boxplot() +
  labs(title = paste("Differential Relative Abundance: Alcohol vs Control"), 
       x = "Intervention", y = "Scale")

#Running ttest for Intervention
x.ss.tt.int <- aldex.ttest(clr_int)

#Running effect
x.ss.e.int <- aldex.effect(clr_int, CI=T)

#combine for a plot
x.ss.all.int <- cbind(x.ss.e.int, x.ss.tt.int)

#save graph
tiff(filename = "viz/intervention_diff_relative_abundance.tiff", width = 1300, height = 1300, units = "px", pointsize = 10, res = 300)

# look at plot to see if positive/negative effect size corresponds to which size
int_dra <- aldex.plot(x.ss.all.int, main = "Differential Relative Abundance by Intervention", type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

#close tiff
dev.off()

#check if there is significant features
sig_aldex_int <- x.ss.all.int %>%
  filter(wi.eBH < 0.05)

View(sig_aldex_int)
#two significant features

#make a table of all correct p-values
all_aldex2_int <- x.ss.all.int %>%
  rownames_to_column(var = "OTU") %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(psf3))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

# add in previously formed taxa information to complete the table
all_aldex2_int <- left_join(all_aldex2_int, taxa_info)
save(all_aldex2_int, file = "R/table_all_taxa_int.RData")

#save in csv or xlsx tables so you can use the table for build a plot
write.table(all_aldex2_int, file = "out/table_all_taxa_int.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)

# make a table of significant corrected p-values
sig_aldex2_int2 <- x.ss.all.int %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# add in previously formed taxa information to complete the table
sig_aldex2_int2 <- left_join(sig_aldex2_int2, taxa_info)

#save in csv or xlsx tables so you can use the table for build a plot
write.table(sig_aldex2_int2, file = "out/table_sig_taxa_int.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)

#---- Diff Rel Ab Background ----
#aldex.clr
clr_bac <- aldex.clr(otu_table_decontam, conds = meta_decontam$Background, 
                     gamma = 0.5, mc.samples = 1000, denom = "none")

clr_bac@scaleSamps

#plot scale
graph.df <- data.frame("scale" = clr_bac@scaleSamps[,1], group = meta_decontam$Background)

ggplot(graph.df, aes(x = group, y = scale))+
  geom_boxplot() +
  labs(title = paste("Differential Relative Abundance: High vs Normal"), 
       x = "Total Bile Acid Level", y = "Scale")

#Running ttest for Intervention
x.ss.tt.bac <- aldex.ttest(clr_bac)

#Running effect
x.ss.e.bac <- aldex.effect(clr_bac, CI=T)

#combine for a plot
x.ss.all.bac <- cbind(x.ss.e.bac, x.ss.tt.bac)

# look at plot to see if positive/negative effect size corresponds to which size
#save graph
tiff(filename = "viz/background_diff_relative_abundance.tiff", width = 1300, height = 1300, units = "px", pointsize = 10, res = 300)

# look at plot to see if positive/negative effect size corresponds to which size
bac_dra <- aldex.plot(x.ss.all.bac, main = "Differential Relative Abundance by Background", type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

#close tiff
dev.off()

#effect size plot
bac_effect <- aldex.plot(x.ss.all.bac, type="MW", test="wicox", main='effect plot', 
                         cutoff.pval = 0.05)

#save graph
ggsave(bac_dra, file = "viz/intervention_diff_relative_abundance.tiff", dpi = 600)

sig_aldex_bac <- x.ss.all.bac %>%
  filter(wi.eBH < 0.05)

View(sig_aldex_bac)
#zero significant features

#make a table of all correct p-vallues
all_aldex2_bac <- x.ss.all.bac %>%
  rownames_to_column(var = "OTU") %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# make a table of significant corrected p-values
#sig_aldex2_bac2 <- x.ss.all.bac %>%
  #rownames_to_column(var = "OTU") %>%
  #filter(wi.eBH < 0.05) %>%
  #arrange(effect, wi.eBH) %>%
  #dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# add in previously formed taxa information to complete the table
all_aldex2_bac <- left_join(all_aldex2_bac, taxa_info)
save(all_aldex2_bac, file = "R/table_all_taxa_bac.RData")

#save in csv or xlsx tables so you can use the table for build a plot
write.table(all_aldex2_bac, file = "out/table_all_taxa_bac.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)

#---- Diff Rel Ab GLM/KW ----

# Check if Background and Intervention are factors
is.factor(meta_decontam$Sample_Type) # prints FALSE
is.factor(meta_decontam$Background) # prints FALSE
is.factor(meta_decontam$Intervention) # prints FALSE

# Convert Background and Intervention to factors if not already
meta_decontam$Sample_Type <- factor(meta_decontam$Sample_Type)
meta_decontam$Background <- factor(meta_decontam$Background)
meta_decontam$Intervention <- factor(meta_decontam$Intervention)

# Recheck if Background and Intervention are factors
is.factor(meta_decontam$Sample_Type) # prints TRUE
is.factor(meta_decontam$Background) # prints TRUE
is.factor(meta_decontam$Intervention) # prints TRUE

# Check the levels of your categorical variable
levels(meta_decontam$Sample_Type)
levels(meta_decontam$Background)
levels(meta_decontam$Intervention)

# Reorder the factor levels to explicitly choose the reference level
meta_decontam$Sample_Type <- relevel(meta_decontam$Sample_Type, ref = "NormCon")

meta_decontam$Background <- relevel(meta_decontam$Background, ref = "normal total bile acid level")

meta_decontam$Intervention <- relevel(meta_decontam$Intervention, ref = "control")

# Check the levels of your categorical variable
levels(meta_decontam$Sample_Type)
levels(meta_decontam$Background)
levels(meta_decontam$Intervention)

# make mm
mm1 <- model.matrix(~ Sample_Type, meta_decontam)
mm2 <- model.matrix(~ Background + Intervention, meta_decontam)

dim(otu_table_decontam)  # Returns features, samples (785, 24)
dim(mm1)             # Returns samples, features (24, 4) 
dim(mm2)             # Returns samples, features (24, 3) 
# NEED samples to be the same

# Step 2: Perform CLR transformation using ALDEx2
# Here, `mc.samples=1000` means we are running 1000 Monte Carlo Dirichlet instances to compute the CLR transformation.
treatment_clr1 <- aldex.clr(otu_table_decontam, mm1,
                           gamma=0.5, mc.samples=128, denom="all", verbose=T)

treatment_clr2 <- aldex.clr(otu_table_decontam, mm2,
                           gamma=0.5, mc.samples=128, denom="all", verbose=T)

save.image(file = file.path("R/aldex2_tmp.RData"))

treatment_results <- aldex.glm(treatment_clr1, mm1, fdr.method = "holm", verbose = TRUE)
# Error in family$family : $ operator is invalid for atomic vectors

treatment_results <- aldex.glm(treatment_clr2, mm2, fdr.method = "holm", verbose = TRUE)
# Error in family$family : $ operator is invalid for atomic vectors

#glm.eff<- aldex.glm.effect(x.glm)

#Kruskal Wallis
treatment.kw <- aldex.kw(treatment_clr1)
# Error in aldex.kw(treatment_clr1) : mismatch btw 'length(conditions)' and 'length(names(clr))'. len(condtitions): 96 len(names(clr)): 24

treatment.kw <- aldex.kw(treatment_clr2)
# Error in aldex.kw(treatment_clr2) : mismatch btw 'length(conditions)' and 'length(names(clr))'. len(condtitions): 72 len(names(clr)): 24

# save image
save.image(file = file.path("R/aldex2kw.RData"))

#---- Pairwise comparisons ----
# Generate all pairwise combinations of Sample_Type levels
is.factor(meta_decontam$Sample_Type) # prints FALSE
meta_decontam$Sample_Type <- factor(meta_decontam$Sample_Type)
is.factor(meta_decontam$Sample_Type) # prints TRUE

sample_types <- levels(meta_decontam$Sample_Type)
pairwise_comparisons <- combn(sample_types, 2, simplify = FALSE)

# List to store results for each comparison
comparison_results <- list()

for (pair in pairwise_comparisons) {
  # Subset the OTU table and Sample_Type for the current pair
  subset_indices <- meta_decontam$Sample_Type %in% pair
  subset_otu_table <- otu_table_decontam[, subset_indices]
  subset_conds <- meta_decontam$Sample_Type[subset_indices]
  
  # Convert conditions to character
  subset_conds <- as.character(subset_conds)
  
  # Run aldex.clr for the current pair
  pw_clr <- aldex.clr(subset_otu_table, 
                      conds = subset_conds, 
                      gamma = 0.5, 
                      mc.samples = 1000)
  
  # Run aldex.ttest and aldex.effect for the current pair
  result_ttest <- aldex.ttest(pw_clr)
  result_effect <- aldex.effect(pw_clr, CI = TRUE)
  
  # Combine results for the pair
  combined_results <- cbind(result_ttest, result_effect)
  
  # Define pair_name
  pair_name <- paste(pair, collapse = "_vs_")
  
  # Open a tiff device
  tiff(filename = paste0("viz/", pair_name, "_diff_relative_abundance.tiff"), width = 1500, height = 1600, units = "px", pointsize = 10, res = 300)
  
  # Create the plot
  pair_plot <- aldex.plot(combined_results, main = paste("Differential Relative Abundance:", pair_name), type = "MW", test = "wilcox", called.cex = 1, cutoff = 0.05)
  
  # Close the tiff device
  dev.off()
  
  # Filter and arrange results for significance
  stat_combined_results <- combined_results %>%
    rownames_to_column(var = "OTU") %>%
    arrange(effect, wi.eBH) %>%
    dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH, wi.ep, wi.eBH)
  
  # Add taxonomy information
  taxa_info <- data.frame(tax_table(ps_decon)) %>% rownames_to_column(var = "OTU")
  combined_results_taxa <- left_join(stat_combined_results, taxa_info)
  
  # Save all results to a file
  write.csv(
    combined_results_taxa,
    file = paste0("out/", pair_name, "_table_all_taxa.csv"),
    row.names = FALSE
  )
  
  # Filter significant results and save
  sig_aldex_pw <- stat_combined_results %>%
    filter(wi.eBH < 0.05) %>%
    left_join(taxa_info)
  
  write.csv(
    sig_aldex_pw,
    file = paste0("out/", pair_name, "_sig_table_all_taxa.csv"),
    row.names = FALSE
  )
  
  # Store results in the list
  comparison_results[[pair_name]] <- combined_results
  
  # Save all combined results to separate CSV files
  write.csv(
    combined_results,
    file = paste0("out/", pair_name, "_aldex_results.csv"),
    row.names = FALSE
  )
}

# View or inspect results
print(comparison_results)



#---- Below is code to create a graph of Significant Taxa----
#mutating the table
sig.HANA <- read.csv(file = "out/HighAlc_vs_NormAlc_sig_table_all_taxa_int.csv", sep = ",")

sig.HANA <- sig.HANA %>%
  select(Phylum, Family, Order, Genus, wi.eBH, effect)

ggplot(sig.HANA, aes(x = wi.eBH, y = Family)) +
  theme_bw() +
  geom_point(alpha = 0.8, size = 6) +
  facet_wrap(~Phylum, nrow = 1) +
  labs(title = "Significant Taxa High Alcohol vs Normal Control", 
       x = "p value, Wilcox + BH", y = "Family") +
  theme(axis.title = element_blank(), axis.text = element_text(size = 12),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24), legend.position = "bottom") #+
  #guides(color = guide_legend(nrow = 1)) #+
  # scale_color_manual(values = c("Actinobacteria" = "#2D3C59",
                                #"Firmicutes" = "#156082"))


