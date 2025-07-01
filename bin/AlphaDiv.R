# Alpha Diversity
# Samantha Seibel
# Singh TBA Alcohol Mouse Study
# November 27th, 2024

# ---- alpha diversity  ----
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(car)
library(MASS)

#get alpha diversity metrics

# check wd
getwd()

# set wd if needed
setwd("/Users/sls6550/work/hTBA_Alc_Singh")

# check wd to confirm change
getwd()

# SET PATHS
RDATA = "R"

# load ps from relative abundance without controls
ps_rel <- readRDS("R/ps_decontam.rds")

ntaxa(ps_rel) #788

# glom taxa
glom <- ps_rel %>%
  tax_glom(taxrank = "Genus")

# make df with SampleID as rownames
sampdf <- sample_data(glom) %>%
  data.frame() %>%
  rownames_to_column(var = "SampleID")

sampdf <- sampdf %>%
  mutate(Group = factor(Group, levels = c("NormCon", "NormAlc", "HighCon", "HighAlc")))

# get diversity
alpha <- estimate_richness(glom, measures = c("Shannon", "Simpson")) %>%
  # make ID
  rownames_to_column(var = "SampleID") %>%
  merge(sampdf, by = "SampleID")

# get summary statistics for Shannon
sum <- alpha %>%
  group_by(Intervention, Background) %>% get_summary_stats(Shannon, type = "mean_sd")
write.table(sum, file = "out/shannon-sum-stats.txt", sep = "\t", row.names = FALSE)

# can already see unequal standard deviations...

# exploratory smoothed histogram
ggdensity(alpha, x = "Shannon", fill = "Intervention")
ggdensity(alpha, x = "Shannon", fill = "Background")
ggdensity(alpha, x = "Shannon", fill = "Sample_Type")

# looks skewed

# exploratory plots
hist(alpha$Shannon) #defintely skewed
ggboxplot(alpha, x = "Intervention", y = "Shannon", facet.by = "Background")
#ggboxplot(alpha, x = "Background", y = "Shannon", facet.by = "Intervention")

# variance again is different

# check on outliers for shannon
# table for determining outlier range
st <- alpha %>% get_summary_stats(Shannon, type = "mean_sd") %>%
  mutate(sdx3 = sd * 3,
         outhi = mean + sdx3,
         outlo = mean - sdx3)

# SD of Shannon is 0.714; Range is 1.555 to 5.839
outliershannon <- alpha$SampleID[alpha$Shannon < st$outlo] # none
alpha$SampleID[alpha$Shannon > st$outhi] # no upper outliers

# confirmed visually on alpha table

# only need to filter if there are outliers
#alpha %>%
#filter(ID %in% outliershannon)

# build test model for linear assumption of Intervention
# Convert Intervention to factor
alpha$Intervention <- factor(alpha$Intervention)
# Convert Background to factor
alpha$Background <- factor(alpha$Background)
# Convert Sample_Type to factor
alpha$Sample_Type <- factor(alpha$Sample_Type)

# set reference level
alpha$Intervention <- relevel(alpha$Intervention, ref = "control")

# run linear model
testint <- lm(Shannon ~ Intervention, data = alpha)

#%>%
# filter(isoutlier == F)) # more normal if outliers removed
summary(testint) #p value = 0.548 Adj R2 = -0.02807
hist(resid(testint))
qqnorm(resid(testint))
qqline(resid(testint))

# build test model for linear assumption of Background
# set reference level
alpha$Background <- relevel(alpha$Background, ref = "normal total bile acid level")
testback <- lm(Shannon ~ Background, data = alpha)
#%>%
#filter(isoutlier == F)) # more normal if outliers removed

summary(testback) #p value =0.345 and Adj R2 =-0.002947
hist(resid(testback))
qqnorm(resid(testback))
qqline(resid(testback))

# build test model for linear assumption of Background
# set reference level
alpha$Sample_Type <- relevel(alpha$Sample_Type, ref = "NormCon")
testST <- lm(Shannon ~ Sample_Type, data = alpha)
#%>%
#filter(isoutlier == F)) # more normal if outliers removed

summary(testST) #p value =0.345 and Adj R2 =-0.002947
hist(resid(testST))
qqnorm(resid(testST))
qqline(resid(testST))


# non normal visually
# shapiro-wilkes to confirm
shapiro.test(alpha$Shannon)
# W = 0.86581, p-value = 0.00436
# not normal, need to either transform or use non-parametric options

shapiro.test(alpha$Simpson)
# W = 0.79984, p-value = 0.0002928
# not normal, need to either transform or use non-parametric options

# homogeneity of variance
leveneTest(Shannon ~ Background, data = alpha) # not significant (0.5276), equal!
leveneTest(Shannon ~ Intervention, data = alpha) # not significant (0.05003), equal! but almost
leveneTest(Shannon ~ Sample_Type, data = alpha) # significant (0.0279), unequal!
leveneTest(Simpson ~ Background, data = alpha) # not significant (0.2763), equal!
leveneTest(Simpson ~ Intervention, data = alpha) # significant (0.02108), unequal!
leveneTest(Simpson ~ Sample_Type, data = alpha) # significant (0.01209), unequal!

# transformation of data
alpha$log_Shannon <- log(alpha$Shannon)
alpha$log_Simpson <- log(alpha$Simpson)

# SW test
shapiro.test(alpha$log_Shannon) #nope still not normal p-value = 4.162e-06
shapiro.test(alpha$log_Simpson) #nope still not normal p-value = 7.446e-07

# Levene's
leveneTest(log_Shannon ~ Background, data = alpha) # not significant (0.2776), equal!
leveneTest(log_Shannon ~ Intervention, data = alpha) # not significant (0.08544), equal!
leveneTest(log_Shannon ~ Sample_Type, data = alpha) # significant (0.04085), unequal!

leveneTest(log_Simpson ~ Background, data = alpha) # not significant (0.1961), equal!
leveneTest(log_Simpson ~ Intervention, data = alpha) # not significant (0.0613), equal!
leveneTest(log_Simpson ~ Sample_Type, data = alpha) # significant (0.03605), unequal!

# box-cox transformation
bc <- boxcox(lm(Shannon ~ 1, data = alpha))
lambda <- bc$x[which.max(bc$y)]
alpha$bc_Shannon <- (alpha$Shannon^lambda - 1) / lambda

# SW test
shapiro.test(alpha$bc_Shannon) #normal; W = 0.96003, p-value = 0.439

# Levene's
leveneTest(bc_Shannon ~ Background, data = alpha) #not significant (0.9285), equal!
leveneTest(bc_Shannon ~ Intervention, data = alpha) # significant (0.03825), unequal!

# not going to do parametric, going to non-parametric

# Mann U test Intervention
wilcox_int <- wilcox.test(Shannon ~ Intervention,  alternative = "two.sided", data = alpha)
wilcox_int #p-value = 0.8428 non-significant


# Mann U test Background
wilcox_back <- wilcox.test(Shannon ~ Background,  alternative = "two.sided", data = alpha)
wilcox_back #p-value = 0.2913 non significant

# Kruskal Wallis
kruskalShannon <- kruskal.test(Shannon ~ Sample_Type, data = alpha)
kruskalShannon
# p-value = 0.528
# non significant

kruskalSimpson <- kruskal.test(Simpson ~ Sample_Type, data = alpha)
kruskalSimpson
# p-value = 0.6429
# non significant

# Calculate pairwise comparisons Dunn's Test if significant KW
# using bonferroni's correct as there are less comparisons
# Dunn <- dunnTest(Shannon ~ Sample_Type, data=alpha, method="bonferroni")
# Dunn 
# p-values are all non-significant for Sample_Type

# colors
pal <- c("#b5e48c", "#76c893", "#168aad", "#1e6091")

# plots

# violin shannon all 4 groups
ggplot(alpha, aes(x=Sample_Type, y = Shannon, fill = Sample_Type)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_jitter(size = 3) +
  labs(x = "Sample Type", y = "Shannon", title = "Alpha Diversity by Alcohol Intervention and Background Bile Acid Level", fill = "Sample Type") +
  scale_fill_manual(values = pal) +
  theme_classic() +
  scale_x_discrete(labels = c(
    "HighAlc" = "High TBA + Alcohol",
    "HighCon" = "High TBA",
    "NormAlc" = "Low TBA + Alcohol", 
    "NormCon" = "Low TBA")) +
  theme(text = element_text(size = 15),
        axis.text= element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(angle = 90, size = 18),
        plot.title = element_text(size = 24, hjust = 0.5),
        legend.position = "none", panel.grid = element_blank())

#save plot
ggsave("viz/alpha_shan_violin.pdf", plot = last_plot(), height = 10, width = 12, units = "in")

# simpson all 4 groups
#ggplot(alpha, aes(x=Sample_Type, y = Simpson, fill = Sample_Type)) +
  #geom_violin(trim = FALSE, alpha = 0.8) +
  #geom_jitter(size = 1) +
  #labs(x = "Sample Type", y = "Simpson", title = "Alpha Diversity by Alcohol Intervention and Background Bile Acid Level", fill = "Sample Type") +
  #scale_fill_manual(values = earth_pal) +
  #theme_bw() +
  #theme(text = element_text(size = 15),
        #axis.text= element_text(size = 24),
       #axis.title.y = element_text(size = 24),
        #axis.title.x = element_text(size = 24),
        #axis.text.x = element_text(size = 18),
        #plot.title = element_text(size = 24, hjust = 0.5),
        #legend.position = "none")

#ggsave("viz/alpha_simpson_ST.pdf", plot = last_plot(), height = 10, width = 12, units = "in")

# shannon 
# shannon all 4 groups
ggplot(alpha, aes(x=Sample_Type, y = Shannon, fill = Sample_Type)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8), size = 2, alpha = 0.8) +
  labs(x = "Sample Type", y = "Shannon", title = "Alpha Diversity by Alcohol Intervention and TBA Level", fill = "Sample Type") +
  scale_fill_manual(values = pal) +
  theme_classic() +
  scale_x_discrete(labels = c(
    "HighAlc" = "High TBA + Alcohol",
    "HighCon" = "High TBA",
    "NormAlc" = "Low TBA + Alcohol", 
    "NormCon" = "Low TBA")) +
  theme(text = element_text(size = 15),
        axis.text= element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        plot.title = element_text(size = 24, hjust = 0.5),
        legend.position = "none")

#save plot
ggsave("viz/alpha_shannon_box.pdf", plot = last_plot(), height = 10, width = 12, units = "in")

# Shannon
ggplot(alpha, aes(x = Intervention, y = Shannon, fill = Background)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 2, alpha = 0.8) +
  labs(x = "Alcohol Intervention", y = "Shannon", title = "Alpha Diversity by Alcohol Intervention and Total Bile Acid Level", fill = "Total Bile Acid") +
  theme_classic() +
  scale_fill_manual(labels = c("high total bile acid level" = "High",
    "normal total bile acid level" = "Low"), values = c("high total bile acid level" = "#006400", "normal total bile acid level" = "#b5e48c")) +
  scale_x_discrete(labels = c(
    "alcohol" = "Yes",
    "control" = "No")) +
  theme(text = element_text(size = 15),
        axis.text= element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        plot.title = element_text(size = 24, hjust = 0.5))

ggsave("viz/alpha_shannon.pdf", plot = last_plot(), height = 10, width = 12, units = "in")

#---- Extra alpha simpson----
ggplot(alpha, aes(x=Intervention, y = Simpson, fill = Background)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8), size = 2, alpha = 0.8) +
  labs(x = "Background", y = "Simpson", title = "Alpha Diversity by Alcohol Intervention and Background Bile Acid Level", fill = "Intervention") +
  scale_fill_manual(values = c("#d4bbff", "#31135e")) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text= element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        plot.title = element_text(size = 20, hjust = 0.5))

ggsave("viz/alpha_simpson_split.pdf", plot = last_plot(), height = 10, width = 12, units = "in")

# three alpha diversity metrics

plot_richness(ps_rel, x = "Sample_Type", measures = c("Chao1", "Shannon", "Simpson")) 

ggsave("viz/alpha.pdf", plot = last_plot(), height = 4, width = 6, units = "in")
