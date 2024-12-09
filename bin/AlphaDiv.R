# Alpha Diversity
# Samantha Seibel
# Singh TBA Alcohol Mouse Study
# November 27th, 2024

# ---- alpha diversity  ----

library(ggpubr)
library(car)
library(MASS)
library(FSA)

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

# make df with SampleID as rownames
sampdf <- sample_data(ps_rel) %>%
  data.frame() %>%
  rownames_to_column(var = "SampleID")

# get diversity
alpha <- estimate_richness(ps_rel, measures = c("Shannon", "Simpson")) %>%
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

# non normal visually
# shapiro-wilkes to confirm
shapiro.test(alpha$Shannon)
# W = 0.81447, p-value = 0.0004993
# not normal, need to either transform or use non-parametric options

# homogeneity of variance
leveneTest(Shannon ~ Background, data = alpha) # not significant (0.6592), equal!
leveneTest(Shannon ~ Intervention, data = alpha) # not significant (0.3703), equal!

# transformation of data
alpha$log_Shannon <- log(alpha$Shannon)

# SW test
shapiro.test(alpha$log_Shannon) #nope still not normal p-value = 1.689e-05

# Levene's
leveneTest(log_Shannon ~ Background, data = alpha) #not significant (0.476), equal!
leveneTest(log_Shannon ~ Intervention, data = alpha) #not significant (0.3133), equal!

# box-cox transformation
bc <- boxcox(lm(Shannon ~ 1, data = alpha))
lambda <- bc$x[which.max(bc$y)]
alpha$bc_Shannon <- (alpha$Shannon^lambda - 1) / lambda

# SW test
shapiro.test(alpha$bc_Shannon) #nope still not normal p-value = 0.009834

# Levene's
leveneTest(bc_Shannon ~ Background, data = alpha) #not significant (0.8985), equal!
leveneTest(bc_Shannon ~ Intervention, data = alpha) #not significant (0.4334), equal!

# not going to do parametric, going to non-parametric

# Mann U test Intervention
wilcox_int <- wilcox.test(Shannon ~ Intervention,  alternative = "two.sided", data = alpha)
wilcox_int #p-value = 0.9774 non-significant


# Mann U test Background
wilcox_back <- wilcox.test(Shannon ~ Background,  alternative = "two.sided", data = alpha)
wilcox_back #p-value = 0.4095 non significant

# Kruskal Wallis
kruskalShannon <- kruskal.test(Shannon ~ Sample_Type, data = alpha)
kruskalShannon
# p-value = 0.3906
# non significant

kruskalSimpson <- kruskal.test(Simpson ~ Sample_Type, data = alpha)
kruskalSimpson
# p-value = 0.3968
# non significant

# Calculate pairwise comparisons Dunn's Test if significant KW
# using bonferroni's correct as there are less comparisons
# Dunn <- dunnTest(Shannon ~ Sample_Type, data=alpha, method="bonferroni")
# Dunn 
# p-values are all non-significant for Sample_Type

#plots
ggplot(alpha, aes(x=Background, y = Shannon, fill = Intervention))+
  geom_violin(trim = FALSE, alpha = 0.8)+
  #geom_jitter(size = 1)+
  scale_fill_manual(values = c("#936639", "#656d4a"))+
  theme_bw()

ggplot(alpha, aes(x=Background, y = Simpson, fill = Intervention))+
  geom_violin(trim = FALSE, alpha = 0.8)+
  #geom_jitter(size = 1)+
  scale_fill_manual(values = c("#936639", "#656d4a"))+
  theme_bw()

ggsave("viz/alpha_shannon.pdf", plot = last_plot(), height = 20, width = 26, units = "in")

plot_richness(ps_rel, x = "Sample_Type", measures = c("Chao1", "Shannon", "Simpson"))
ggsave("viz/alpha.pdf", plot = last_plot(), height = 20, width = 26, units = "in")
