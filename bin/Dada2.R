# check wd
getwd()

# set wd if needed
setwd("/Users/sls6550/work/hTBA_Alc_Singh")

#check wd to confirm change
getwd()

# install packages if needed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.19")

#---- dada2 ----
library(dada2)

# SET VARIABLES

## modify these variables:

# path to filtered and cleaned reads
CLEANEDPATH = "data/fastp_filt" # CHANGE ME to the directory containing the fastp filtered

# path to output tables
OUTPATH = "out" # CHANGE ME to desired directory for out

RDATA = "R" # CHANGE ME to desired directory for R files

# database to silva training set
# downloaded from: https://benjjneb.github.io/dada2/training.html
SILVADB = "ref/silva_nr_v132_train_set.fa.gz"
SILVASPASSIGN = "ref/silva_species_assignment_v132.fa.gz"

# paired end patterns
FILTEREDF = "_fastp_1.fastq"
FILTEREDR = "_fastp_2.fastq"

# Forward and reverse fastq filenames have format: SAMPLENAME_fastp_1.fastq and SAMPLENAME_fastp_2.fastq
fnFs <- sort(list.files(CLEANEDPATH, pattern=FILTEREDF, full.names = TRUE))
fnRs <- sort(list.files(CLEANEDPATH, pattern=FILTEREDR, full.names = TRUE))

# extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# check to make sure that the lengths of both files are the same and that they match
fwdNames <- sapply(strsplit(basename(fnFs), fnFs), `[`, 1)
revNames <- sapply(strsplit(basename(fnRs), fnRs), `[`, 1)

# error catch
if(length(fwdNames) != length(revNames)) {
  stop("The number of forward and reverse files do not match.")
}

# inspect read quality profiles
plotQualityProfile(fnFs[2:3])
plotQualityProfile(fnRs[2:3])
# still some issues with read quality loss on reverse reads

# filtering via Dada2
filtFs <- file.path(CLEANEDPATH, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(CLEANEDPATH, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
# trimming fwd reads after 240 bp, rev after 160 bp, removing all Ns, removing any reads with an expected error count > 2, remove PhiX standard, 

head(out)

#reinspect
plotQualityProfile(filtFs[2:3])
plotQualityProfile(filtRs[2:3])

#looks better!

#error rates for each read (forward and reverse)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plot error rates for visualization purposes
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# save progress
save.image(file = file.path(RDATA, "error_learning.RData"))

#dereplicate
#can use to decrease computational time
#not necessary
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

#sample interference denosing
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#view dada2 class objects for fwd and rev reads
dadaFs[[1]]
dadaRs[[1]]

# save image
save.image(file = file.path(RDATA, "denoised.RData"))

#merge paired end reads together for the entire denoised seq to form contigs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
View(mergers[[1]])

#construct ASV table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # 25 1273

#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab))) 
# 252 253 254 255 
# 413 837  22   1 

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) # 25 794

# save ASV table as RData object
saveRDS(seqtab.nochim, file = file.path(RDATA, "asv_table.rds"))

# write out cleaned count table
if (!is.data.frame(seqtab.nochim)) {
  seqtab.nochim_df <- as.data.frame(seqtab.nochim)
}

# save table
seq_path <- file.path(OUTPATH, "seqtab_counts.csv")
write.csv(seqtab.nochim_df, seq_path, row.names = TRUE)

#total percentage of chimera abundance
(1 - (sum(seqtab.nochim)/sum(seqtab))) * 100 # in percent

#track progress
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("Raw_Reads", "Trimmed_Reads", "Dada2_FWD", "Dada2_REV", "Merged_Read_Count", "Final_Read_Count")
View(track)

# write the reads lost to file
track_path <- file.path(OUTPATH, "track-reads-filtering.txt")
write.table(track, track_path, sep = "\t")

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
#rownames(track) <- sample.names
#head(track)

# assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, SILVADB, multithread=TRUE)

# species assignment with exact matches ONLY
taxa <- addSpecies(taxa, SILVASPASSIGN)

# write taxonomy table to file
write.table(taxa, file = file.path(OUTPATH, "tax-table.txt"), sep = "\t", row.names = FALSE)

# save tax table as RData object
saveRDS(taxa, file = file.path(RDATA, "tax_table.rds"))

# print finished message
cat("done!")
