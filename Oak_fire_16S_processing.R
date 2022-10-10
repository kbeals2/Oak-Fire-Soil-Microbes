#### BACTERIAL & ARCHAEAL 16S AMPLICON PROCESSING PIPELINE USING DADA2 based off of official DADA2 tutorial (https://benjjneb.github.io/dada2/ITS_workflow.html) and detailed tutorial by Happy Belly Bioinformatics (https://astrobiomike.github.io/amplicon/dada2_workflow_ex#the-data) ####

# Fist, follow the instructions at https://benjjneb.github.io/dada2/dada-installation.html to install DADA2.


(.packages())
install.packages("pacman") 
library(pacman)
p_load("ShortRead", "Biostrings", "dada2")


#### 1) Set working directory ####
# set to wherever the fastq files are stored on your computer.
path <- "/Users/kendallb/Documents/Research/Collaborations/Oak_fire_2019_prj/Oak_fire_2019_microbiome_data/Oak_fire_2019_16Sseqs_no_lowmod"
list.files(path)


#### 2) Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz ####
FWD_reads <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
REV_reads <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
(sample.names <- sapply(strsplit(basename(FWD_reads), "_"), `[`, 1)) #sample ID occurs before the first _ in the fastq file


#### 3) Inspect quality profiles ####
# quality plots allow for visualizing the quality of the sequences of each sample; want quality score of 30 or above
# bases on the x-axis
# quality score on the y-axis
# green line = median quality score at that base position
# orange line = quartiles
plotQualityProfile(FWD_reads[5])
# Quality drops below q30 after ~260 bases

plotQualityProfile(REV_reads[1])
# Quality drops below q30 at about 220 bases

#### 4) Filter and Trim ####
# For this dataset, we will use standard filtering paraments: maxN = 0 (DADA2 requires sequences contain no Ns), truncQ = 2,  rm.phix = TRUE and maxEE = 2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note that the primers are still are on the forward & reverse reads, so we remove them with the trimLeft parameter(length of FWD primer, length of REV primer). We used the 341 forward primer: CCTACGGGNGGCWGCAG (17 bp in length) & 785 reverse primer: GACTACHVGGGTATCTAATCC (21 bps in length).
filt_FWD <- file.path(path, "filtered", paste0(sample.names, "_F_filtered.fastq.gz"))
filt_REV <- file.path(path, "filtered", paste0(sample.names, "_R_filtered.fastq.gz"))
names(filt_FWD) <- sample.names
names(filt_REV) <- sample.names
out <- filterAndTrim(FWD_reads, filt_FWD, REV_reads, filt_REV, trimLeft = c(17, 21), truncLen = c(260, 220),
                     maxN = 0, maxEE = c(2,5), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE, verbose = TRUE)

# Display # of reads present in each sample after filtering; aim for retention of 85-90% of reads
out


#### 5) Inspect forward read quality profiles after trimming. ####
plotQualityProfile(filt_FWD[1])
plotQualityProfile(filt_REV[1])

#### 6) Learn the error rates ####
err_FWD <- learnErrors(filt_FWD, multithread = TRUE) # forward filtered reads
err_REV <- learnErrors(filt_REV, multithread = TRUE) # reverse filtered reads

# Visualize the estimated error rates as a sanity check
# red line = expected error rate based on quality score (note that this decreases as the quality increases)
# black line = estimated error rate after convergence of the machine-learning algorithm
# black dots = observed error frequency in our samples
# black dots should align with black line
plotErrors(err_FWD, nominalQ = TRUE) # looking for the points to align with the black line
plotErrors(err_REV, nominalQ = TRUE)


#### 7) Dereplicate and denoise ####
derep_FWD <- derepFastq(filt_FWD, verbose = TRUE)
derep_REV <- derepFastq(filt_REV, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derep_FWD) <- sample.names
names(derep_REV) <- sample.names


#### 8) Sample inference ####
# Here is where the real power of DADA2 comes in.
# DADA2 infers biological sequences by incorporating the consensus quality profiles & abundances of each unique sequence, and figures out which sequence is more likely to be of biological origin or more likely to be spurious.
# The default of the dada function processes each sample independently, but pseudo-pooling allows samples to be processed independently after sharing information between samples. Pooling information across samples can increase sensitivity to sequence variants that may be present at lower frequencies. 
dada_FWD <- dada(derep_FWD, err = err_FWD, multithread = TRUE, pool = "pseudo") #how many unique sequences? can also do pool=TRUE to pool across all samples to detect low abundance variants
dada_REV <- dada(derep_REV, err = err_REV, multithread = TRUE, pool = "pseudo")


#### 9) Merge paired reads ####
# aligning denoised forward reads with reverse-complement of the corresponding denoised reverse reads
mergers <- mergePairs(dada_FWD, derep_FWD, dada_REV, derep_REV, verbose = TRUE)

# BGB: 59% merged
# FCM: 87% merged
# LCM: 81% merged
# PL441: 69% merged
# RGB: 85% merged
# RS441: 73% merged
# SCM: 87% merged

# Inspect the merger data.frame from the first sample
head(mergers[[1]]) # Can blast some of the sequences to double check they are what you think they are: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome


#### 10) Construct a sequence table ####
seq_table <- makeSequenceTable(mergers)
dim(seq_table) # (number of samples, number of amplicon sequence variants ASVs) (found 3,376 ASV across the 7 samples)


#### 11) Remove chimeras ####
# chimeras are generally considered contaminants because they can be interpreted as a novel sequence when they are in fact an error
seq_table_nochim <- removeBimeraDenovo(seq_table, method = "consensus", multithread = TRUE, verbose = TRUE)
# found 464 chimeras out of 3,376 sequences (only 13.7% identified as chimeric. Good!)

# We don't know if these 464 chimeras held a lot in terms of abundance. To find out, we can calculate the proportion of sequences retained after chimeras removed.
sum(seq_table_nochim)/sum(seq_table)
# We retained 97% of sequence abundance. Great!

# ** Most of the reads should remain after chimera removal. If most of the reads were removed as chimeric, it likely requires revisiting the upstream processing (it usually means the primers were not properly removed).


#### 12) Track reads through the pipeline to verify everything worked as expected ####
getN <- function(x) sum(getUniques(x))
(summary_table <- data.frame(row.names = sample.names, 
                            input = out[, 1],
                            filtered = out[, 2], 
                            denoised_FWD = sapply(dada_FWD, getN),
                            denoised_REV = sapply(dada_REV, getN), 
                            merged = sapply(mergers, getN),
                            non_chim = rowSums(seq_table_nochim), 
                            final_perc_reads_retained = round(rowSums(seq_table_nochim)/out[, 1]*100, 1)))

# Looks good! The most amount of reads were removed during filtering, as expected. No large removal of reads in any sample after filtering. Retained on average 60% of reads.
# Besides filtering, there should be no step in which a majority of reads are lost.


#### 13) Assign bacteria and archaea taxonomy with the Ribosomal Database Projects (RDP) training set database from https://zenodo.org/record/801828#.XNHerqZ7lSM ####
taxa <- assignTaxonomy(seq_table_nochim, "rdp_train_set_16.fa.gz", multithread = TRUE) # if filled with NAs; try including tryRC = TRUE
taxa <- addSpecies(taxa, "rdp_species_assignment_16.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


#### 14) Generate a fasta file, a count table, and a taxonomy table ####
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seq_table_nochim)
asv_headers <- vector(dim(seq_table_nochim)[2], mode = "character")

for (i in 1:dim(seq_table_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# make count table:
asv_tab <- t(seq_table_nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "Oak_fire_2019_Bact_ASVs_counts_no_lowmod.tsv", sep = "\t", quote = F, col.names = NA)

# make taxa table:
asv_taxa <- taxa
row.names(asv_taxa) <- sub(">", "", asv_headers)
write.table(asv_taxa, "Oak_fire_2019_Bact_ASVs_taxonomy_no_lowmod.tsv", sep = "\t", quote = F, col.names = NA)
