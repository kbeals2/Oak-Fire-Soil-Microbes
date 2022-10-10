
#### FUNGAL ITS AMPLICON PROCESSING PIPELINE USING DADA2 based off of official DADA2 tutorial (https://benjjneb.github.io/dada2/ITS_workflow.html) ####

# The front end of the ITS processing pipeline is slightly different from the 16S pipeline. The official DADA2 tutorial has a great summary explaining reasons for this that should be read beforehand. See above link.


# Fist, follow the instructions at https://benjjneb.github.io/dada2/dada-installation.html to install DADA2.

(.packages())
install.packages("pacman") 
library(pacman)
p_load("ShortRead", "Biostrings", "dada2")


#### 1) Set working directory ####
# set to wherever you the fastq files are stored on your computer.
path <- "/Users/kendallb/Documents/Research/Collaborations/Oak_fire_2019_prj/Oak_fire_2019_microbiome_data/Oak_fire_2019_ITSseqs_no_lowmod"
list.files(path)

#### 2) Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz ####
FWD_reads <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
REV_reads <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))


#### 3) Identify fungal ITS primers (5.8S-Fun and ITS4-Fun; can be found in the paper: https://aem.asm.org/content/aem/82/24/7217.full.pdf) ####
FWD_primer <- "AACTTTYRRCAAYGGATCWCT"
REV_primer <- "AGCCTCCGCTTATTGATATGCTTAART"


#### 4) Verify presence & orientation of primers ####
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

(FWD_orients <- allOrients(FWD_primer))
(REV_orients <- allOrients(REV_primer))

# Create new orientation with modified primers. Then re-run step 6 (primers should be in correct places). 
(FWD_orients <- allOrients(new_FWD_primer))
(REV_orients <- allOrients(new_REV_primer))


#### 5) The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. We are going to “pre-filter” the sequences just to remove those with Ns. ####
FWD_filt <- file.path(path, "filtered", basename(FWD_reads)) # Put N-filterd files in filtN/ subdirectory
REV_filt <- file.path(path, "filtered", basename(REV_reads))
filterAndTrim(FWD_reads, FWD_filt, REV_reads, REV_filt, maxN = 0, multithread = TRUE)


#### 6) Now we are going to count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations. Identifying and counting the primers on one set of paired end FASTQ files is sufficient, assuming all the files were created using the same library preparation, so we’ll just process the first sample. ####
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD_ForwardReads = sapply(FWD_orients, primerHits, fn = FWD_filt[[1]]), 
      FWD_ReverseReads = sapply(FWD_orients, primerHits, fn = REV_filt[[1]]), 
      REV_ForwardReads = sapply(REV_orients, primerHits, fn = FWD_filt[[1]]), 
      REV_ReverseReads = sapply(REV_orients, primerHits, fn = REV_filt[[1]]))

# Since the FWD primer is found on the REV reads in its FWD orientation & the REV primer is found on the FWD reads in its FWD orientation, you can just create new objects called new_FWD_primer (original REV primer) & new_REV_primer (original FWD primer). 
new_FWD_primer <- REV_primer
new_REV_primer <- FWD_primer

# Now, return to line 43. If the FWD primer of the forward reads were found in its forward orientation, and the REV primer were found in the reverse reads in its forward orientation, then we would not need to do any modification step here. Orientation mixups are common, so do not worry if you have to modify. 


#### 7) Remove primers using cutadapt (https://cutadapt.readthedocs.io/en/stable/index.html) ####
cutadapt <- "/Users/kendallb/Library/Python/3.8/bin/cutadapt" # this worked; houses version 3.0
system2(cutadapt, args = "--version") # Run shell commands from R

# Now create output file names for the cutadapt-ed files, and define the parameters to give the cutadapt command. The critical parameters are the primers, and they need to be in the right orientation, i.e. the FWD primer should have been matching the forward reads in its forward orientation, and the REV primer should have been matching the reverse reads in its forward orientation.
path_cut <- file.path(path, "cutadapt")
if(!dir.exists(path_cut)) dir.create(path_cut) # creates "cutadapt" file

FWD_cut <- file.path(path_cut, basename(new_FWD_primer))
REV_cut <- file.path(path_cut, basename(new_REV_primer))

FWD_RC <- dada2:::rc(new_FWD_primer)
REV_RC <- dada2:::rc(new_REV_primer)

# Trim FWD primer and the rc of REV primer off of R1 (forward reads)
R1_flags <- paste("-g", new_FWD_primer, "-a", REV_RC) 
# Trim REV primer and the rc of FWD primer off of R2 (reverse reads)
R2_flags <- paste("-G", new_REV_primer, "-A", FWD_RC) 

# Run Cutadapt
for(i in seq_along(FWD_reads)) {
  system2(cutadapt, args = c(R1_flags, R2_flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", FWD_cut[i], "-p", REV_cut[i], # output files
                             FWD_filt[i], REV_filt[i])) # input files
}


# Above loop resulted in a bug. Instead, run cutadapt on each sample individually in the shell, using the following:
cd /Users/kendallb/Documents/Research/Collaborations/Oak_fire_2019_prj/Oak_fire_2019_microbiome_data/Oak_fire_2019_ITSseqs_no_lowmod/filtered # set wd

# replace "BGB" with each sample name
cutadapt SCEA-BGB_R1.fastq.gz SCEA-BGB_R2.fastq.gz -o SCEA-BGB_R1_cut.fastq.gz -p SCEA-BGB_R2_cut.fastq.gz -g ^AGCCTCCGCTTATTGATATGCTTAART -a AGWGATCCRTTGYYRAAAGTT$ -G ^AACTTTYRRCAAYGGATCWCT -A AYTTAAGCATATCAATAAGCGGAGGCT$
  
  # -o: output file name of FWD reads
  # -p: output file name of REV reads
  # -g: FWD primer in normal orientation (to be trimmed off R1)
  # -a: reverse complement of REV primer (to be trimmed off R1)
  # -G: REV primer in normal orientation (to be trimmed off R2)
  # -A: reverse complement of FWD primer (to be trimmed off R2)
  
# Once this is complete, manually move cutadapted files to "cutadapt" folder  
  
# Re-assign FWD_cut & REV_cut to new cutadapted files
FWD_cut <- sort(list.files(path_cut, pattern = "_R1_cut.fastq.gz", full.names = TRUE))
REV_cut <- sort(list.files(path_cut, pattern = "_R2_cut.fastq.gz", full.names = TRUE))

# As a sanity check, we will count the presence of primers in the first cutadapted sample.
rbind(FWD_ForwardReads = sapply(FWD_orients, primerHits, fn = FWD_cut[[1]]), 
      FWD_ReverseReads = sapply(FWD_orients, primerHits, fn = REV_cut[[1]]), 
      REV_ForwardReads = sapply(REV_orients, primerHits, fn = REV_cut[[1]]), 
      REV_ReverseReads = sapply(REV_orients, primerHits, fn = REV_cut[[1]]))
# Great! Primers are no longer detected in the cutadapted reads.

# The primer-free sequence files are now ready to be analyzed through the DADA2 pipeline. Similar to the earlier steps of reading in FASTQ files, we read in the names of the cutadapt-ed FASTQ files and did some string manipulation to get the matched lists of forward and reverse fastq files.

# Forward and reverse fastq filenames have the format:
cut_FWD <- sort(list.files(path_cut, pattern = "_R1_cut.fastq.gz", full.names = TRUE))
cut_REV <- sort(list.files(path_cut, pattern = "_R2_cut.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cut_FWD, get.sample.name))
sample.names


#### 8) Inspect read quality profiles ####
# quality plots allow for visualizing the quality of the sequences of each sample; want quality score of 30 or above
# bases on the x-axis
# quality score on the y-axis
# green line = median quality score at that base position
# orange line = quartiles

# forward reads first
plotQualityProfile(cut_FWD[1])
# Quality drops below Q30 at around 230 bases; the sequences were cutadapt-ed to about 250 nucleotides in length (shown by red line of bottom of figures)

# then reverse reads
plotQualityProfile(cut_REV[1])
# Quality of reverse reads is not as good, drops below Q30 at around 200 bases. Reverse seqs are also cutadapt-ed to about 250 ncts and have exact same number of reads as their respective forward seqs!


#### 9) Filter and trim ####

# Assign the file names for the output of the filtered reads to be stored as fastq.gz files
filt_FWD <- file.path(path_cut, "filtered", basename(cut_FWD))
filt_REV <- file.path(path_cut, "filtered", basename(cut_REV))

# For this dataset, we will use standard filtering paraments: maxN = 0 (DADA2 requires sequences contain no Ns), truncQ = 2,  rm.phix = TRUE and maxEE = 2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Note: We enforce a minLen here, to get rid of spurious very low-length sequences. This was not needed for processing 16S sequences because truncLen already served that purpose.

out <- filterAndTrim(cut_FWD, filt_FWD, cut_REV, filt_REV, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

# View reads in compared to reads out   
out

#### 10) Learn the error rates (From here, the processing steps are the same as the 16S tutorial) ####
err_FWD <- learnErrors(filt_FWD, multithread = TRUE)
err_REV <- learnErrors(filt_REV, multithread = TRUE)

# Visualize the estimated error rates as a sanity check
# red line = expected error rate based on quality score (note that this decreases as the quality increases)
# black line = estimated error rate after convergence of the machine-learning algorithm
# black dots = observed error frequency in our samples
# black dots should align with black line
plotErrors(err_FWD, nominalQ = TRUE)
plotErrors(err_REV, nominalQ = TRUE)


#### 11) Dereplicate identical reads ####
derep_FWD <- derepFastq(filt_FWD, verbose = TRUE)
derep_REV <- derepFastq(filt_REV, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derep_FWD) <- sample.names
names(derep_REV) <- sample.names


#### 12) Sample inference ####
# Here is where the real power of DADA2 comes in.
# DADA2 infers biological sequences by incorporating the consensus quality profiles & abundances of each unique sequence, and figures out which sequence is more likely to be of biological origin or more likely to be spurious.
# The default of the dada function processes each sample independently, but pseudo-pooling allows samples to be processed independently after sharing information between samples. Pooling information across samples can increase sensitivity to sequence variants that may be present at lower frequencies. 
dada_FWD <- dada(derep_FWD, err = err_FWD, multithread = TRUE, pool = "pseudo")
dada_REV <- dada(derep_REV, err = err_REV, multithread = TRUE, pool = "pseudo")


#### 13) Merge paired reads ####
# aligning denoised forward reads with reverse-complement of the corresponding denoised reverse reads
mergers <- mergePairs(dada_FWD, derep_FWD, dada_REV, derep_REV, verbose = TRUE)

# BGB: 91% merged
# FCM: 96% merged
# LCM: 93% merged
# PL441: 93% merged
# RGB: 94% merged
# RS441: 94% merged
# SCM: 98% merged


#### 14) Construct sequence/count table ####
seq_table <- makeSequenceTable(mergers)
dim(seq_table)
# 1,676 ASVs in 7 samples


#### 15) Remove chimeras ####
# chimeras are generally considered contaminants because they can be interpreted as a novel sequence when they are in fact an error
seq_table_nochim <- removeBimeraDenovo(seq_table, method = "consensus", multithread = TRUE, verbose = TRUE)
# found 19 chimeras out of 1,676 sequences (1% identified as chimeric. Great!)

# although we only lost 19 sequences, we don't know if they held a lot in terms of abundance. To find out, we can calculate the proportion of sequences retained after chimeras removed
sum(seq_table_nochim)/sum(seq_table)
# We retained 99% of sequence abundance. Great!

# ** Most of the reads should remain after chimera removal. If most of the reads were removed as chimeric, it likely requires revisiting the upstream processing (it usually means the primers were not properly removed).

# Inspect distribution of sequence lengths
table(nchar(getSequences(seq_table_nochim)))
# As expected, there is a quite a bit of length variability in the amplified ITS region.


#### 16) Track reads through the pipeline to verify everything worked as expected ####
getN <- function(x) sum(getUniques(x))
(summary_table <- data.frame(row.names = sample.names, 
                             input = out[, 1],
                             filtered = out[, 2], 
                             denoised_FWD = sapply(dada_FWD, getN),
                             denoised_REV = sapply(dada_REV, getN), 
                             merged = sapply(mergers, getN),
                             non_chim = rowSums(seq_table_nochim), 
                             final_perc_reads_retained = round(rowSums(seq_table_nochim)/out[, 1]*100, 1)))
# Looks good! The most amount of reads were removed during filtering, as expected. No large removal of reads in any sample after filtering. 
# Besides filtering, there should be no step in which a majority of reads are lost.


#### 17) Assign taxonomy with the UNITE ITS database ####
# DADA2 supports fungal taxonmic assignment using the UNITE database. For fungal taxonomy, the General Fasta release files from the UNITE ITS database can be downloaded and used as the reference. The file is called "UNITE_ref.10.10.2017.fasta". It is in the Microbiome analysis folder on the Drive.

unite.ref <- "/Users/kendallb/Documents/Research/Collaborations/Oak_fire_2019_prj/Oak_fire_2019_microbiome_data/UNITE_ref.10.10.2017.fasta" # set this to wherever the file is located on your computer.

taxa <- assignTaxonomy(seq_table_nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

# Inspecting the taxonomic assignments
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# giving seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seq_table_nochim)
asv_headers <- vector(dim(seq_table_nochim)[2], mode = "character")

for (i in 1:dim(seq_table_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# make ASV count table
asv_tab <- t(seq_table_nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "Oak_fire_2019_fungi_ASVs_counts_no_lowmod.tsv", sep = "\t", quote = F, col.names = NA)

# make taxa table
asv_taxa <- taxa
row.names(asv_taxa) <- sub(">", "", asv_headers)
write.table(asv_taxa, "Oak_fire_2019_fungi_ASVs_taxonomy_no_lowmod.tsv", sep = "\t", quote = F, col.names = NA)
