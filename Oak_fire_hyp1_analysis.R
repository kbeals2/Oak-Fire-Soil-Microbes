### OAK FIRE 2019 PSF EXPERIMENT: ANALYSIS OF 16S & ITS2 AMPLICON SEQUENCES FROM FIELD SOIL (based off Mike Lee Happy Belly Bioinformatics tutorial: https://astrobiomike.github.io/amplicon/dada2_workflow_ex#the-data, .....) ###

## Project Questions
# Q1: Does fire affect the soil microbiome? 
# Hypothesis 1: Bacterial and fungal diversity is reduced with increasing burn severity.
# Q2: If fire affects the soil microbiome, are there particular microbial taxa associated with levels of burn severity?

# First, need to establish working directory (where R pulls files from).
getwd()
setwd("/Users/kendallb/Documents/Documents_KK_Macbook_Pro/Research/Collaborations/Oak_fire_2019_prj")

(.packages())
install.packages("pacman")
library(pacman)
p_load("fBasics", "rcompanion", "tidyverse", "lme4", "car", "vegan", "phyloseq", "DESeq2", "dendextend")

# First, test if dNDVI varies by unburned and burned sites
soil_site_metadata <- read.table("Oak_fire_site_metadata.csv", header = T, sep = ",")

ndvi_mod <- lmer(dNDVI ~ burn_status + (1|dd_lat) + (1|dd_long) + (1|elevation), data = soil_site_metadata)
Anova(ndvi_mod)
# burn status: Chi sq = 36.3, p < 0.0001
# Great! dNDVI is significantly different between unburned & burned sites

(t.test(dNDVI ~ burn_severity, data = soil_site_metadata))
# t = 7.0841, p = 0.001987

# Then test if soil pH varies by unburned and burned sites
ph_mod <- lmer(soil_pH ~ burn_status + (1|dd_lat) + (1|dd_long) + (1|elevation), data = soil_site_metadata)
Anova(ph_mod)
# burn status: Chi sq = 26.03, p < 0.0001
# pH is significantly greater in unburned sites 

(t.test(soil_pH ~ burn_severity, data = soil_site_metadata))
# t = -6.4789, p = 0.0179

setwd("/Users/kendallb/Documents/Documents_KK_Macbook_Pro/Research/Collaborations/Oak_fire_2019_prj/Oak_fire_2019_microbiome_data")

#### 1) LOAD IN DATA #### 
bact_count_table <- read.table("Oak_fire_Bact_ASVs_counts.tsv", header = T, row.names = 1, check.names = F, sep = "\t")
fungi_count_table <- read.table("Oak_fire_Fungi_ASVs_counts.tsv", header = T, row.names = 1, check.names = F, sep = "\t")

bact_taxa_table <- read.table("Oak_fire_Bact_ASVs_taxonomy.tsv", header = T, check.names = F, sep = "\t")
fungi_taxa_table <- read.table("Oak_fire_Fungi_ASVs_taxonomy.tsv", header = T, check.names = F, sep = "\t")

sample_info_table <- read.table("Oak_fire_sample_info.csv", header = T, row.names = 1, sep = ",")
# Note: since the row names I want are in column 1, we can set the row names to be those that are in column 1.


#### 2) NORMALIZE FOR SAMPLING DEPTH USING A VARIANCE STABILIZING TRANSFORMATION RATHER THAN RAREFYING BY LOWEST SEQUENCING DEPTH ####
bact_deseq_counts <- DESeqDataSetFromMatrix(bact_count_table, colData = sample_info_table, design = ~burn_status) 
bact_deseq_counts_vst <- varianceStabilizingTransformation(bact_deseq_counts)

fungi_deseq_counts <- DESeqDataSetFromMatrix(fungi_count_table, colData = sample_info_table, design = ~burn_status) 
fungi_deseq_counts_vst <- varianceStabilizingTransformation(fungi_deseq_counts)
# If receive error "every gene contains at least one zero, cannot compute log geometric means," the count table is sparse with many zeroes. In this case, use "deseq_counts <- estimateSizeFactors(deseq_counts)", then continue to varianceStabilizingTransformation step.

fungi_deseq_counts <- estimateSizeFactors(fungi_deseq_counts, type = "poscounts")
fungi_deseq_counts_vst <- varianceStabilizingTransformation(fungi_deseq_counts)

bact_vst_trans_count_table <- assay(bact_deseq_counts_vst) # extracting transformed data table
fungi_vst_trans_count_table <- assay(fungi_deseq_counts_vst) 

# Because some ASVs had a count of 0, they became negative after variance stabilizing transformation. So now we need to change all ASV counts that are less than 0 (i.e. negative) to 0. 
bact_vst_trans_count_table[bact_vst_trans_count_table < 0.0] <- 0.0
fungi_vst_trans_count_table[fungi_vst_trans_count_table < 0.0] <- 0.0

# Bray-Curtis distance doesn't like working with zeros. To fix this, add small number to cells with zero. This step would be necessary even if we were not using Bray-Curtis distance in order to comply with the capscale function assumptions needed to perform db-RDA.
bact_vst_trans_count_table_001 <- (bact_vst_trans_count_table + 0.001)
fungi_vst_trans_count_table_001 <- (fungi_vst_trans_count_table + 0.001)


#### 3) RICHNESS & DIVERSITY ESTIMATES (addresses Q1: Does fire affect the soil microbiome?) ####
# Use hilldiv package to estimate diversity with effective number of ASVs; order of diversity (q) specifies sensitivity towards abundant and rare ASVs.

# Based off of: https://github.com/anttonalberdi/hilldiv 
install.packages("hilldiv")
library(hilldiv)

# Hilldiv requires ASV table in which each column is a sample and each row is an ASV. Data can be raw counts or relative proportions.
# q = 0; raw richness; only considers whether an ASV is present or absent; same weight is given to all ASVs regardless of abundance (weighs rare taxa same as abundant taxa)
hill_div(bact_count_table, qvalue = 0)
hill_div(fungi_count_table, qvalue = 0)

# q = 1; weighs ASVs by their frequency without disproportionately favoring either rare or abundant ones
hill_div(bact_count_table, qvalue = 1)
hill_div(fungi_count_table, qvalue = 1)

# q = 2; overweighs abundant ASVs
hill_div(bact_count_table, qvalue = 2)
hill_div(fungi_count_table, qvalue = 2)

# create hierarchy tables that has burn status assigned to each sample
sample <- c("SCEA-BGB", "SCEA-PL441", "SCEA-RS441", "SCEA-RGB", "SCEA-FCM", "SCEA-LCM", "SCEA-SCM")
burn_status <- c(rep("unburned", 3), rep("burned", 4))
hierarchy_table <- cbind(sample, burn_status)
colnames(hierarchy_table) <- c("Sample", "Burn_status")

# test whether diversity values differ by burn severity
(div_test(bact_count_table, qvalue = 0, hierarchy = hierarchy_table))
# t = 2.26, p = 0.073
(div_test(bact_count_table, qvalue = 1, hierarchy = hierarchy_table))
# t = 1.42, p = 0.22
(div_test(bact_count_table, qvalue = 2, hierarchy = hierarchy_table))
# t = 0.35, p = 0.75

(div_test(fungi_count_table, qvalue = 0, hierarchy = hierarchy_table))
# t = -5.79, p = 0.011 *
# On average, unburned sites had 2.7x greater fungal diversity than burned sites
(div_test(fungi_count_table, qvalue = 1, hierarchy = hierarchy_table))
# t = -3.08, p = 0.077
(div_test(fungi_count_table, qvalue = 2, hierarchy = hierarchy_table))
# t = -2.99, p = 0.07


#### 4) DOES THE AMOUNT OF EACH FUNGAL GUILD DIFFER BETWEEN UNBURNED & BURNED SOIL? ####

# load in count table with row ID as first column and fungal_guild column added
funguild_count_table <- read.table("Oak_fire_FunGuild_count_for_R.csv", header = T, row.names = 1, sep = ",")

# Sidebar: Identify proportions of each fungal guild by counting number of rows for each guild
funguild_count_table %>% group_by(fungal_guild) %>% tally()
# Saprotrophic: 65% (401 ASVs)
# Other: 11.5 % (71 ASVs)
# ECM: 6.7% (41 ASVs)
# AMF: 6.2% (38 ASVs)
# Plant pathogen: 6% (37 ASVs)
# Endophyte: 4.1% (25 ASVs)
# Ericoid: 0.5% (3 ASVs)

# melt table so that count is its own column
library(reshape)
funguild_count_melted <- melt(funguild_count_table, id = c("ASV_ID", "fungal_guild"))

# rename columns (variable = site, value = count) (note: rename function does not work when reshape packages is installed; uninstall reshape before running rename function)
funguild_count_melted <-rename(funguild_count_melted,
                               "variable" = "site",
                               "value" = "count")

# add column that assigns burn category based on condition in site column
funguild_count_melted <- funguild_count_melted %>% mutate(burn_category = case_when(
  site == "BGB" ~ "Unburned",
  site == "PL441" ~ "Unburned",
  site == "RS441" ~ "Unburned",
  site == "RGB" ~ "Burned",
  site == "FCM" ~ "Burned",
  site == "LCM" ~ "Burned",
  site == "SCM" ~ "Burned"))

# add column that assigns sample_count_total based on condition in site column
funguild_count_melted <- funguild_count_melted %>% mutate(sample_count_total = case_when(
  site == "BGB" ~ 7063,
  site == "PL441" ~ 7528,
  site == "RS441" ~ 7177,
  site == "RGB" ~ 3001,
  site == "FCM" ~ 7950,
  site == "LCM" ~ 9591,
  site == "SCM" ~ 20073))

# create smaller dataframes by subsetting by each fungal guild
AMF <- funguild_count_melted %>% filter(fungal_guild == "Arbuscular_mycorrhizal")
ECM <- funguild_count_melted %>% filter(fungal_guild == "Ectomycorrhizal")
Endophyte <- funguild_count_melted %>% filter(fungal_guild == "Endophyte")
Plant_path <- funguild_count_melted %>% filter(fungal_guild == "Plant_pathogen")
Saprotrophs <- funguild_count_melted %>% filter(fungal_guild == "Saprotrophic")

## AMF ##
AMF_abundance_tab <- AMF %>% group_by(burn_category, site) %>% tally(count) # create small table that has the count sum for each site and grouped by burn category
site_total_count <- c(3001, 7950, 9591, 20073, 7063, 7528, 7177) # create vector of total count for each site
AMF_abundance_tab$site_total_count <- site_total_count # create new column in table and assign as vector made above
AMF_abundance_tab$rel_abundance <- AMF_abundance_tab$n/AMF_abundance_tab$site_total_count # create new column of relative abundance

plotNormalHistogram(AMF_abundance_tab$rel_abundance)
normalTest(AMF_abundance_tab$rel_abundance, method = c("sw")) 

AMF_abundance_mod <- lm(rel_abundance ~ burn_category, data = AMF_abundance_tab)
Anova(AMF_abundance_mod)
# F = 4.7354, p = 0.08151

# Average relative abundance by burn severity
(AMF_abund_stats <- AMF_abundance_tab %>%
    group_by(burn_category) %>%
    summarize(N = length(rel_abundance),
              Mean = mean(rel_abundance),
              SD = sd(rel_abundance),
              SE = SD/sqrt(N),
              Lower_ci = Mean - SE,
              Upper_ci = Mean + SE))
# 2.98x more AMF in unburned soil than burned soil


## ECM ##
ECM_abundance_tab <- ECM %>% group_by(burn_category, site) %>% tally(count)
ECM_abundance_tab$site_total_count <- site_total_count 
ECM_abundance_tab$rel_abundance <- ECM_abundance_tab$n/ECM_abundance_tab$site_total_count

plotNormalHistogram(ECM_abundance_tab$rel_abundance)
normalTest(ECM_abundance_tab$rel_abundance, method = c("sw")) # need to normalize

plotNormalHistogram(logit(ECM_abundance_tab$rel_abundance))
normalTest(logit(ECM_abundance_tab$rel_abundance), method = c("sw")) # recommended to use logit transform rather than arcsine

ECM_abundance_tab$rel_abund_logit <- logit(ECM_abundance_tab$rel_abundance)

ECM_abundance_mod <- lm(rel_abund_logit ~ burn_category, data = ECM_abundance_tab)
Anova(ECM_abundance_mod)
# F = 1.0337, p = 0.3559


## Endophyte ##
Endophyte_abundance_tab <- Endophyte %>% group_by(burn_category, site) %>% tally(count)
Endophyte_abundance_tab$site_total_count <- site_total_count 
Endophyte_abundance_tab$rel_abundance <- Endophyte_abundance_tab$n/Endophyte_abundance_tab$site_total_count

plotNormalHistogram(Endophyte_abundance_tab$rel_abundance)
normalTest(Endophyte_abundance_tab$rel_abundance, method = c("sw"))

Endophyte_abundance_mod <- lm(rel_abundance ~ burn_category, data = Endophyte_abundance_tab)
Anova(Endophyte_abundance_mod)
# F = 0.0787, p = 0.7903


## Plant pathogen ##
Plant_path_abundance_tab <- Plant_path %>% group_by(burn_category, site) %>% tally(count)
Plant_path_abundance_tab$site_total_count <- site_total_count 
Plant_path_abundance_tab$rel_abundance <- Plant_path_abundance_tab$n/Plant_path_abundance_tab$site_total_count
Plant_path_abundance_tab$proportion <- Plant_path_abundance_tab$rel_abundance*100

plotNormalHistogram(Plant_path_abundance_tab$rel_abundance)
normalTest(Plant_path_abundance_tab$rel_abundance, method = c("sw"))

Plant_path_abundance_mod <- lm(rel_abundance ~ burn_category, data = Plant_path_abundance_tab)
Anova(Plant_path_abundance_mod)
# F = 10.966, p = 0.0212 *

# Average relative abundance by burn status
(plant_path_abund_stats <- Plant_path_abundance_tab %>%
    group_by(burn_category) %>%
    summarize(N = length(proportion),
              Mean = mean(proportion),
              SD = sd(proportion),
              SE = SD/sqrt(N),
              Lower_ci = Mean - SE,
              Upper_ci = Mean + SE))
# 2.54x more plant pathogenic fungi in unburned soil than in burned soil


## Saprotroph ##
Sapro_abundance_tab <- Saprotrophs %>% group_by(burn_category, site) %>% tally(count)
Sapro_abundance_tab$site_total_count <- site_total_count 
Sapro_abundance_tab$rel_abundance <- Sapro_abundance_tab$n/Sapro_abundance_tab$site_total_count

plotNormalHistogram(Sapro_abundance_tab$rel_abundance)
normalTest(Sapro_abundance_tab$rel_abundance, method = c("sw"))

Sapro_abundance_mod <- lm(rel_abundance ~ burn_category, data = Sapro_abundance_tab)
Anova(Sapro_abundance_mod)
# F = 0.0348, p = 0.8593


#### 5) DIVERSITY ESTIMATES OF EACH FUNGAL GUILD BETWEEN UNBURNED & BURNED SOIL ####

library(hilldiv)

# load in count table with row ID as first column and fungal_guild column added
funguild_count_table <- read.table("Oak_fire_FunGuild_count_for_R.csv", header = T, row.names = 1, sep = ",")

# make separate dataframes for each fungal guild
AMF_count_tab <- funguild_count_table %>% filter(fungal_guild == "Arbuscular_mycorrhizal")
AMF_count_tab_new <- AMF_count_tab %>% remove_rownames() %>% column_to_rownames(var = "ASV_ID")
AMF_count_tab_new <- AMF_count_tab_new[-c(8)]

ECM_count_tab <- funguild_count_table %>% filter(fungal_guild == "Ectomycorrhizal")
ECM_count_tab_new <- ECM_count_tab %>% remove_rownames() %>% column_to_rownames(var = "ASV_ID")
ECM_count_tab_new <- ECM_count_tab_new[-c(8)]

Endophyte_count_tab <- funguild_count_table %>% filter(fungal_guild == "Endophyte")
Endophyte_count_tab_new <- Endophyte_count_tab %>% remove_rownames() %>% column_to_rownames(var = "ASV_ID")
Endophyte_count_tab_new <- Endophyte_count_tab_new[-c(8)]

Plant_path_count_tab <- funguild_count_table %>% filter(fungal_guild == "Plant_pathogen")
Plant_path_count_tab_new <- Plant_path_count_tab %>% remove_rownames() %>% column_to_rownames(var = "ASV_ID")
Plant_path_count_tab_new <- Plant_path_count_tab_new[-c(8)]

Saprotroph_count_tab <- funguild_count_table %>% filter(fungal_guild == "Saprotrophic")
Saprotroph_count_tab_new <- Saprotroph_count_tab %>% remove_rownames() %>% column_to_rownames(var = "ASV_ID")
Saprotroph_count_tab_new <- Saprotroph_count_tab_new[-c(8)]

# create hierarchy tables that has burn severity assigned to each sample
site <- c("BGB", "PL441", "RS441", "RGB", "FCM", "LCM", "SCM")
burn_status <- c(rep("unburned", 3), rep("burned", 4))
hierarchy_table <- cbind(site, burn_status)
colnames(hierarchy_table) <- c("Site", "Burn_status")

# q = 0; raw richness; only considers whether an ASV is present or absent; same weight is given to all ASVs regardless of abundance (weighs rare taxa same as abundant taxa)
hill_div(AMF_count_tab_new, qvalue = 0)
hill_div(ECM_count_tab_new, qvalue = 0)
hill_div(Endophyte_count_tab_new, qvalue = 0)
hill_div(Plant_path_count_tab_new, qvalue = 0)
hill_div(Saprotroph_count_tab_new, qvalue = 0)

# test whether diversity values differ by burn status
(div_test(AMF_count_tab_new, qvalue = 0, hierarchy = hierarchy_table))
# t = -1.50, p = 0.20
(div_test(ECM_count_tab_new, qvalue = 0, hierarchy = hierarchy_table))
# t = 0.77, p = 0.48
(div_test(Endophyte_count_tab_new, qvalue = 0, hierarchy = hierarchy_table))
# t = -2.40, p = 0.063 
(div_test(Plant_path_count_tab_new, qvalue = 0, hierarchy = hierarchy_table))
# t = -5.47, p = 0.003 **
(div_test(Saprotroph_count_tab_new, qvalue = 0, hierarchy = hierarchy_table))
# t = -3.28, p = 0.048 **

hill_div(AMF_count_tab_new, qvalue = 1)
hill_div(ECM_count_tab_new, qvalue = 1)
hill_div(Endophyte_count_tab_new, qvalue = 1)
hill_div(Plant_path_count_tab_new, qvalue = 1)
hill_div(Saprotroph_count_tab_new, qvalue = 1)

# test whether diversity values differ by burn status
(div_test(AMF_count_tab_new, qvalue = 1, hierarchy = hierarchy_table))
# t = -1.34, p = 0.24
(div_test(ECM_count_tab_new, qvalue = 1, hierarchy = hierarchy_table))
# t = 1.24, p = 0.30
(div_test(Endophyte_count_tab_new, qvalue = 1, hierarchy = hierarchy_table))
# t = -1.08, p = 0.34
(div_test(Plant_path_count_tab_new, qvalue = 1, hierarchy = hierarchy_table))
# t = -7.83, p = 0.001 **
(div_test(Saprotroph_count_tab_new, qvalue = 1, hierarchy = hierarchy_table))
# t = -2.50, p = 0.077 


hill_div(AMF_count_tab_new, qvalue = 2)
hill_div(ECM_count_tab_new, qvalue = 2)
hill_div(Endophyte_count_tab_new, qvalue = 2)
hill_div(Plant_path_count_tab_new, qvalue = 2)
hill_div(Saprotroph_count_tab_new, qvalue = 2)

# test whether diversity values differ by burn status
(div_test(AMF_count_tab_new, qvalue = 2, hierarchy = hierarchy_table))
# t = -0.99, p = 0.39
(div_test(ECM_count_tab_new, qvalue = 2, hierarchy = hierarchy_table))
# t = 1.32, p = 0.25
(div_test(Endophyte_count_tab_new, qvalue = 2, hierarchy = hierarchy_table))
# t = -0.95, p = 0.40
(div_test(Plant_path_count_tab_new, qvalue = 2, hierarchy = hierarchy_table))
# t = -6.48, p = 0.002 **
(div_test(Saprotroph_count_tab_new, qvalue = 2, hierarchy = hierarchy_table))
# t = -1.50, p = 0.22
