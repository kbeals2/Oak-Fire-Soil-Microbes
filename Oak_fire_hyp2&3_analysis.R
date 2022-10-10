getwd()
setwd("/Users/kendallb/Documents/Documents_KK_Macbook_Pro/Research/Collaborations/Oak_fire_2019_prj") # change this to where your files are stored

# check what packages are installed and loaded
(.packages())

# if no packages are loaded:
install.packages("pacman")
library(pacman)
p_load("fBasics", "rcompanion", "tidyverse", "lme4", "lmerTest", "car",  "agricolae", "bbmle")

#### Load in data ####
master_rev_data <- read.csv("oak_pine_master_revised_r.csv") # does not contain LGB or NEON sites ("low-mod" burn)

# Check structure of data sets
str(master_data) 

# First, need to subset by just pine-oak & oak-oak pots.
oak_rev_data <- master_rev_data %>% filter(plant_community != "pinus_pinus")

# remove all rows that are control pots (i.e. where burn_severity is labelled as "none")
oak_rev_no_control <- subset(oak_rev_data, burn_status != "none")


#### Hypothesis 2a: Fire-induced changes to the soil microbiome will reduce oak seedling growth ####

## ABOVEGROUND BIOMASS
# create smaller dataframe that just contains aboveground_L as only data column
oak_rev_abvgrd <- oak_rev_no_control %>% select(pot_id, plant_community, plant_community_id, burn_status, burn_status_id, site, site_id, abvgrd_L)

# Before continuing, remove any NAs from data column.
oak_rev_abvgrd <- oak_rev_abvgrd %>% drop_na(abvgrd_L)

# Check distribution of abvgrd biomass.
plotNormalHistogram(oak_rev_abvgrd$abvgrd_L)
normalTest(oak_rev_abvgrd$abvgrd_L, method = c("sw"))
# If data follow a normal distribution (p value > 0.05), then there is no need to transform data. Proceed to building and running the model. If data are NOT normally distributed, evaluate distributions from transformations below. Pick the transformation that gives the highest p value. While not listed below, other transformations you can try are square (^2), cube (^3), logarithmic (log), and reciprocal (1/).

normalTest(transformTukey(oak_rev_abvgrd$abvgrd_L), method = c("sw"))

# Now we need to create an aboveground_L column of Tukey transformed data that we will use in the model in the next step.
oak_rev_abvgrd$abvgrd_L_tukey <- transformTukey(oak_rev_abvgrd$abvgrd_L)

# Now build and run the model. We can actually run a mixed effects model here using plant community and site as random effects. Even though we're not fundamentally interested in how the soil of individual field sites affects plant growth, it is likely that field sites differ somewhat in their soil composition which could influence plant growth. So we want to be able to account for that variation contributed to sites.
oak_rev_abvgrd_hyp2_mod <- lmer(abvgrd_L_tukey ~ burn_status + (1|plant_community) + (1|site), data = oak_rev_abvgrd)

# View the output of the model. Note that instead of a F statistic, a Chi square value is reported. Chi square is essentially the same thing as the F statistic (It's still a ratio of the between group variance to within group variance). It will probably be useful to write a note with the results for the model. Enter the Chi square values and p values into the Results Excel file.
Anova(oak_rev_abvgrd_hyp2_mod)
# burn severity: Chi sq = 0.39, p = 0.53

## RELATIVE GROWTH RATE
oak_rgr_rev <- oak_rev_no_control %>% select(pot_id, plant_community, plant_community_id, burn_status, burn_status_id, site, site_id, rgr_L)
oak_rgr_rev <- oak_rgr_rev %>% filter(rgr_L > 0)

plotNormalHistogram(oak_rgr_rev$rgr_L)
normalTest(oak_rgr_rev$rgr_L, method = c("sw"))

normalTest(transformTukey(oak_rgr_rev$rgr_L), method = c("sw"))

oak_rgr_rev$rgr_L_tukey <- transformTukey(oak_rgr_rev$rgr_L)

oak_rev_rgr_hyp2_mod <- lmer(rgr_L_tukey ~ burn_status + (1|plant_community) + (1|site), data = oak_rgr_rev)
Anova(oak_rev_rgr_hyp2_mod)
# burn severity: Chi sq = 0.24, p = 0.62


## BELOWGROUND BIOMASS
# create smaller dataframe that just contains belowground_L as only data column
oak_blwgrd_rev <- oak_rev_no_control %>% select(pot_id, plant_community, plant_community_id, burn_status, burn_status_id, site, site_id, blwgrd_L)

oak_blwgrd_rev <- oak_blwgrd_rev %>% drop_na(blwgrd_L)

plotNormalHistogram(oak_blwgrd$blwgrd_L)
normalTest(oak_blwgrd_rev$blwgrd_L, method = c("sw"))

oak_rev_blwgrd_hyp2_mod <- lmer(blwgrd_L ~ burn_status + (1|plant_community) + (1|site), data = oak_blwgrd_rev)
Anova(oak_rev_blwgrd_hyp2_mod)
# burn severity: Chi sq = 3.66, p = 0.056

# Average values of root bimass L individuals by burn severity
(oak_blwgrd_stats <- ddply(oak_blwgrd_rev, c("burn_status"),
                           summarise, 
                           N = length(blwgrd_L),
                           mean = mean(blwgrd_L),
                           sd = sd(blwgrd_L),
                           se = sd/sqrt(N),
                           lower_ci = mean - se,
                           upper_ci = mean + se))
# 24.7% more root biomass in burned soil inoculum than unburned soil inoculum


## SPECIFIC ROOT LENGTH
# create smaller dataframe that just contains srl_L as only data column
oaks_srl_rev <- oak_rev_no_control %>% select(pot_id, plant_community, plant_community_id, burn_status, burn_status_id, site, site_id, srl_L)

oaks_srl_rev <- oaks_srl_rev %>% drop_na(srl_L)

plotNormalHistogram(oaks_srl$srl_L)
normalTest(oaks_srl_rev$srl_L, method = c("sw"))

normalTest(transformTukey(oaks_srl_rev$srl_L), method = c("sw"))

oaks_srl_rev$srl_L_tukey <- transformTukey(oaks_srl_rev$srl_L)

oaks_rev_srl_hyp2_mod <- lmer(srl_L_tukey ~ burn_status + (1|plant_community) + (1|site), data = oaks_srl_rev)
Anova(oaks_rev_srl_hyp2_mod)
# burn severity: Chi sq = 0.32, p = 0.57


## SPECIFIC LEAF AREA
# create smaller dataframe that just contains srl_L as only data column
oak_sla_rev <- oak_rev_no_control %>% select(pot_id, plant_community, plant_community_id, burn_status, burn_status_id, site, site_id, sla_L)

oak_sla_rev <- oak_sla_rev %>% drop_na(sla_L)

plotNormalHistogram(oak_sla$sla_L)
normalTest(oak_sla_rev$sla_L, method = c("sw"))

oak_rev_sla_hyp2_mod <- lmer(sla_L ~ burn_status + (1|plant_community) + (1|site), data = oak_sla_rev)
Anova(oak_rev_sla_hyp2_mod)
# burn severity: Chi sq = 1.46, p = 0.23


## CHLOROPHYLL CONTENT
# create smaller dataframe that just contains srl_L as only data column
oak_chl_rev <- oak_rev_no_control %>% select(pot_id, plant_community, plant_community_id, burn_status, burn_status_id, site, site_id, Chl_L)

plotNormalHistogram(oak_chl$Chl_L)
normalTest(oak_chl_rev$Chl_L, method = c("sw"))

normalTest(transformTukey(oak_chl_rev$Chl_L), method = c("sw"))

oak_chl_rev$Chl_L_tukey <- transformTukey(oak_chl_rev$Chl_L)

oak_rev_chl_hyp2_mod <- lmer(Chl_L_tukey ~ burn_status + (1|plant_community) + (1|site), data = oak_chl_rev)
Anova(oak_rev_chl_hyp2_mod)
# burn severity: Chi sq = 1.62, p = 0.20


#### Hypothesis 2b: Fire-induced changes to fungal diversity and abundance of specific fungal guilds will alter oak seedling root biomass ####

# Assessing the effect of microbial diversity and fungal guild relative abundance 
master_rev_fungi <- read.csv("oak_pine_revised_with_fungi_hill_numbers_and_guilds.csv")

master_rev_fungi <- master_rev_fungi %>% filter(plant_community != "pinus_pinus")

master_rev_fungi <- subset(master_rev_fungi, burn_status != "none")

oak_root_fungi <- master_rev_fungi %>% select(pot_id, plant_community, plant_community_id, burn_status, burn_status_id, site, site_id, hill_q0_fungi, hill_q0_path, hill_q0_sap, hill_q1_path, hill_q2_path, path_proportion, blwgrd_L)

oak_root_fungi <- oak_root_fungi %>% drop_na(blwgrd_L)

plotNormalHistogram(oak_root_fungi$blwgrd_L)
normalTest(oak_root_fungi$blwgrd_L, method = c("sw"))

# Effect of difference in overall fungal diversity caused by burn  
oak_root_hyp2_fungi_q0_mod <- lmer(blwgrd_L ~ hill_q0_fungi + (1|plant_community) + (1|site), data = oak_root_fungi)
Anova(oak_root_hyp2_fungi_q0_mod)
# hill_q0_fungi: Chi sq = 1.79, p = 0.18

# To get coefficient of determination from mixed effects model, use r.squaredGLMM function in MuMIn package
install.packages("MuMIn")
library(MuMIn)

r.squaredGLMM(oak_root_hyp2_fungi_q0_mod)
# R sq.m = 0.06; R sq.c = 0.37
# R sq.m: variance explained by fixed effect (i.e. fungal diversity)
# R sq.c: variance explained by entire model including random effect(s)


# Effect of difference in plant pathogen diversity caused by burn
oak_root_hyp2_path_q0_mod <- lmer(blwgrd_L ~ hill_q0_path + (1|plant_community) + (1|site), data = oak_root_fungi)
Anova(oak_root_hyp2_path_q0_mod)
# hill_q0_path: Chi sq = 1.46, p = 0.23

oak_root_hyp2_path_q1_mod <- lmer(blwgrd_L ~ hill_q1_path + (1|plant_community) + (1|site), data = oak_root_fungi)
Anova(oak_root_hyp2_path_q1_mod)
# hill_q0_path: Chi sq = 3.95, p = 0.047 *

r.squaredGLMM(oak_root_hyp2_path_q1_mod)
# R sq.m = 0.106; R sq.c = 0.362
# R sq.m: variance explained by fixed effect (i.e. q1 diversity)
# R sq.c: variance explained by entire model including random effect(s)

oak_root_hyp2_path_q2_mod <- lmer(blwgrd_L ~ hill_q2_path + (1|plant_community) + (1|site), data = oak_root_fungi)
Anova(oak_root_hyp2_path_q2_mod)
# hill_q0_path: Chi sq = 5.64, p = 0.018 *

r.squaredGLMM(oak_root_hyp2_path_q2_mod)
# R sq.m = 0.129; R sq.c = 0.359
# R sq.m: variance explained by fixed effect (i.e. q1 diversity)
# R sq.c: variance explained by entire model including random effect(s)


# Effect of difference in saprotroph diversity caused by burn
oak_root_hyp2_sap_q0_mod <- lmer(blwgrd_L ~ hill_q0_sap + (1|plant_community) + (1|site), data = oak_root_fungi)
Anova(oak_root_hyp2_sap_q0_mod)
# hill_q0_path: Chi sq = 0.68, p = 0.41


# Effect of difference in relative abundance of plant pathogenic fungi
oak_root_hyp2_path_abund_mod <- lmer(blwgrd_L ~ path_proportion + (1|plant_community) + (1|site), data = oak_root_fungi)
Anova(oak_root_hyp2_path_abund_mod)
# path proportion: Chi sq = 0.5166, p = 0.4723


#### Hypothesis 3: The interaction of competition with pine seedlings and fire-induced changes to the soil microbiome influences oak seedling success ####

## ABOVEGROUND BIOMASS
oak_abvgrd_hyp4_mod <- lmer(abvgrd_L_tukey ~ plant_community + burn_status + plant_community*burn_status + (1|site), data = oak_abvgrd)
Anova(oak_abvgrd_hyp4_mod)
# plant community: Chi sq = 2.58, p = 0.11
# burn severity: Chi sq = 0.39, p = 0.53
# plant community x burn severity: Chi sq = 0.63, p = 0.43


## RELATIVE GROWTH RATE
oak_rgr_hyp4_mod <- lmer(rgr_L_tukey ~ plant_community + burn_status + plant_community*burn_status + (1|site), data = oak_rgr)
Anova(oak_rgr_hyp4_mod)
# plant community: Chi sq = 5.63, p = 0.018
# burn severity: Chi sq = 0.25, p = 0.61
# plant community x burn severity: Chi sq = 1.27, p = 0.26


## BELOWGROUND BIOMASS
oak_blwgrd_hyp4_mod <- lmer(blwgrd_L ~ plant_community + burn_status + plant_community*burn_status + (1|site), data = oak_blwgrd)
Anova(oak_blwgrd_hyp4_mod)
# plant community: Chi sq = 7.06, p = 0.0079
# burn severity: Chi sq = 4.1, p = 0.042
# plant community x burn severity: Chi sq = 3.43, p = 0.064

oak_blwgrd_hyp4_posthoc_mod <- aov(blwgrd_L ~plant_community + burn_status + plant_community*burn_status + (1|site_id), data = oak_blwgrd)
TukeyHSD(oak_blwgrd_hyp4_posthoc_mod)
# oak-oak/burned soil vs. oak-pine/burned soil: p = 0.013
# oak-pine/unburned soil vs. oak-pine/burned soil: p = 0.019

# Average values of belowground biomass L individuals by plant neighbor and soil treatment
(oak_blwgrd_hyp4_stats <- ddply(oak_blwgrd, c("plant_community", "burn_status"),
                                summarise, 
                                N = length(blwgrd_L),
                                mean = mean(blwgrd_L),
                                sd = sd(blwgrd_L),
                                se = sd/sqrt(N),
                                lower_ci = mean - se,
                                upper_ci = mean + se))
# oaks produced 39% more roots with burned soil microbiome when growing with pines relative to growing with other oaks (i.e. soil microbiome only has effect when oaks are growing with pines)

# oaks produced 39% more roots with pines when interacting with burned soil microbiome relative to unburned soil microbiome (i.e. advantage with pines only present when growing with burned soil microbiome)


## SPECIFIC ROOT LENGTH
oak_srl_hyp4_mod <- lmer(srl_L_tukey ~ plant_community + burn_status + plant_community*burn_status + (1|site), data = oak_srl)
Anova(oak_srl_hyp4_mod)
# plant community: Chi sq = 0.91, p = 0.34
# burn severity: Chi sq = 0.28, p = 0.60
# plant community x burn severity: Chi sq = 0.31, p = 0.58


## SPECIFIC LEAF AREA
oak_sla <- oak_sla %>% drop_na(sla_L)

oak_sla_hyp4_mod <- lmer(sla_L ~ plant_community + burn_status + plant_community*burn_status + (1|site), data = oak_sla)
Anova(oak_sla_hyp4_mod)
# plant community: Chi sq = 0.37, p = 0.55
# burn severity: Chi sq = 1.36, p = 0.24
# plant community x burn severity: Chi sq = 2.05, p = 0.15


## CHLOROPHYLL CONTENT
oak_chl_hyp4_mod <- lmer(Chl_L_tukey ~ plant_community + burn_status + plant_community*burn_status + (1|site), data = oak_chl)
Anova(oak_chl_hyp4_mod)
# plant community: Chi sq = 0.05, p = 0.82
# burn severity: Chi sq = 2.11, p = 0.15
# plant community x burn severity: Chi sq = 0.006, p = 0.94



