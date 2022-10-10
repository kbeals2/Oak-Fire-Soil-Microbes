
library(pacman)
p_load("ggplot2", "patchwork")
(.packages())

#### Hypothesis 2: Fire-induced changes to the soil microbiome will reduce success of oak seedlings ####

# Hyp 2a: ROOT BIOMASS BY SOIL BURN SEVERITY
burn_severity_labels <- c("Unburned", "Burned")

(root_soil_burn_fig <- ggplot(oak_blwgrd_rev, aes(x = factor(burn_severity, level = c('unburned', 'burned')), y = blwgrd_L, fill = burn_severity)) +
    geom_boxplot(width = 0.30, lwd = 1.3) +
    scale_fill_manual(values = c("#f58231", "#4363d8")) +
    scale_y_continuous(limits = c(0.6, 5.2), breaks = c(1, 2, 3, 4, 5)) +
    scale_x_discrete(labels = burn_severity_labels) +
    xlab("Soil burn severity") +
    ylab("Root biomass (g)") +
    theme_classic(base_size = 10, base_family = "") +
    theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2),
          plot.title = element_text(face = "bold", size = 22, hjust = -0.1, vjust = 0.2),
          axis.title = element_text(face = "bold", color = "black", size = 22),
          axis.text = element_text(face = "plain", color = "black", size = 20),
          axis.ticks.length = unit(0.15, "cm"),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          legend.position = "none"))

# Hyp 2b: ROOT BIOMASS BY PLANT PATHOGENIC FUNGI ABUNDANCE
(root_path_abund_fig <- ggplot(oak_root_fungi, aes(x = path_proportion, y = blwgrd_L)) + 
    geom_point(size = 3.75, shape = 21, fill = "white", color = "black", stroke = 1.5) +
    stat_smooth(method = "lm", color = "orange1") +
    xlab("Plant pathogen proportion (%)") +
    ylab("Root biomass (g)") +
    theme_classic(base_size = 10, base_family = "") +
    theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2),
          axis.title = element_text(face = "bold", color = "black", size = 20),
          axis.text = element_text(face = "plain", color = "black", size = 18),
          axis.ticks.length = unit(0.20, "cm"),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_blank()))

# Hyp 2b: ROOT BIOMASS BY PLANT PATHOGENIC FUNGI DIVERSITY (q1)
(root_pathogen_q1_fig <- ggplot(oak_root_fungi, aes(x = hill_q1_path, y = blwgrd_L)) + 
    geom_point(size = 3.75, shape = 21, fill = "white", color = "black", stroke = 1.5) +
    stat_smooth(method = "lm", color = "orange1") +
    xlab("Plant pathogen diversity (q1)") +
    ylab("Root biomass (g)") +
    theme_classic(base_size = 10, base_family = "") +
    theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2),
          axis.title = element_text(face = "bold", color = "black", size = 20),
          axis.text = element_text(face = "plain", color = "black", size = 18),
          axis.ticks.length = unit(0.20, "cm"),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_blank()))

# Hyp 2b: ROOT BIOMASS BY PLANT PATHOGEN DIVERSITY (q2)
(root_pathogen_q2_fig <- ggplot(oak_root_fungi, aes(x = hill_q2_path, y = blwgrd_L)) + 
    geom_point(size = 3.75, shape = 21, fill = "white", color = "black", stroke = 1.5) +
    stat_smooth(method = "lm", color = "orange1") +
    xlab("Plant pathogen diversity (q2)") +
    ylab("Root biomass (g)") +
    theme_classic(base_size = 10, base_family = "") +
    theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2),
          axis.title = element_text(face = "bold", color = "black", size = 20),
          axis.text = element_text(face = "plain", color = "black", size = 18),
          axis.ticks.length = unit(0.20, "cm"),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_blank()))

(patch <- (root_soil_burn_fig/root_path_abund_fig) | (root_pathogen_q1_fig/root_pathogen_q2_fig))



#### Hypothesis 3: Competition with pine seedlings reduces oak seedling success ####

# RELATIVE GROWTH RATE
rgr_label <- c(expression("Relative growth rate" ~ (cm^{2}*wk^{-1})))

(rgr_plant_neighbor_fig <- ggplot(oak_rgr, aes(x = factor(plant_community, level = c('quercus_quercus', 'quercus_pinus')), y = rgr_L)) +
    geom_boxplot(width = 0.5, lwd = 1.3) +
    ylab(rgr_label) +
    scale_y_continuous(limits = c(0, 0.049), breaks = c(0.01, 0.02, 0.03, 0.04)) +
    theme_classic(base_size = 10, base_family = "") +
    theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2),
          plot.title = element_text(face = "bold", size = 22, hjust = -0.1, vjust = 0.2),
          axis.title = element_blank(),
          axis.text.y = element_text(face = "plain", color = "black", size = 22),
          axis.text.x = element_blank(),
          axis.ticks.length = unit(0.15, "cm"),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5)))



# BELOWGROUND BIOMASS
(root_biomass_plant_neighbor_fig <- ggplot(oak_blwgrd, aes(x = factor(plant_community, level = c('quercus_quercus', 'quercus_pinus')), y = blwgrd_L)) +
    geom_boxplot(width = 0.5, lwd = 1.3) +
    ylab("Root biomass (g)") +
    scale_y_continuous(limits = c(0.5, 5), breaks = c(1, 2, 3, 4)) +
    theme_classic(base_size = 10, base_family = "") +
    theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2),
          plot.title = element_text(face = "bold", size = 22, hjust = -0.1, vjust = 0.2),
          axis.title = element_blank(),
          axis.text.y = element_text(face = "plain", color = "black", size = 22),
          axis.text.x = element_blank(),
          axis.ticks.length = unit(0.15, "cm"),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5)))


# Assembling figures using patchwork package
(patch <- rgr_plant_neighbor_fig /
    root_biomass_plant_neighbor_fig)

root_soil_burn_fig + patch


#### Hypothesis 4: The interaction of competition with pine seedlings and fire-induced changes to the soil microbiome influences oak seedling success ####

# BELOWGROUND BIOMASS
# remove control soil treatment
oak_blwgrd_no_control <- subset(oak_blwgrd, burn_severity != "none")

(interactive_effect_fig <- ggplot(oak_blwgrd_no_control, aes(x = factor(plant_community, level = c('quercus_quercus', 'quercus_pinus')), y = blwgrd_L, fill = factor(burn_severity, level = c('unburned', 'burned')))) +
    geom_boxplot(width = 0.5, lwd = 1.3) +
    scale_fill_manual(values = c("#4363d8", "#f58231"), name = "Burn severity", breaks = c("unburned", "burned"), labels = c("Unburned", "Burned")) +
    ylab("Root biomass (g)\n") +
    scale_y_continuous(limits = c(0.65, 4.85), breaks = c(1, 2, 3, 4)) +
    theme_classic(base_size = 10, base_family = "") +
    theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2),
          axis.title = element_blank(),
          axis.text.y = element_text(face = "plain", color = "black", size = 22),
          axis.text.x = element_blank(),
          axis.ticks.length = unit(0.15, "cm"),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5), 
          legend.position = "none"))



# Combining root biomass response to soil burn severity & microbial diversity by root biomass figures 
(patch1 <- root_soil_burn_fig + root_fungal_div_fig)
(patch2 <- Endo_prop_root_fig + Sapro_prop_root_fig)

(patch_one_and_two <- patch1 / patch2)


# Combining Hyp3 & Hyp4 figs
(patch3 <- rgr_plant_neighbor_fig + root_biomass_plant_neighbor_fig + interactive_effect_fig)


