
install.packages("pacman")
library(pacman)

p_load("ggplot2", "patchwork", "devtools")
(.packages())


#### MICROBIAL DIVERSITY & FUNGAL GUILD RELATIVE ABUNDANCE ####

## Bacteria ##
bact_q0_values <- c(263, 369, 607, 1022, 559, 566, 1044)
bact_q1_values <- c(168.8, 277.9, 433.5, 596.6, 279.2, 352.3, 518.1)
bact_q2_values <- c(118.8, 212.2, 276.1, 317.8, 134.8, 215.5, 223)
sample_names <- c("BGB", "PL441", "RS441", "RGB", "FCM", "LCM", "SCM")
burn_severity <- c(rep("unburned", 3), rep("burned", 4))
bact_q_table <- cbind(sample_names, burn_severity, bact_q0_values, bact_q1_values, bact_q2_values)
colnames(bact_q_table) <- c("sample", "burn_severity", "q0", "q1", "q2")

# bacteria q0 unburned mean
mean(c(263, 369, 607))
# 413

# bacteria q0 unburned sd
sd(c(263, 369, 607))
# 176.17

# bacteria q1 unburned mean
mean(c(168.8, 277.9, 433.5))
# 293.4

# bacteria q1 unburned sd
sd(c(168.8, 277.9, 433.5))
# 133.03

# bacteria q2 unburned mean
mean(c(118.8, 212.2, 276.1))
# 202.37

# bacteria q2 unburned sd
sd(c(118.8, 212.2, 276.1))
# 79.11

# bacteria q0 burned mean
mean(c(1022, 559, 566, 1044))
# 797.75

# bacteria q0 burned sd
sd(c(1022, 559, 566, 1044))
# 271.81

# bacteria q1 burned mean
mean(c(596.6, 279.2, 352.3, 518.1))
# 436.55

# bacteria q1 burned sd
sd(c(596.6, 279.2, 352.3, 518.1))
# 146.20

# bacteria q2 burned mean
mean(c(317.8, 134.8, 215.5, 223))
# 222.775

# bacteria q2 burned sd
sd(c(317.8, 134.8, 215.5, 223))
# 74.88

bact_mean <- c(413, 293.4, 202.37, 797.75, 436.55, 222.775)
bact_sd <- c(176.17, 133.03, 79.11, 271.81, 146.20, 74.88)
burn <- c(rep("unburned", 3), rep("burned", 3))
diversity_order <- c("q0", "q1", "q2", "q0", "q1", "q2")
bact_q_table <- cbind(burn, diversity_order, bact_mean, bact_sd)
bact_q_table <- data.frame(bact_q_table)

bact_q_table$bact_mean <- as.numeric(bact_q_table$bact_mean)
bact_q_table$bact_sd <- as.numeric(bact_q_table$bact_sd)

bact_q_table$se <- bact_q_table$bact_sd/sqrt(3)
bact_q_table$lower_ci <- bact_q_table$bact_mean - bact_q_table$se
bact_q_table$upper_ci <- bact_q_table$bact_mean + bact_q_table$se

pd <- position_dodge(0.3)
bact_q_table$burn <- factor(bact_q_table$burn, levels = c("unburned", "burned"))

(bact_div_fig <- ggplot(bact_q_table, aes(x = diversity_order, y = bact_mean, group = burn)) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), color = "black", width = 0.15, position = pd) +
    geom_point(aes(fill = burn), position = pd, size = 6, shape = 21, stroke = 1.5) +
    scale_fill_manual(values = c("#4363d8", "#f58231")) +
    scale_y_continuous(limits = c(150, 1070), breaks = c(250, 500, 750, 1000)) +
    xlab("Order of diversity (q)") +
    ylab("Hill number") +
    theme_classic(base_size = 10, base_family = "") +
    theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2),
          axis.title.y = element_text(face = "bold", color = "black", size = 26),
          axis.title.x = element_blank(),
          axis.text.y = element_text(face = "plain", color = "black", size = 24),
          axis.text.x = element_blank(),
          axis.ticks.length = unit(0.15, "cm"),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          legend.position = "none"))


## Fungi ##
fungi_q0_values <- c(565, 413, 441, 141, 158, 239, 159)
fungi_q1_values <- c(192, 103, 92.3, 22.9, 21.5, 53.4, 16.5)
fungi_q2_values <- c(59.2, 37.9, 28.2, 9.8, 11.1, 22.4, 6.5)
fungi_q_table <- cbind(sample_names, burn_severity, fungi_q0_values, fungi_q1_values, fungi_q2_values)
colnames(fungi_q_table) <- c("sample", "burn_severity", "q0", "q1", "q2")


# fungi q0 unburned mean 
mean(c(565, 413, 441))
# 473

# fungi q0 unburned sd 
sd(c(565, 413, 441))
# 80.89

# fungi q1 unburned mean 
mean(c(192, 103, 92.3))
# 129.1

# fungi q1 unburned sd 
sd(c(192, 103, 92.3))
# 54.74

# fungi q2 unburned mean 
mean(c(59.2, 37.9, 28.2))
# 41.77

# fungi q2 unburned sd 
sd(c(59.2, 37.9, 28.2))
# 15.86

# fungi q0 burned mean 
mean(c(141, 158, 239, 159))
# 174.25

# fungi q0 burned sd 
sd(c(141, 158, 239, 159))
# 43.95

# fungi q1 burned mean 
mean(c(22.9, 21.5, 53.4, 16.5))
# 28.58

# fungi q1 burned sd 
sd(c(22.9, 21.5, 53.4, 16.5))
# 16.78

# fungi q2 burned mean 
mean(c(9.8, 11.1, 22.4, 6.5))
# 12.45 

# fungi q2 burned sd 
sd(c(9.8, 11.1, 22.4, 6.5))
# 6.91

fungi_mean <- c(473, 129.1, 41.77, 174.25, 28.58, 12.45)
fungi_sd <- c(80.89, 54.74, 15.86, 43.95, 16.78, 6.91)
fungi_q_table <- cbind(burn, diversity_order, fungi_mean, fungi_sd)
fungi_q_table <- data.frame(fungi_q_table)

fungi_q_table$fungi_mean <- as.numeric(fungi_q_table$fungi_mean)
fungi_q_table$fungi_sd <- as.numeric(fungi_q_table$fungi_sd)

fungi_q_table$se <- fungi_q_table$fungi_sd/sqrt(3)
fungi_q_table$lower_ci <- fungi_q_table$fungi_mean - fungi_q_table$se
fungi_q_table$upper_ci <- fungi_q_table$fungi_mean + fungi_q_table$se

pd <- position_dodge(0.3)
fungi_q_table$burn <- factor(fungi_q_table$burn, levels = c("unburned", "burned"))

(fungi_div_fig <- ggplot(fungi_q_table, aes(x = diversity_order, y = fungi_mean, group = burn)) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), color = "black", width = 0.15, position = pd) +
    geom_point(aes(fill = burn), position = pd, size = 6, shape = 21, stroke = 1.5) +
    scale_fill_manual(values = c("#4363d8", "#f58231")) +
    scale_y_continuous(limits = c(0, 590), breaks = c(100, 200, 300, 400, 500)) +
    xlab("Order of diversity (q)") +
    ylab("Hill number") +
    theme_classic(base_size = 10, base_family = "") +
    theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2),
          axis.title = element_text(face = "bold", color = "black", size = 26),
          axis.text = element_text(face = "plain", color = "black", size = 24),
          axis.ticks.length = unit(0.15, "cm"),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          legend.position = "none"))


## Plant pathogenic fungi ##
path_q0_values <- c(13, 14, 11, 5, 5, 8, 3)
path_q1_values <- c(7.35, 8.71, 7.23, 2.54, 3.59, 3.66, 2.13)
path_q2_values <- c(5.54, 6.58, 5.23, 2.04, 3.32, 2.32, 1.79)

# Path q0 unburned mean 
mean(c(13, 14, 11))
# 12.67

# Path q0 unburned sd 
sd(c(13, 14, 11))
# 1.53

# Path q0 burned mean 
mean(c(5, 5, 8, 3))
# 5.25

# Path q0 burned sd 
sd(c(5, 5, 8, 3))
# 2.06

# Path q1 unburned mean 
mean(c(7.35, 8.71, 7.23))
# 7.76

# Path q1 unburned sd 
sd(c(7.35, 8.71, 7.23))
# 0.82

# Path q1 burned mean 
mean(c(2.54, 3.59, 3.66, 2.13))
# 2.98

# Path q1 burned sd 
sd(c(2.54, 3.59, 3.66, 2.13))
# 0.76

# Path q2 unburned mean
mean(c(5.54, 6.58, 5.23))
# 5.78

# Path q2 unburned sd
sd(c(5.54, 6.58, 5.23))
# 0.71

# Path q2 burned mean
mean(c(2.04, 3.32, 2.32, 1.79))
# 2.37

# Path q2 burned sd
sd(c(2.04, 3.32, 2.32, 1.79))
# 0.67

path_mean <- c(12.67, 7.76, 5.78, 5.25, 2.98, 2.37)
path_sd <- c(1.53, 0.82, 0.71, 2.06, 0.76, 0.67)
path_q_table <- cbind(burn, diversity_order, path_mean, path_sd)
path_q_table <- data.frame(path_q_table)

path_q_table$path_mean <- as.numeric(path_q_table$path_mean)
path_q_table$path_sd <- as.numeric(path_q_table$path_sd)

path_q_table$se <- path_q_table$path_sd/sqrt(3)
path_q_table$lower_ci <- path_q_table$path_mean - path_q_table$se
path_q_table$upper_ci <- path_q_table$path_mean + path_q_table$se

pd <- position_dodge(0.3)
path_q_table$burn <- factor(path_q_table$burn, levels = c("unburned", "burned"))

(path_div_fig <- ggplot(path_q_table, aes(x = diversity_order, y = path_mean, group = burn)) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), color = "black", width = 0.15, position = pd) +
    geom_point(aes(fill = burn), position = pd, size = 6, shape = 21, stroke = 1.5) +
    scale_fill_manual(values = c("#4363d8", "#f58231")) +
    scale_y_continuous(limits = c(1.75, 14), breaks = c(2, 4, 6, 8, 10, 12)) +
    xlab("Order of diversity (q)") +
    ylab("Hill number") +
    theme_classic(base_size = 10, base_family = "") +
    theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2),
          axis.title = element_text(face = "bold", color = "black", size = 26),
          axis.text = element_text(face = "plain", color = "black", size = 24),
          axis.ticks.length = unit(0.15, "cm"),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          legend.position = "none"))


# Plant pathogenic fungi relative abundance
pd <- position_dodge(0.3)
plant_path_abund_stats$burn_category <- factor(plant_path_abund_stats$burn_category, levels = c("Unburned", "Burned"))

(plant_path_abundance_fig <- ggplot(plant_path_abund_stats, aes(x = burn_category, y = Mean)) +
    geom_errorbar(aes(ymin = Lower_ci, ymax = Upper_ci), color = "black", width = 0.15, position = pd) +
    geom_point(aes(fill = burn_category), position = pd, size = 6, shape = 21, stroke = 1.5) +
    scale_fill_manual(values = c("#4363d8", "#f58231")) +
    scale_y_continuous(limits = c(1.25, 6.75), breaks = c(2, 3, 4, 5, 6)) +
    xlab("Soil burn severity") +
    ylab("Proportion (%)") +
    theme_classic(base_size = 10, base_family = "") +
    theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2),
          axis.title = element_text(face = "bold", color = "black", size = 26),
          axis.text = element_text(face = "plain", color = "black", size = 24),
          axis.ticks.length.y = unit(0.15, "cm"),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          legend.position = "none"))

(combined_patch <- (bact_div_fig/fungi_div_fig) | (path_div_fig/plant_path_abundance_fig))

