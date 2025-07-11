################################################################################
################################## Figure S3 ###################################
################################################################################

# Loading the R packages
library(ggplot2)
library(openxlsx)
library(dplyr)
library(patchwork)
library(vegan)
library(car)
library(phytools)

# Custom Style
mytheme <- theme(panel.background = element_rect(fill = 'transparent', colour='black'),
                 legend.position = 'none',
                 panel.grid=element_blank(), 
                 legend.title = element_blank(),
                 legend.text = element_text(size = 8),
                 legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
                 axis.ticks = element_line(color='black'),
                 axis.line = element_line(colour = "black"), 
                 axis.title.x = element_text(colour='black', size=13),
                 axis.title.y = element_text(colour='black', size=13),
                 axis.text = element_text(colour='black',size=11),
                 plot.tag = element_text(size = 14, face = "bold")) 

family_color <- c("Poaceae"="#332500","Cyperaceae"="#542D20","Solanaceae" ="#994240",
                  "Verbenaceae" = "#D75C4D", "Acanthaceae" = "#E68C51", "Lamiaceae" ="#F59D52", 
                  "Asteraceae" = "#EFBA55", "Polygonaceae" ="#FCD170", "Phytolaccaceae"="#FEE1B0", 
                  "Caryophyllaceae" = "#C5E0D0", "Amaranthaceae" = "#ABDCE0", "Euphorbiaceae" ="#7FC5D8",
                  "Onagraceae"="#73BCD5", "Urticaceae" = "#528FAC", "Fabaceae" = "#376694", "Malvaceae" = "#1F466F")

family_order <- c("Poaceae", "Cyperaceae", "Solanaceae", "Verbenaceae", "Acanthaceae", "Lamiaceae", 
                  "Asteraceae", "Polygonaceae", "Phytolaccaceae", "Caryophyllaceae", "Amaranthaceae", 
                  "Euphorbiaceae", "Onagraceae", "Urticaceae", "Fabaceae", "Malvaceae")

# Loading the greenhouse experiment data
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "green_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)

Green_otu_raw <- read.xlsx("Greenhouse_data_raw_ASVs.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
Green_otu_raw <- Green_otu_raw[ ,Green_group$Sample_ID]
dim(Green_otu_raw)

# Richness of overall fungi
Green_richness <- as.data.frame(specnumber(t(Green_otu_raw)))
Green_richness$Sample_ID <- rownames(Green_richness); colnames(Green_richness)[1] <- "Green_SR"

# Effect of family and species on fungal richness
Green_exp_data <- Green_group %>% left_join(Green_richness)
mod <- lm(Green_SR ~ Family/Species , data = Green_exp_data)
Table_Fig_S3 <- as.data.frame(car::Anova(mod, type = 2))
Table_Fig_S3 <- Table_Fig_S3[-which(rownames(Table_Fig_S3) == "Residuals"),]
Table_Fig_S3$p.adj <- p.adjust(Table_Fig_S3$`Pr(>F)`, method = "holm")
Table_Fig_S3$p.adj <- sprintf("%.3f", Table_Fig_S3$p.adj)
Table_Fig_S3$`Pr(>F)` <- round(Table_Fig_S3$`Pr(>F)`, 3)
Table_Fig_S3$`F` <- round(Table_Fig_S3$`F value`, 2)
rownames(Table_Fig_S3) <- c("Family", "Species")
Table_Fig_S3$df <- c("15,108", "38,108")
print(Table_Fig_S3[ ,c("df","F","Pr(>F)","p.adj")])

# loading phylogenetic tree
tree <- read.newick("IQ_tree_plant_2025.NEWICK")
to_drop <- c("Amborella_trichopoda","")
tree <- drop.tip(as.phylo(tree), to_drop) 

################################# Figure S2A ###################################
Green_exp_data$Family <- factor(Green_exp_data$Family, levels = rev(family_order))

ggplot(Green_exp_data,aes(x = Family, y = Green_SR, fill = Family)) + 
  geom_boxplot(width = 0.6,alpha = 0.75) + 
  #scale_color_manual(values = family_color) +
  scale_fill_manual(values = family_color) +
  annotate("text", x = 16, y = 120,
           label = expression("F"[15*","*108] == paste("1.36")), size = 4) + 
  annotate("text", x = 15.5, y = 120,
           label = expression(italic(p) == 0.216), size = 4) +
  mytheme + 
  coord_flip() + 
  labs(x = NULL, y = "Overall fungal richness", tag = "(a)") -> p1; p1

################################# Figure S2B ###################################
Green_exp_data$Species <- gsub("_", " ", Green_exp_data$Species)
Green_exp_data$Species <- factor(Green_exp_data$Species, levels = gsub("_", " ", tree$tip.label))
#str(Green_exp_data)
#Green_exp_data$sig <- ifelse(Green_exp_data$Origin == "Native"," ","*")

ggplot(Green_exp_data,aes(x = Species, y = Green_SR, fill = Family)) + 
  geom_boxplot(width = 0.6, alpha = 0.75) + 
  #scale_color_manual(values = family_color) +
  scale_fill_manual(values = family_color) +
  annotate("text", x = 53, y = 120,
           label = expression("F"[38*","*108] == paste("1.37")), size = 4) + 
  annotate("text", x = 51, y = 120,
           label = expression(italic(p) == 0.216), size = 4) +
  mytheme + 
  theme(axis.text.y = element_text(face = "italic",size=10)) + 
  coord_flip() + 
  labs(x = NULL, y = "Overall fungal richness", tag = "(b)") -> p2; p2

# Combination 
p1|p2 -> Fig_S3; Fig_S3

