################################################################################
############################# Figure S7 & Table S3 #############################
################################################################################

# Loading the R packages
library(openxlsx)
library(car)
library(ggplot2)
library(effectsize)
library(vegan)
library(glmm.hp)
library(dplyr)
library(patchwork)
library(ggtext)

# Loading field survey data
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)
Field_group$Years <- as.factor(Field_group$Years)
Field_group$Site <- factor(Field_group$Site, levels = unique(Field_group$Site[order(Field_group$Latitude)]))
Field_group$Species <- as.factor(Field_group$Species)

# Soil sample abundance information
Field_otu_raw <- read.xlsx("Field_data_raw_ASVs.xlsx", sheet = "raw_otu", colNames = T, rowNames = T)
Field_otu_raw <- Field_otu_raw[ ,Field_group$Sample_ID]
dim(Field_otu_raw)

# Richness of overall fungi
Field_richness <- as.data.frame(specnumber(t(Field_otu_raw)))
Field_richness$Sample_ID = rownames(Field_richness); colnames(Field_richness)[1] = "Field_SR"
Rmisc::summarySE(Field_richness, measurevar = "Field_SR")

# Data Transformation
Field_group$RS <- sqrt(Field_group$RS)
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont*100)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Funct_Di_log <- log10(Field_group$Funct_Di)
Field_group$Phylo_Di_log <- log10(Field_group$Phylo_Di)

colnames(Field_group)
# Consider normalizing your data
var_select <- c("Site_pool","Phylo_Di","Funct_Di","Phylo_Di_log","Funct_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation","PCoA1")
pd_attributes_variable <- attributes(scale(Field_group[var_select]))
total_data <- Field_group
total_data[var_select] <- scale(total_data[var_select])

# add dataset of fungal richness 
total_data <- total_data %>% left_join(Field_richness)
total_data$Field_SR <- as.numeric(total_data$Field_SR)

######################### Table S3 (Field survey part) #########################
# build linear model
fm1_SR <- lm(Field_SR ~ Family/Species + Site + Years +
               Family:Site + Family:Years + Family:Years:Site + 
               Family/Species:Site + Family/Species:Years + Years:Site, data = total_data)

anova(fm1_SR)

# Calculate the proportion of variance explained by each factor
anova_result <- anova(fm1_SR)
anova_table <- anova_result
total_SS <- sum(anova_table$`Sum Sq`)

# Create a results table containing the proportion of variance explained
variance_table <- data.frame(
  Effect = rownames(anova_table),
  Df = anova_table$Df,
  Sum_Sq = anova_table$`Sum Sq`,
  Mean_Sq = anova_table$`Mean Sq`,
  F_value = anova_table$`F value`,
  P_value = anova_table$`Pr(>F)`,
  Variance_Explained = (anova_table$`Sum Sq` / total_SS))

print(variance_table)


################################## Figure S7 ###################################
# performed the global linear regression model
fm1 <- lm(Field_SR ~ Phylo_Di_log + Funct_Di_log + Soil_ph + Wcont + Soil_N + Precipitation + Tave +
            Phylo_Di_log:Tave + Phylo_Di_log:Precipitation + Phylo_Di_log:Soil_ph + Phylo_Di_log:Wcont + Phylo_Di_log:Soil_N + 
            Funct_Di_log:Tave + Funct_Di_log:Precipitation + Funct_Di_log:Soil_ph + Funct_Di_log:Wcont + Funct_Di_log:Soil_N, data = total_data)


options(na.action = "na.fail")
dd12 <- dredge(fm1, subset = 
                 dc(Funct_Di_log, Tave, Funct_Di_log:Tave) & 
                 dc(Funct_Di_log, Precipitation, Funct_Di_log:Precipitation) & 
                 dc(Funct_Di_log, Soil_ph, Funct_Di_log:Soil_ph) & 
                 dc(Funct_Di_log, Wcont, Funct_Di_log:Wcont) & 
                 dc(Funct_Di_log, Soil_N, Funct_Di_log:Soil_N) &
                 ##
                 dc(Phylo_Di_log, Tave, Phylo_Di_log:Tave) & 
                 dc(Phylo_Di_log, Precipitation, Phylo_Di_log:Precipitation) & 
                 dc(Phylo_Di_log, Soil_ph, Phylo_Di_log:Soil_ph) & 
                 dc(Phylo_Di_log, Wcont, Phylo_Di_log:Wcont) & 
                 dc(Phylo_Di_log, Soil_N, Phylo_Di_log:Soil_N), trace = 2, rank = "AICc")

# Identified the best-fitting models for each response variable based on Corrected Akaike Information Criterion (AICc)
Final_model <- get.models(dd12,1)[[1]] 
performance::r2(Final_model)
AICc(Final_model)

############################ Table of Figure S7  ###############################
Table_Fig_S7 <- as.data.frame(car::Anova(Final_model, type = 2))
Table_Fig_S7 <- Table_Fig_S7[-which(rownames(Table_Fig_S7) == "Residuals"),]
Table_Fig_S7$p.adj <- p.adjust(Table_Fig_S7$`Pr(>F)`, method = "holm")
Table_Fig_S7$p.adj <- sprintf("%.3f", Table_Fig_S7$p.adj)
Table_Fig_S7$`Pr(>F)` <- round(Table_Fig_S7$`Pr(>F)`, 3)
Table_Fig_S7$`F value` <- round(Table_Fig_S7$`F value`, 2)
Table_Fig_S7$Parameter <- rownames(Table_Fig_S7)
Table_Fig_S7$VIF <- round(car::vif(Final_model), 2)
Table_Fig_S7$df <- paste0("1", ",", "377")
rownames(Table_Fig_S7) <- c("Funct-Dist", "Phylo-Dist", "Precipitation", "Soil N", "Temperature",
                            "Funct-Dist × Precipitation", "Funct-Dist × Temperature")
print(Table_Fig_S7[c(1,2,3,5,4,6,7) ,c("df","F value","Pr(>F)","p.adj","VIF")])


# Obtaining standardized regression coefficients and their 95% CI
MegaModelSummary <- as.data.frame(effectsize(Final_model))[-1,]
MegaModelSummary$Parameter2 <- c("Funct-Dist","Phylo-Dist","Precipitation", "Soil N", "Temperature","Funct-Dist × Precipitation","Funct-Dist × Temperature")
MegaModelSummary$Parameter2 <- factor(MegaModelSummary$Parameter2, levels = rev(c("Funct-Dist","Phylo-Dist","Precipitation", "Temperature", "Soil N",
                                                                                  "Funct-Dist × Precipitation","Funct-Dist × Temperature")))
MegaModelSummary <- MegaModelSummary %>% left_join(Table_Fig_S7[,c("Parameter","p.adj")], by = "Parameter")

MegaModelSummary$Group <- c("Plant attributes", "Plant attributes", "Climate", "Soil properties", "Climate", "Interaction", "Interaction")
MegaModelSummary$Group <- factor(MegaModelSummary$Group, levels = c("Plant attributes", "Climate", "Soil properties", "Interaction"))

ggplot(MegaModelSummary, aes(x = Parameter2, y = Std_Coefficient, fill = Group))+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, size = 0.8, color = "black")+
  geom_point(size = 3.5, pch = 21)+
  geom_segment(aes(y = 0, yend = 0, x = 0.5, xend = 7.35), color = "black", linetype = "dashed") + 
  geom_text(aes(y = CI_high, label = paste("italic(p)==", `p.adj`)),
            parse = TRUE, hjust = -0.4, vjust = 0.4, size = 3.6) + 
  annotate("text", x = 5.0, y = -0.10,
           label = expression(italic(p) * " = 0.022"), colour = "black", size = 3.6) +
  labs(x = '', 
       y = 'Standard regression coefficients', 
       title = "Best model: <i>R</i><sup>2</sup> = 0.082") +  
  theme_classic() + coord_flip() +  
  scale_fill_manual(values = c("#49A3A4","#80B09A","#C2887D","#F0C986")) +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x =  element_text(color = "black", size = 14),
        legend.text = element_text(size = 9, color = "black"),
        plot.title = element_textbox(size = 12, color = "black", fill = "white",     
                                     box.color = "black", width = grid::unit(1, "npc"),padding = margin(5, 5, 5, 5),  
                                     margin = margin(b = 5), halign = 0.5,linetype = "solid"),
        plot.margin = margin(0.5,1.5,0.5,1.5, unit = "cm"),
        legend.position = 'none',legend.title = element_blank(),legend.key = element_blank(),
        plot.tag = element_text(size = 14, face = "bold")) +
  scale_x_discrete(expand = expansion(mult = c(0.00, 0.12))) + 
  scale_shape_manual(values = c(16,21)) -> Figure_S7a; Figure_S7a

# Relative contribution of variables
hierarchical_data <- as.data.frame(glmm.hp(Final_model, type = "R2")$hierarchical.partitioning)
hierarchical_data$Group <- c("Plant attributes", "Plant attributes", "Climate", "Soil properties", "Climate", "Interaction", "Interaction")
hierarchical_data$Parameter2 <- c("Funct-Dist","Phylo-Dist","Precipitation", "Temperature", "Soil N",
                                  "Funct-Dist × Precipitation","Funct-Dist × Temperature")
hierarchical_data$Group <- factor(hierarchical_data$Group, levels = c("Plant attributes", "Climate", "Soil properties", "Interaction"))
hierarchical_data$`I.perc(%)` <- (hierarchical_data$`I.perc(%)`/sum(hierarchical_data$`I.perc(%)`)*100)

hierarchical_data2 <- hierarchical_data %>% dplyr::group_by(Group) %>% 
  dplyr::summarise(`I.perc(%)` = sum(`I.perc(%)`)) %>% 
  mutate(label = paste0(round(`I.perc(%)`, 1), "%")) 

ggplot()+
  geom_bar(data = hierarchical_data2, aes(x = "", y = `I.perc(%)`, fill = Group), 
           stat = "identity", width = 0.5, color = "black")+
  geom_text(data = hierarchical_data2, 
            aes(x = "", y = `I.perc(%)`, label = label, group = Group), 
            position = position_stack(vjust = 0.5), color = "black", size = 4) +
  scale_y_continuous(position = "right", limits = c(0, 100), expand = c(0, 0)) +
  scale_fill_manual(values = c("#49A3A4","#80B09A","#C2887D","#F0C986")) +
  theme_classic()+
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.y = element_text(color = "black", size = 13),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_line(linewidth = NA))+
  labs(y = "Relative effect of estimates (%)") -> Figure_S7b; Figure_S7b

# combination
(Figure_S7a + Figure_S7b) + plot_layout(widths = c(0.7,0.3)) -> Figure_S7; Figure_S7


####################### Table S3 (Greenhouse exp. part) ########################
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
Green_group <- Green_group %>% left_join(Green_richness)
Green_SR <- lm(Green_SR ~ Family/Species , data = Green_group)
anova(Green_SR)

Table_Fig_S3_green <- as.data.frame(anova(Green_SR))
#Table_Fig_S3_green <- Table_Fig_S3_green[-which(rownames(Table_Fig_S3_green) == "Residuals"),]
Table_Fig_S3_green$`Pr(>F)` <- round(Table_Fig_S3_green$`Pr(>F)`, 3)
Table_Fig_S3_green$`F value` <- round(Table_Fig_S3_green$`F value`, 2)
rownames(Table_Fig_S3_green)[1:2] <- c("Family", "Species")
Table_Fig_S3_green$Variance_explained <- round(Table_Fig_S3_green$`Sum Sq`/sum(Table_Fig_S3_green$`Sum Sq`), 3)
print(Table_Fig_S3_green)

