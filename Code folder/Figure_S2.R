################################################################################
################################## Figure S2 ###################################
################################################################################

# Loading the R packages
library(openxlsx)
library(car)
library(ggplot2)
library(effectsize)
library(vegan)
library(glmm.hp)
library(dplyr)

# Loading field survey data
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)
Field_group$Years <- as.factor(Field_group$Years)
Field_group$Site <- factor(Field_group$Site, levels = unique(Field_group$Site[order(Field_group$Latitude)]))

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
Field_group$Wcont <- sqrt(Field_group$Wcont)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Phy_Di_log <- log10(Field_group$Phy_Di)
Field_group$Fun_Di_log <- log10(Field_group$Fun_Di)

# Consider normalizing your data
pd_attributes_variable <- attributes(scale(Field_group[c("Site_pool", "Phy_Di", "Fun_Di", "Phy_Di_log", "Fun_Di_log", "Soil_ph", "Wcont", "Soil_N", "Tave", "Precipitation")]))
total_data <- Field_group
total_data[c("Site_pool", "Phy_Di", "Fun_Di", "Phy_Di_log", "Fun_Di_log", "Soil_ph", "Wcont", "Soil_N", "Tave", "Precipitation")] <-  
  scale(total_data[c("Site_pool", "Phy_Di", "Fun_Di", "Phy_Di_log", "Fun_Di_log", "Soil_ph", "Wcont", "Soil_N", "Tave", "Precipitation")])

total_data <- total_data %>% left_join(Field_richness)
##
fm1 <- lm(Field_SR ~ Family/Species + Site + Years + 
            Family:Site + Family:Years + Family:Years:Site + 
            Family/Species:Site + Family/Species:Years + Years:Site, data = total_data)

#car::Anova(fm1, type = 2)
anova(fm1)


##
fm1 <- lm(Field_SR ~ Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Precipitation + Tave +
             Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + Phy_Di_log:Soil_N + 
             Fun_Di_log:Tave + Fun_Di_log:Precipitation + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont + Fun_Di_log:Soil_N, data = total_data)


options(na.action = "na.fail")
dd12 <- dredge(fm1, subset =
                 dc(Fun_Di_log, Tave, Fun_Di_log:Tave) & 
                 dc(Fun_Di_log, Precipitation, Fun_Di_log:Precipitation) & 
                 dc(Fun_Di_log, Soil_ph, Fun_Di_log:Soil_ph) & 
                 dc(Fun_Di_log, Wcont, Fun_Di_log:Wcont) & 
                 dc(Fun_Di_log, Soil_N, Fun_Di_log:Soil_N) &
                 ##
                 dc(Phy_Di_log, Tave, Phy_Di_log:Tave) & 
                 dc(Phy_Di_log, Precipitation, Phy_Di_log:Precipitation) & 
                 dc(Phy_Di_log, Soil_ph, Phy_Di_log:Soil_ph) & 
                 dc(Phy_Di_log, Wcont, Phy_Di_log:Wcont) & 
                 dc(Phy_Di_log, Soil_N, Phy_Di_log:Soil_N), trace = 2)

# Identified the best-fitting models for each response variable based on Corrected Akaike Information Criterion (AICc)
Final_model <- get.models(dd12,1)[[1]] 
summary(Final_model)

cor.test(total_data$Soil_N, total_data$Field_SR)
cor.test(total_data$Precipitation, total_data$Field_SR)

################################ Table_Fig_S2  #################################
Table_Fig_S2 <- as.data.frame(car::Anova(Final_model, type = 2))
Table_Fig_S2 <- Table_Fig_S2[-which(rownames(Table_Fig_S2) == "Residuals"),]
Table_Fig_S2$p.adj <- p.adjust(Table_Fig_S2$`Pr(>F)`, method = "holm")
Table_Fig_S2$p.adj <- sprintf("%.3f", Table_Fig_S2$p.adj)
Table_Fig_S2$`Pr(>F)` <- round(Table_Fig_S2$`Pr(>F)`, 3)
Table_Fig_S2$`F` <- round(Table_Fig_S2$`F value`, 2)
Table_Fig_S2$Parameter <- rownames(Table_Fig_S2)
Table_Fig_S2$VIF <- round(car::vif(Final_model), 2)
Table_Fig_S2$df <- paste0("1", ",", "377")
rownames(Table_Fig_S2) <- c("Funct-Dist", "Phylo-Dist", "Precipitation", "Soil N", "Temperature",
                        "Funct-Dist × Precipitation", "Funct-Dist × Temperature")
print(Table_Fig_S2[c(1,2,3,5,4,6,7) ,c("df","F","Pr(>F)","p.adj","VIF")])
performance::r2(Final_model)

# Obtaining standardized regression coefficients and their 95% CI
MegaModelSummary <- as.data.frame(effectsize(Final_model))[-1,]
MegaModelSummary$Parameter2 <- c("Funct-Dist","Phylo-Dist","Precipitation", "Soil N", "Temperature","Funct-Dist × Precipitation","Funct-Dist × Temperature")
MegaModelSummary$Parameter2 <- factor(MegaModelSummary$Parameter2, levels = rev(c("Funct-Dist","Phylo-Dist","Precipitation", "Temperature", "Soil N",
                                                                                  "Funct-Dist × Precipitation","Funct-Dist × Temperature")))
MegaModelSummary <- MegaModelSummary %>% left_join(Table_Fig_S2[,c("Parameter","p.adj")], by = "Parameter")

MegaModelSummary$Group <- c("Plant attributes", "Plant attributes", "Climate", "Soil properties", "Climate", "Interaction", "Interaction")
MegaModelSummary$Group <- factor(MegaModelSummary$Group, levels = c("Plant attributes", "Climate", "Soil properties", "Interaction"))

ggplot(MegaModelSummary, aes(x = Parameter2, y = Std_Coefficient, fill = Group))+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, size = 0.8, color = "black")+
  geom_point(size = 3.5, pch = 21)+
  geom_segment(aes(y = 0, yend = 0, x = 0.5, xend = 7.35), color = "black", linetype = "dashed") + 
  geom_text(aes(y = CI_high, label = paste("italic(p)==", `p.adj`)),
            parse = TRUE, hjust = -0.4, vjust = 0.4, size = 3.6) + 
  annotate("text", x = 7.6, y = 0.08,
           label = expression("Best model: " * italic(R)^2 * "= 0.08"), colour = "black", size = 4) + 
  annotate("text", x = 5.0, y = -0.10,
           label = expression(italic(p) * " = 0.022"), colour = "black", size = 3.6) + 
  labs(x = '', y = 'Standard regression coefficients') +  
  theme_classic() + coord_flip() +  
  scale_fill_manual(values = c("#436C88","#49A3A4","#C2887D","#F0C986")) +
  theme(axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.x =  element_text(color = "black", size = 13),
        legend.text = element_text(size = 9, color = "black"),
        plot.margin = margin(0.5,1.5,0.5,1.5, unit = "cm"),
        legend.position = 'none',legend.title = element_blank(),legend.key = element_blank(),
        plot.tag = element_text(size = 14, face = "bold")) +
  scale_x_discrete(expand = expansion(mult = c(0.00, 0.12))) + 
  scale_shape_manual(values = c(16,21)) -> p_a; p_a

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
  scale_fill_manual(values = c("#436C88","#49A3A4","#C2887D","#F0C986")) +
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
  labs(y = "Relative effect of estimates (%)") -> p_b; p_b

# combination
(p_a+p_b) + plot_layout(widths = c(0.7,0.3)) -> Fig_S2; Fig_S2

