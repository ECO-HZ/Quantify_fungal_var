################################################################################
################################## Figure S1 ###################################
################################################################################

# Loading the R packages
library(openxlsx)
library(nlme)
library(car)
library(MuMIn)
library(ggplot2)
library(effectsize)
library(ggtext)
library(glmm.hp)
library(dplyr)

# Loading field survey data
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)
Field_group$Years <- as.factor(Field_group$Years)
Field_group$Site <- factor(Field_group$Site, levels = unique(Field_group$Site[order(Field_group$Latitude)]))

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

fm1 <- lme(Fungal_SR ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Precipitation + Tave +
            Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + Phy_Di_log:Soil_N + 
            Fun_Di_log:Tave + Fun_Di_log:Precipitation + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont + Fun_Di_log:Soil_N, 
          random = ~1|Years, method = "ML", data = total_data)

options(na.action = "na.fail")
dd12 <- dredge(fm1, subset = ~ Site_pool &
                 dc(Fun_Di_log, Tave, Fun_Di_log:Tave) & 
                 dc(Fun_Di_log, Precipitation, Fun_Di_log:Precipitation) & 
                 dc(Fun_Di_log, Wcont, Fun_Di_log:Wcont) & 
                 dc(Fun_Di_log, Tave, Fun_Di_log:Tave) & 
                 dc(Fun_Di_log, Soil_N, Fun_Di_log:Soil_N) &
                 ##
                 dc(Phy_Di_log, Tave, Phy_Di_log:Tave) & 
                 dc(Phy_Di_log, Precipitation, Phy_Di_log:Precipitation) & 
                 dc(Phy_Di_log, Wcont, Phy_Di_log:Wcont) & 
                 dc(Phy_Di_log, Tave, Phy_Di_log:Tave) & 
                 dc(Phy_Di_log, Soil_N, Phy_Di_log:Soil_N), trace = 2)

# Identified the best-fitting models for each response variable based on Corrected Akaike Information Criterion (AICc)
Final_model <- get.models(dd12,1)[[1]] 
Final_model <- update(Final_model, ~., method = "REML") # Update to REML to extract estimates
summary(Final_model)
r.squaredGLMM(Final_model)
Anova_table <- as.data.frame(car::Anova(Final_model, type = 3))
Anova_table <- Anova_table[-which(rownames(Anova_table) == "(Intercept)"),]
Anova_table <- as.data.frame(Anova_table)
Anova_table$`q-vaules` <- p.adjust(Anova_table$`Pr(>Chisq)`, method = "holm")
Anova_table$`q-vaules` <- sprintf("%.3f", Anova_table$`q-vaules`)
Anova_table$`p-value` <- round(Anova_table$`Pr(>Chisq)`, 3)
Anova_table$Parameter <- rownames(Anova_table)

# Obtaining standardized regression coefficients and their 95% CI
MegaModelSummary <- as.data.frame(effectsize(Final_model))[-1,]
MegaModelSummary$Parameter2 <- c("Funct-Dist","Phylo-Dist","Precipitation", "Site pool","Soil N","Tave","Funct-Dist × Soil N","Phylo-Dist × Tave")
MegaModelSummary$Parameter2 <- factor(MegaModelSummary$Parameter2, levels = rev(c("Site pool","Soil N", "Precipitation","Tave",
                                                                                  "Funct-Dist","Phylo-Dist","Phylo-Dist × Tave", "Funct-Dist × Soil N")))
MegaModelSummary <- MegaModelSummary %>% left_join(Anova_table[,c("Parameter","q-vaules")], by = "Parameter")

MegaModelSummary$Group <- c("Plant attributes", "Plant attributes", "Climate", "Site pool", "Soil properties", "Climate", "Interaction", "Interaction")
MegaModelSummary$Group <- factor(MegaModelSummary$Group, levels = c("Site pool", "Plant attributes", "Soil properties", "Climate", "Interaction"))

ggplot(MegaModelSummary, aes(x = Parameter2, y = Std_Coefficient, fill = Group))+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, size = 0.8, color = "black")+
  geom_point(size = 3.5, pch = 21)+
  geom_segment(aes(y = 0, yend = 0, x = 0.5, xend = 8.35), color = "black", linetype = "dashed") + 
  geom_text(aes(y = CI_high, label = paste("italic(p)==", `q-vaules`)),
            parse = TRUE, hjust = -0.4, vjust = 0.4, size = 3.5) + 
  annotate("text", x = 8.6, y = 0.08,
           label = expression("Best model: " * italic(R)^2 * "m = 0.09, " * italic(R)^2 * "c = 0.13"), colour = "black", size = 4) + 
  annotate("text", x = 6.0, y = 0.03,
           label = expression(italic(p) * " < 0.001"), colour = "black", size = 3.5) + 
  labs(x = '', y = 'Standard regression coefficients') +  
  theme_classic() + coord_flip() +  
  scale_fill_manual(values = c("Site pool" = "#949698","Plant attributes" = "#DE7963", "Soil properties" = "#3F425A",
                               "Climate" = "#82B19D","Interaction" = "#F0C986")) +
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
hierarchical_data$Group <- c("Plant attributes", "Plant attributes", "Climate", "Site pool", "Soil properties", "Climate", "Interaction", "Interaction")
hierarchical_data$Parameter2 <- c("Funct-Dist","Phylo-Dist","Precipitation", "Site pool","Soil N","Tave","Funct-Dist × Soil N","Phylo-Dist × Tave")
hierarchical_data$Group <- factor(hierarchical_data$Group, levels = c("Site pool", "Plant attributes", "Soil properties", "Climate", "Interaction"))
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
  scale_y_continuous(expand = c(0, 0), position = "right")+
  scale_fill_manual(values = c("Site pool" = "#949698","Plant attributes" = "#DE7963", "Soil properties" = "#3F425A",
                               "Climate" = "#82B19D","Interaction" = "#F0C986")) +
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
(p_a+p_b) + plot_layout(widths = c(0.7,0.3)) -> Fig_S1; Fig_S1
