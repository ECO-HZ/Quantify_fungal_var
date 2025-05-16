################################################################################
############################ Figure 3 & Table S5 ###############################
################################################################################

# Loading the R packages
library(openxlsx)
library(nlme)
library(car)
library(MuMIn)
library(ggplot2)
library(effectsize)
library(glmm.hp)
library(effects)
library(patchwork)

# Custom style
mytheme = theme(panel.background = element_rect(fill='white', colour='black'),
                legend.position = "none",
                legend.key = element_blank(),
                #legend.background = element_blank(),   
                legend.box.background = element_blank(),
                panel.grid=element_blank(), 
                legend.title = element_text(size = 10),
                legend.text = element_text(size = 9),
                legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
                axis.ticks = element_line(color='black'),
                axis.line = element_line(colour = "black"), 
                axis.title.x = element_text(colour='black', size=13),
                axis.title.y = element_text(colour='black', size=13),
                axis.text = element_text(colour='black',size=11),
                plot.tag = element_text(size = 14, face = "bold")) 

# Loading field survey data
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)
#Field_group$Years <- as.factor(Field_group$Years)
Field_group$Site <- factor(Field_group$Site, levels = unique(Field_group$Site[order(Field_group$Latitude)]))
climate_data <- unique(Field_group[,c("Site", "Years", "Tave", "Precipitation")])
colnames(Field_group)

# calculate Environmental effects
Field_group$Effect_size <- log(Field_group$Fungal_field_Di/Field_group$Fungal_green_Di)

# Data Transformation
Field_group$RS <- sqrt(Field_group$RS)
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Phy_Di_log <- log10(Field_group$Phy_Di)
Field_group$Fun_Di_log <- log10(Field_group$Fun_Di)
#unique(Field_group$Species)
#total_data$Years = as.factor(total_data$Years)

# Consider normalizing your data
pd_attributes_variable <- attributes(scale(Field_group[c("Site_pool","Phy_Di","Fun_Di","Phy_Di_log","Fun_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation")]))
total_data = Field_group
colnames(total_data)
total_data[c("Site_pool","Phy_Di","Fun_Di","Phy_Di_log","Fun_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation")] = 
  scale(total_data[c("Site_pool","Phy_Di","Fun_Di","Phy_Di_log","Fun_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation")])

fm1 = lme(Effect_size ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Precipitation + Tave +
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


# Identified the best-fitting models for each response variable based on 
# Corrected Akaike Information Criterion (AICc)
Final_model <- get.models(dd12,1)[[1]] 
Final_model <- update(Final_model, ~., method = "REML") # Update to REML to extract estimates
summary(Final_model)
r.squaredGLMM(Final_model)

################################### Table S6 ###################################
Table_S6 <- as.data.frame(car::Anova(Final_model, type = 3))
Table_S6 <- Table_S6[-which(rownames(Table_S6) == "(Intercept)"),]
Table_S6 <- as.data.frame(Table_S6)
Table_S6$p.adj <- p.adjust(Table_S6$`Pr(>Chisq)`, method = "holm")
Table_S6$p.adj <- sprintf("%.3f", Table_S6$p.adj)
Table_S6$`Pr(>Chisq)` <- round(Table_S6$`Pr(>Chisq)`, 3)
Table_S6$Chisq <- round(Table_S6$Chisq, 2)
Table_S6$Parameter <- rownames(Table_S6)
Table_S6$VIF <- round(car::vif(Final_model), 2)
print(Table_S6[c(3,2,1,6,4,5,9,8,7), c(5,1:4,6)])

################################## Figure 3A ###################################
# Obtaining standardized regression coefficients and their 95% CI
MegaModelSummary <- as.data.frame(effectsize::effectsize(Final_model))[-1,]
MegaModelSummary$Parameter2 <- c("Funct-Dist", "Phylo-Dist", "Site pool", "Soil N", "Temperature", "Wcont",
                               "Funct-Dist × Soil N", "Funct-Dist × Temperature", "Phylo-Dist × Temperature")
MegaModelSummary$Parameter2 <- factor(MegaModelSummary$Parameter2, levels = rev(c("Site pool","Phylo-Dist", "Funct-Dist", "Wcont", "Soil N", "Temperature", 
                                                                                 "Phylo-Dist × Temperature", "Funct-Dist × Temperature", "Funct-Dist × Soil N")))
MegaModelSummary <- MegaModelSummary %>% left_join(Table_S6[,c("Parameter","p.adj")], by = "Parameter")

MegaModelSummary$Group <- c("Plant attributes", "Plant attributes", "Site pool", "Soil properties", "Climate", 
                            "Soil properties", "Interaction", "Interaction", "Interaction")
MegaModelSummary$Group <- factor(MegaModelSummary$Group, levels = c("Site pool", "Plant attributes", "Soil properties", "Climate", "Interaction"))

ggplot(MegaModelSummary, aes(x = Parameter2, y = Std_Coefficient, fill = Group))+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0.2, size = 0.8, color = "black")+
  geom_point(size = 3.5, pch = 21)+
  geom_segment(aes(y = 0, yend = 0, x = 0.5, xend = 9.3), color = "black", linetype = "dashed") + 
  geom_text(aes(y = CI_high, label = paste("italic(p)==", `p.adj`)),
            parse = TRUE, hjust = -0.4, vjust = 0.4, size = 3.5) + 
  #annotate("text", x = 9.5 , y = 0.08,label = "Best model: R2m = 0.30, R2c = 0.32", colour="black", size = 4) +  
  annotate("text", x = 9.7, y = 0.08,
           label = expression("Best model: " * italic(R)^2 * "m = 0.31, " * italic(R)^2 * "c = 0.33"), colour = "black", size = 4) + 
  annotate("text", x = 4.0, y = 0.12,
           label = expression(italic(p) * " < 0.001"), colour = "black", size = 3.5) + 
  labs(x = '', y = 'Standard regression coefficients', color = '', tag = "(a)") +  
  theme_classic() + coord_flip() +  
  scale_fill_manual(values = c("#949698","#DE7963", "#3F425A","#82B19D","#F0C986")) +
  theme(axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.x =  element_text(color = "black", size = 13),
        legend.text = element_text(size = 9, color = "black"),
        plot.margin = margin(0.5,1.5,0.5,1.5, unit = "cm"),
        legend.position = 'none',legend.title = element_blank(),legend.key = element_blank(),
        plot.tag = element_text(size = 14, face = "bold")) +
  scale_x_discrete(expand = expansion(mult = c(0.00, 0.11))) + 
  scale_shape_manual(values = c(16,21)) -> p3a; p3a

# Relative contribution of variables
hierarchical_data <- as.data.frame(glmm.hp(Final_model, type = "R2")$hierarchical.partitioning)
hierarchical_data$Group <- c("Plant attributes", "Plant attributes", "Site pool", "Soil properties", "Climate", 
                            "Soil properties", "Interaction", "Interaction", "Interaction")
hierarchical_data$Parameter2 <- c("Funct-Dist", "Phylo-Dist", "Site pool", "Soil N", "Temperature", "Wcont",
                                 "Funct-Dist × Soil N", "Funct-Dist × Temperature", "Phylo-Dist × Temperature")
hierarchical_data$Group <- factor(hierarchical_data$Group, levels = c("Site pool", "Plant attributes", "Soil properties", "Climate", "Interaction"))

hierarchical_data2 <- hierarchical_data %>% dplyr::group_by(Group) %>% 
  dplyr::summarise(`I.perc(%)` = sum(`I.perc(%)`))  %>% 
  mutate(label = paste0(round(`I.perc(%)`, 1), "%")) 

ggplot()+
  geom_bar(data = hierarchical_data2, aes(x = "", y = `I.perc(%)`, fill = Group), 
           stat = "identity", width = 0.5, color = "black") +
  geom_text(data = hierarchical_data2, 
            aes(x = "", y = `I.perc(%)`, label = label, group = Group), 
            position = position_stack(vjust = 0.5), color = "black", size = 4) +
  scale_y_continuous(expand = c(0, 0), position = "right")+
  scale_fill_manual(values = c("#949698","#DE7963", "#3F425A","#82B19D","#F0C986")) +
  theme_classic()+
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.y = element_text(color = "black", size = 13),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_line(linewidth = NA)) +
  labs(y = "Relative effect of estimates (%)") -> p3b; p3b


################################## Figure 3B ###################################
# Phylo-Dist × Temperature
eff_mod <- effect("Phy_Di_log:Tave", Final_model, xlevels = 5)
eff_mod_data <- data.frame(eff_mod)

# backtransform
eff_mod_data["Tave"] <- pd_attributes_variable$`scaled:center`["Tave"] + pd_attributes_variable$`scaled:scale`["Tave"]*eff_mod_data["Tave"]
eff_mod_data["Phy_Di_log"] <- pd_attributes_variable$`scaled:center`["Phy_Di_log"] + pd_attributes_variable$`scaled:scale`["Phy_Di_log"]*eff_mod_data["Phy_Di_log"]

total_data["Phy_Di_row"] <- pd_attributes_variable$`scaled:center`["Phy_Di_log"] + pd_attributes_variable$`scaled:scale`["Phy_Di_log"]*total_data["Phy_Di_log"]
total_data["Tave_row"] <- pd_attributes_variable$`scaled:center`["Tave"] + pd_attributes_variable$`scaled:scale`["Tave"]*total_data["Tave"]
eff_mod_data$Phy_Di_log <- round(eff_mod_data$Phy_Di_log, 2)

ggplot()+
  #geom_point(data = total_data, mapping = aes(x = Tave_row, y = Effect_size), pch = 21) + 
  geom_line(data = eff_mod_data, mapping = aes(Tave, fit, color = factor(Phy_Di_log)), size = 1.25) +
  labs(x = expression("Annual average temperature (°C)"), 
       y = bquote(atop("Environmental effects", 
                       Ln ~ "(" ~ frac(Fungi-dist["estimated in field survey"], 
                                       Fungi-dist["estimated in greenhouse experiment"]) ~ ")")),
       tag = "(b)", color = expression("Phylo Di(log"[10]*"(10)")) +
  geom_ribbon(data = eff_mod_data, mapping = aes(x = Tave, ymin=fit-se,ymax=fit+se,fill = factor(Phy_Di_log)), alpha = 0.3, colour = NA, show.legend = F) +
  scale_fill_manual(values = c("#184C3F","#6CB2AA", "#E4CB8F", "#BF8E43", "#57320F")) +
  scale_color_manual(values = c("#184C3F","#6CB2AA", "#E4CB8F", "#BF8E43", "#57320F"), name = expression("Phylo-Dist (log"[10]*")")) +
  annotate("text", label = expression(italic(p) == 0.021), x = 15, y = 0.25, size = 4) + 
  mytheme + theme(legend.position = c(0.8,0.28))-> P3b; P3b

# combination
(p3a+p3b) + plot_layout(widths = c(0.7,0.3)) -> P3a

(P3a|P3b) + plot_layout(widths = c(0.6,0.4))

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.