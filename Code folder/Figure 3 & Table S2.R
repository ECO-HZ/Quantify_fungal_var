################################################################################
############################ Figure 3 & Table S2 ###############################
################################################################################

# Loading the R packages
library(openxlsx)
library(car)
library(MuMIn)
library(ggplot2)
library(glmm.hp)
library(ggeffects)
library(patchwork)
library(ggtext)
library(dplyr)
library(lme4)
library(lmerTest)

# Custom style
mytheme = theme(panel.background = element_rect(fill='white', colour='black'),
                legend.position = "none",
                legend.key = element_blank(),
                legend.box.background = element_blank(),
                panel.grid=element_blank(), 
                legend.title = element_text(size = 11),
                legend.text = element_text(size = 10),
                legend.background = element_rect(fill = NA), 
                axis.ticks = element_line(color='black'),
                axis.line = element_line(colour = "black"), 
                axis.title.x = element_text(colour='black', size=13),
                axis.title.y = element_text(colour='black', size=13),
                axis.text = element_text(colour='black',size=11),
                plot.tag = element_text(size = 14, face = "bold")) 

# Loading field survey data
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)
Field_group$Years <- as.factor(Field_group$Years)
Field_group$Site <- factor(Field_group$Site, levels = unique(Field_group$Site[order(Field_group$Latitude)]))

# Calculate environmental effects
Field_group$Effect_size <- log(Field_group$Fungal_Di_field_all/Field_group$Fungal_Di_green_all)

# Data Transformation
Field_group$RS <- sqrt(Field_group$RS)
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont*100)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Funct_Di_log <- log10(Field_group$Funct_Di)
Field_group$Phylo_Di_log <- log10(Field_group$Phylo_Di)

# ==============================================================================
# STEP 1: DECOMPOSE MACRO-CLIMATE VARIABLES (Tave & Precipitation)
# As requested by Reviewer 1 (temp = grand_mean + site_effect + year_effect + residual)
# ==============================================================================

# Calculate Grand Means
grand_Tave   <- mean(Field_group$Tave, na.rm = TRUE)
grand_Precip <- mean(Field_group$Precipitation, na.rm = TRUE)

# Calculate Site Means (Site Effect, 5 df)
site_means <- Field_group %>%
  group_by(Site) %>%
  summarise(Tave_site = mean(Tave, na.rm = TRUE),
            Precip_site = mean(Precipitation, na.rm = TRUE),
            .groups = "drop")

# Calculate Year Means (Year Effect, 2 df)
year_means <- Field_group %>%
  group_by(Years) %>%
  summarise(Tave_year = mean(Tave, na.rm = TRUE),
            Precip_year = mean(Precipitation, na.rm = TRUE),
            .groups = "drop")

# Merge decomposed components back into Field_group
Field_group <- Field_group %>%
  left_join(site_means, by = "Site") %>%
  left_join(year_means, by = "Years") %>%
  mutate(
    # Site x Year Residual Deviation (10 df)
    Tave_dev   = Tave - Tave_site - Tave_year + grand_Tave,
    Precip_dev = Precipitation - Precip_site - Precip_year + grand_Precip)

unique(Field_group$Tave_dev)
unique(Field_group$Precip_dev)

Field_group_cor = Field_group[,c("Funct_Di_log","Phylo_Di_log",
                                 "Tave_site","Tave_year","Tave_dev",
                                 "Precip_site","Precip_year","Precip_dev", 
                                 "Soil_N", "Soil_ph", "Wcont")]
colnames(Field_group_cor) = c("Funct-Dist","Phylo-Dist","Spatial temperature","Interannual temperature","Temperature anomaly", 
                              "Spatial precipitation " ,"Interannual precipitation", "Precipitation anomaly",
                              "Soil N","Soil pH", "Wcont")

# ==============================================================================
# STEP 2: NORMALIZE DATA
# ==============================================================================
var_select <- c("Effect_size","Site_pool","Phylo_Di","Funct_Di","Phylo_Di_log","Funct_Di_log",
                "Soil_ph", "Wcont","Soil_N",
                "Tave_site", "Tave_year", "Tave_dev",
                "Precip_site", "Precip_year", "Precip_dev",
                "PCoA1")

pd_attributes_variable <- attributes(scale(Field_group[var_select]))
total_data <- Field_group
total_data[var_select] <- scale(total_data[var_select])

# ==============================================================================
# STEP 3: GLOBAL LINEAR MIXED MODEL WITH DECOMPOSED CLIMATE TERMS
# ==============================================================================

fm1 <- lme4::lmer(Effect_size ~ PCoA1 + Phylo_Di_log + Funct_Di_log + Soil_ph + Wcont + Soil_N + 
                    
                    # Decomposed Climate Terms
                    Tave_site + Tave_year + Tave_dev +
                    Precip_site + Precip_year + Precip_dev +
                    
                    # Interactions with Plant Traits
                    Phylo_Di_log:Tave_site + Phylo_Di_log:Tave_dev + 
                    Phylo_Di_log:Precip_site + Phylo_Di_log:Precip_dev + 
                    Phylo_Di_log:Soil_ph + Phylo_Di_log:Wcont + Phylo_Di_log:Soil_N + 
                    
                    Funct_Di_log:Tave_site + Funct_Di_log:Tave_dev + 
                    Funct_Di_log:Precip_site + Funct_Di_log:Precip_dev + 
                    Funct_Di_log:Soil_ph + Funct_Di_log:Wcont + Funct_Di_log:Soil_N + 
                    
                    # Explicit Design Random Factors (Site, Year, Site x Year)
                    (1|Site) + (1|Years) + (1|Site:Years), 
                  data = total_data, REML = FALSE)

summary(fm1)
vif(fm1)
ranova(fm1)

# The initial model was fitted using maximum likelihood and simplified through backward elimination based on likelihood-ratio tests
drop1(fm1, test = "Chi") 
f.sbs1 <- update(fm1, ~. -Funct_Di_log:Precip_site)
drop1(f.sbs1, test = "Chi") 
f.sbs2 <- update(f.sbs1, ~. -Phylo_Di_log:Soil_ph)
drop1(f.sbs2, test = "Chi") 
f.sbs3 <- update(f.sbs2, ~. -Phylo_Di_log:Precip_dev)
drop1(f.sbs3, test = "Chi") 
f.sbs4 <- update(f.sbs3, ~. -Phylo_Di_log:Wcont)
drop1(f.sbs4, test = "Chi") 
f.sbs5 <- update(f.sbs4, ~. -Phylo_Di_log:Soil_N)
drop1(f.sbs5, test = "Chi") 
f.sbs6 <- update(f.sbs5, ~. -Phylo_Di_log:Precip_site)
drop1(f.sbs6, test = "Chi") 
f.sbs7 <- update(f.sbs6, ~. -Funct_Di_log:Tave_site)
drop1(f.sbs7, test = "Chi") 
f.sbs8 <- update(f.sbs7, ~. -Precip_site)
drop1(f.sbs8, test = "Chi") 
f.sbs9 <- update(f.sbs8, ~. -Funct_Di_log:Precip_dev)
drop1(f.sbs9, test = "Chi") 
f.sbs10 <- update(f.sbs9, ~. -Precip_year)
drop1(f.sbs10, test = "Chi") 
f.sbs11 <- update(f.sbs10, ~. -Phylo_Di_log:Tave_dev)
drop1(f.sbs11, test = "Chi") 
f.sbs12 <- update(f.sbs11, ~. -Funct_Di_log:Soil_ph)
drop1(f.sbs12, test = "Chi") 
f.sbs13 <- update(f.sbs12, ~. -Soil_ph)
drop1(f.sbs13, test = "Chi") 
f.sbs14 <- update(f.sbs13, ~. -Funct_Di_log:Wcont)
drop1(f.sbs14, test = "Chi") 
f.sbs15 <- update(f.sbs14, ~. -Funct_Di_log:Tave_dev)
drop1(f.sbs15, test = "Chi") 
f.sbs16 <- update(f.sbs15, ~. -Wcont)
drop1(f.sbs16, test = "Chi") 
AIC(f.sbs16)

f.sbs.fin <- update(f.sbs16, REML = TRUE) # Update to REML to extract estimates
fm1_test <- as_lmerModLmerTest(f.sbs.fin)
vif(fm1_test)
anova(fm1_test, ddf = "Kenward-Roger")


Table_S2 <- as.data.frame(anova(fm1_test, ddf = "Kenward-Roger"))
#Table_S2 <- Table_S4[-which(rownames(Table_S2) == "Residuals"),]
Table_S2$`Pr(>F)` <- round(Table_S2$`Pr(>F)`, 3)
Table_S2$`F` <- round(Table_S2$`F value`, 2)
Table_S2$Parameter <- rownames(Table_S2)
Table_S2$VIF <- round(car::vif(fm1_test), 2)

################################## Figure 4a ###################################
# Obtaining standardized regression coefficients and their 95% CI
Final_model = lmer(Effect_size ~ PCoA1 + Tave_site + Tave_year + Tave_dev + Precip_dev + Soil_N + Phylo_Di_log + Funct_Di_log + 
                     Funct_Di_log:Soil_N + Phylo_Di_log:Tave_site + (1|Site) + (1|Years) + (1|Site:Years), data = total_data)
model_aov_results = anova(Final_model, ddf = "Kenward-Roger")
summary(Final_model)
MuMIn::r.squaredGLMM(Final_model)

as.data.frame(vif(Final_model))

model_aov_results_df = as.data.frame(model_aov_results)

model_aov_results_df$Parameter = rownames(model_aov_results_df)

model_aov_results_df <- model_aov_results_df %>%
  mutate(Parameter = recode(as.character(Parameter),
                            "Soil_N:Funct_Di_log" = "Funct_Di_log:Soil_N",
                            "Tave_site:Phylo_Di_log" = "Phylo_Di_log:Tave_site"))

summary_Data <- as.data.frame(summary(Final_model))

summary_Data <- as.data.frame(summary(Final_model)$coefficients)
summary_Data$Parameter <- rownames(summary_Data)
summary_Data = summary_Data[-1, ]

summary_Data <- summary_Data %>%
  mutate(Parameter = recode(as.character(Parameter),
                            "Soil_N:Funct_Di_log" = "Funct_Di_log:Soil_N",
                            "Tave_site:Phylo_Di_log" = "Phylo_Di_log:Tave_site"))

################################################################################
MegaModelSummary <- as.data.frame(effectsize::effectsize(Final_model))[-1,]

MegaModelSummary <- MegaModelSummary %>%
  mutate(Parameter = recode(as.character(Parameter),
                            "Soil_N:Funct_Di_log" = "Funct_Di_log:Soil_N",
                            "Tave_site:Phylo_Di_log" = "Phylo_Di_log:Tave_site"))

# Relative contribution of variables
hierarchical_data <- glmm.hp::glmm.hp(Final_model, type = "R2")$hierarchical.partitioning
hierarchical_data_df = as.data.frame(hierarchical_data)
hierarchical_data_df$Parameter = rownames(hierarchical_data_df)
print(hierarchical_data_df)

MegaModelSummary_all = MegaModelSummary %>% left_join(hierarchical_data_df) %>%
  left_join(model_aov_results_df[,c("Parameter", "Pr(>F)")]) %>%
  left_join(summary_Data[,c("Parameter", "Estimate", "Std. Error")])

################################## Figure 3a ###################################
# Obtaining standardized regression coefficients and their 95% CI
MegaModelSummary_all$Term_display = c("Field fungal composition", "Spatial temperature", "Interannual temperature", "Temperature anomaly",
                                      "Precipitation anomaly", "Soil N", "Phylo-Dist", "Funct-Dist", 
                                      "Funct-Dist × Soil N", "Phylo-Dist × Spatial temperature")

MegaModelSummary_all$Term_display = factor(MegaModelSummary_all$Term_display, levels = rev(unique(MegaModelSummary_all$Term_display)))

# add group information
MegaModelSummary_all$Group <- c("Field com", rep("Climate", 4), rep("Soil properties", 1), 
                                rep("Plant attributes", 2), rep("Interaction", 2))

MegaModelSummary_all$Group <- factor(MegaModelSummary_all$Group, levels = unique(MegaModelSummary_all$Group))

ggplot(MegaModelSummary_all, aes(x = Term_display, y = Std_Coefficient, fill = Group))+
  geom_hline(yintercept = 0, linetype = 1, color = "grey") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width=0, size = 0.8, color = "black")+
  geom_point(size = 3.5, pch = 21) +
  #geom_segment(aes(y = 0, yend = 0, x = 0.5, xend = 12.3), color = "black", linetype = "dashed") + 
  geom_text(aes(y = CI_high, label = paste("italic(p)==", round(`Pr(>F)`, 3))),
            parse = TRUE, hjust = -0.4, vjust = 0.4, size = 4) + 
  labs(x = NULL, 
       y = 'Parameter estimates', 
       #title = "Best model: <i>R</i><sup>2</sup> = 0.320", 
       tag = "(a)") +  
  theme_classic() + coord_flip() +  
  scale_fill_manual(values = c("Field com" = "#BC5546", "Climate" = "#2F4590", "Soil properties" = "#6EA3C5",
                               "Plant attributes" = "#8E333A", "Interaction" = "#EFA961")) +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title =  element_text(color = "black", size = 14),
        legend.text = element_text(size = 9, color = "black"),
        plot.title = element_textbox(size = 12, color = "black", fill = "white",     
                                     box.color = "black", width = grid::unit(1, "npc"),padding = margin(5, 5, 5, 5),  
                                     margin = margin(b = 5), halign = 0.5,linetype = "solid"),
        plot.margin = margin(0.5,1.5,0.5,1.5, unit = "cm"),
        legend.position = 'none',
        plot.tag = element_text(size = 14, face = "bold")) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) + 
  scale_shape_manual(values = c(16,21)) -> Figure_3a1; Figure_3a1

################################################################################
MegaModelSummary_deal2 = MegaModelSummary_all %>% group_by(Group) %>%
  summarise(explained_all = sum(`I.perc(%)`))

MegaModelSummary_deal2$explained_all2 = (MegaModelSummary_deal2$explained_all)/sum(MegaModelSummary_deal2$explained_all)*100

ggplot(MegaModelSummary_deal2, aes(x = "Importance", y = explained_all2, fill = Group, color = Group)) +
  geom_col(width = 1) +
  #scale_fill_viridis(option = "D", direction = -1) + # + 
  theme_classic()+ 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), position = "right") +
  scale_color_manual(values = c("Field com" = "#BC5546", "Climate" = "#2F4590", "Soil properties" = "#6EA3C5",
                                "Plant attributes" = "#8E333A", "Interaction" = "#EFA961")) +
  scale_fill_manual(values = c("Field com" = "#BC5546", "Climate" = "#2F4590", "Soil properties" = "#6EA3C5",
                               "Plant attributes" = "#8E333A", "Interaction" = "#EFA961")) +
  theme(panel.grid = element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), 
        axis.title =  element_text(color = "black", size = 14),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none") + 
  labs(x = '', y = "Relative effect of estimates (%)") -> Figure_3a2; Figure_3a2

# 9.11 x 10.10
Figure_3a1+Figure_3a2 + plot_layout(widths = c(0.9,0.1)) -> Figure_3a; Figure_3a
#ggsave("Figure_3a.pdf", plot = Figure_3a, width = 9.11, height = 10.10, units = "in", dpi = 300)

################################## Figure 3b ###################################
pred_mode <- ggeffect(Final_model, terms = c("Soil_N","Funct_Di_log"))
eff_mod_data <- data.frame(pred_mode)
colnames(eff_mod_data)[1] <- "Soil_N"
colnames(eff_mod_data)[2] <- "Effect_size"
colnames(eff_mod_data)[6] <- "Funct_Di_log"

# back transform attributes_variable
eff_mod_data["Soil_N"] <- pd_attributes_variable$`scaled:center`["Soil_N"] + 
  pd_attributes_variable$`scaled:scale`["Soil_N"]*eff_mod_data["Soil_N"]

eff_mod_data["Effect_size"] <- pd_attributes_variable$`scaled:center`["Effect_size"] + 
  pd_attributes_variable$`scaled:scale`["Effect_size"]*eff_mod_data["Effect_size"]

eff_mod_data$Funct_Di_log <- ifelse(eff_mod_data$Funct_Di_log == "-1", "Low Funct−Dist (-1 SD)", 
                                    ifelse(eff_mod_data$Funct_Di_log == "0", "Mean Funct−Dist", "High Funct−Dist (+1 SD)"))
eff_mod_data$Funct_Di_log <- factor(eff_mod_data$Funct_Di_log, levels = c("Low Funct−Dist (-1 SD)", "Mean Funct−Dist", "High Funct−Dist (+1 SD)"))

ggplot()+
  geom_line(data = eff_mod_data, mapping = aes(Soil_N, Effect_size, color = factor(Funct_Di_log)), size = 1.25) +
  labs(x = expression("Soil total nitrogen content (%, sqrt)"), 
       y = bquote(atop("Environmental effects", 
                       Ln ~ "(" ~ frac(Fungi-dist[" estimated in field"], 
                                       Fungi-dist[" estimated in greenhouse"]) ~ ")")),
       tag = "(b)", color = expression("Phylo Di(log"[10]*"(10)")) +
  geom_hline(yintercept = 0, linetype = 1, color = "grey") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_fill_manual(values = c("#184C3F", "#E4CB8F", "#57320F")) +
  scale_color_manual(values = c("#184C3F", "#E4CB8F", "#57320F"), name = "Funct-Dist") +
  annotate("text", label = expression(italic(p) == 0.049), x = 0.6, y = -0.12, size = 4) + 
  mytheme + theme(legend.position = c(0.4,0.80)) -> Figure_3b; Figure_3b

################################## Figure 4c ###################################
pred_mode <- ggeffect(Final_model, terms = c("Tave_site","Phylo_Di_log"))
eff_mod_data <- data.frame(pred_mode)
colnames(eff_mod_data)[1] = "Tave_site"
colnames(eff_mod_data)[2] <- "Effect_size"
colnames(eff_mod_data)[6] = "Phylo_Di_log"

# back transform attributes_variable
eff_mod_data["Tave_site"] <- pd_attributes_variable$`scaled:center`["Tave_site"] + 
  pd_attributes_variable$`scaled:scale`["Tave_site"]*eff_mod_data["Tave_site"]

eff_mod_data["Effect_size"] <- pd_attributes_variable$`scaled:center`["Effect_size"] + 
  pd_attributes_variable$`scaled:scale`["Effect_size"]*eff_mod_data["Effect_size"]

eff_mod_data$Phylo_Di_log <- ifelse(eff_mod_data$Phylo_Di_log == "-1", "Low Phylo−Dist (- 1 SD)", 
                                    ifelse(eff_mod_data$Phylo_Di_log == "0", "Mean Phylo−Dist", "High Phylo−Dist (+ 1 SD)"))
eff_mod_data$Phylo_Di_log <- factor(eff_mod_data$Phylo_Di_log, levels = c("Low Phylo−Dist (- 1 SD)", "Mean Phylo−Dist", "High Phylo−Dist (+ 1 SD)"))


ggplot()+
  geom_line(data = eff_mod_data, mapping = aes(Tave_site, Effect_size, color = factor( Phylo_Di_log)), size = 1.25) +
  labs(x = expression("Spatial temperature (°C)"), 
       y = bquote(atop("Environmental effects", 
                       Ln ~ "(" ~ frac(Fungi-dist[" estimated in field"], 
                                       Fungi-dist[" estimated in greenhouse"]) ~ ")")),
       tag = "(c)", color = expression("Phylo Di(log"[10]*"(10)")) +
  geom_hline(yintercept = 0, linetype = 1, color = "grey") +
  scale_fill_manual(values = c("#184C3F", "#E4CB8F", "#57320F")) +
  scale_color_manual(values = c("#184C3F", "#E4CB8F", "#57320F"), name = "Phylo-Dist") +
  annotate("text", label = expression(italic(p) == 0.026), x = 21, y = -0.03, size = 4) + 
  mytheme + theme(legend.position = c(0.4,0.80)) -> Figure_3c; Figure_3c

# 
Figure_3a
Figure_3b/Figure_3c -> Figure_3_right
#ggsave("Figure_3a_0730.pdf", plot = Figure_3a, width = 8.5, height = 8, units = "in", dpi = 300)
#ggsave("Figure_3b_0730.pdf", plot = Figure_3_right, width = 4.5, height = 7.8, units = "in", dpi = 300)

