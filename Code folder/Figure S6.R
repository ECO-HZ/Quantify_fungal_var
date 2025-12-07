################################################################################
################################## Figure S6 ###################################
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

# Custom style
mytheme <- theme_bw() + 
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none",
        legend.key = element_blank(),
        panel.grid=element_blank(), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11.38),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_text(colour='black', size=13),
        axis.title.y = element_text(colour='black', size=13),
        axis.text = element_text(colour='black',size=11),
        plot.tag = element_text(size = 14, face = "bold"),
        plot.title = element_textbox(
          size = 14, color = "black", fill = "grey90",
          box.color = "grey50",padding = margin(5, 5, 5, 5), margin = margin(b = 0),       
          halign = 0.5, width = grid::unit(1, "npc"))) 

################################# Figure S6a ###################################
# Loading field survey data
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

Field_group$Effect_size_all <- log(Field_group$Fungal_Di_field_all/Field_group$Fungal_Di_green_all)
Field_group$Effect_size_com <- log(Field_group$Fungal_Di_field_com/Field_group$Fungal_Di_green_com)

cor.test(Field_group$Effect_size_all, Field_group$Effect_size_com, method = "spearman")

mod <- lm(Effect_size_all ~ Effect_size_com, data = Field_group)
anova(mod)
AIC(mod)
summary(mod)

ggplot(Field_group, aes(Effect_size_com, Effect_size_all))+
  geom_point(size = 2.5, pch = 21, fill = "grey") + 
  geom_abline(intercept=0,slope=1 , linetype = 1, color = "#96383E", size = 0.8)+
  geom_smooth(data = Field_group, aes(x =Effect_size_com , y = Effect_size_all), 
              method = "lm", formula = y ~ x, color = "black", se = F, linetype = 1, size = 0.8) + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(-0.8,0.35)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(-0.8,0.35)) + 
  theme_bw() + mytheme + 
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm")) + 
  annotate("segment", x = -0.75, xend = -0.68, y = 0.35, yend = 0.35, size = 1, color = "#96383E") + 
  annotate("text", x = -0.60, y = 0.35, label = "1:1 line", size = 4) +
  annotate("segment", x = -0.75, xend = -0.68, y = 0.30, yend = 0.30, size = 1, color = "black") + 
  annotate("text", x = -0.58, y = 0.30, label = "Model fit", size = 4) + 
  labs(x = "Environmental effects estimated\nbased on common fungal taxa",
       y = "Environmental effects estimated\nbased on complete fungal taxa", 
       tag = "(a)") -> Figure_S6a; Figure_S6a

################################# Figure S6b ###################################
# Loading field survey data
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)
Field_group$Years <- as.factor(Field_group$Years)
Field_group$Site <- factor(Field_group$Site, levels = unique(Field_group$Site[order(Field_group$Latitude)]))

# calculate environmental effects
Field_group$Effect_size <- log(Field_group$Fungal_Di_field_com/Field_group$Fungal_Di_green_com)

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

# performed the global linear regression model
fm1 <- lm(Effect_size ~ PCoA1 + Phylo_Di_log + Funct_Di_log + Soil_ph + Wcont + Soil_N + Precipitation + Tave +
            Phylo_Di_log:Tave + Phylo_Di_log:Precipitation + Phylo_Di_log:Soil_ph + Phylo_Di_log:Wcont + Phylo_Di_log:Soil_N + 
            Funct_Di_log:Tave + Funct_Di_log:Precipitation + Funct_Di_log:Soil_ph + Funct_Di_log:Wcont + Funct_Di_log:Soil_N, data = total_data)

car::vif(fm1)

options(na.action = "na.fail")
dd12 <- dredge(fm1, subset = ~ PCoA1 &
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

# Identified the best-fitting models for each response variable based on 
# Corrected Akaike Information Criterion (AICc)
Final_model <- get.models(dd12,1)[[1]] 
performance::r2(Final_model)
AICc(Final_model)

################################ Table of model ################################
Table_mod <- as.data.frame(Anova(Final_model))
Table_mod <- Table_mod[-which(rownames(Table_mod) == "Residuals"),]
Table_mod$p.adj <- p.adjust(Table_mod$`Pr(>F)`, method = "holm")
Table_mod$p.adj <- sprintf("%.3f", Table_mod$p.adj)
Table_mod$`Pr(>F)` <- round(Table_mod$`Pr(>F)`, 3)
Table_mod$`F` <- round(Table_mod$`F value`, 2)
Table_mod$Parameter <- rownames(Table_mod)
Table_mod$VIF <- round(car::vif(Final_model), 2)
Table_mod$df <- paste0("1", ",", "374")
rownames(Table_mod) <- c("Funct-Dist", "Field fungal composition", "Phylo-Dist", "Soil N", "Soil pH", "Temperature", "Wcont",
                        "Funct-Dist × Soil N", "Funct-Dist × Wcont", "Phylo-Dist × Temperature")
print(Table_mod[c(2,1,3,6,5,4,7,8,9,10) ,c("df","F","Pr(>F)","p.adj","VIF")])


################################## Figure 4a ###################################
# Obtaining standardized regression coefficients and their 95% CI
Final_model <- lm(Effect_size ~ Funct_Di_log + PCoA1 + Phylo_Di_log + 
                    Soil_N + Soil_ph + Tave + Wcont + Funct_Di_log:Soil_N + Funct_Di_log:Wcont + 
                    Phylo_Di_log:Tave, data = total_data)

MegaModelSummary <- as.data.frame(effectsize::effectsize(Final_model))[-1,]
MegaModelSummary$Parameter2 <- c("Funct-Dist", "Field fungal composition", "Phylo-Dist", "Soil N", "Soil pH", "Temperature", "Wcont",
                                 "Funct-Dist × Soil N", "Funct-Dist × Wcont", "Phylo-Dist × Temperature")
MegaModelSummary$Parameter2 <- factor(MegaModelSummary$Parameter2, levels = rev(c("Field fungal composition", "Funct-Dist", "Phylo-Dist", "Temperature", "Soil pH", "Soil N", "Wcont", 
                                                                                  "Funct-Dist × Soil N", "Funct-Dist × Wcont", "Phylo-Dist × Temperature")))
MegaModelSummary <- MegaModelSummary %>% left_join(Table_mod[,c("Parameter","p.adj")], by = "Parameter")

MegaModelSummary$Group <- c("Plant attributes", "Field com", "Plant attributes", "Soil properties", "Soil properties", "Climate", "Soil properties",
                            "Interaction", "Interaction", "Interaction")
MegaModelSummary$Group <- factor(MegaModelSummary$Group, levels = c("Field com", "Plant attributes", "Climate", "Soil properties", "Interaction"))

ggplot(MegaModelSummary, aes(x = Parameter2, y = Std_Coefficient, fill = Group))+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0.2, size = 0.8, color = "black")+
  geom_point(size = 3.5, pch = 21)+
  geom_segment(aes(y = 0, yend = 0, x = 0.5, xend = 10.3), color = "black", linetype = "dashed") + 
  geom_text(aes(y = CI_high, label = paste("italic(p)==", `p.adj`)),
            parse = TRUE, hjust = -0.4, vjust = 0.4, size = 4) + 
  labs(x = '', 
       y = 'Standard regression coefficients', 
       title = "Best model: <i>R</i><sup>2</sup> = 0.262", tag = "(a)") +  
  theme_classic() + coord_flip() +  
  scale_fill_manual(values = c("#436C88","#49A3A4","#80B09A","#C2887D","#F0C986")) +
  #scale_fill_manual(values = c("#408EA8","#58C5BF","#FFD15C","#FF8942")) +
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
  scale_x_discrete(expand = expansion(mult = c(0.00, 0.05))) + 
  scale_shape_manual(values = c(16,21)) -> Figure_S6b1; Figure_S6b1

# Relative contribution of variables
hierarchical_data <- as.data.frame(glmm.hp::glmm.hp(Final_model, type = "R2")$hierarchical.partitioning)
hierarchical_data$Group <- c("Plant attributes", "Field com", "Plant attributes", "Soil properties", "Soil properties", "Climate", "Soil properties",
                             "Interaction", "Interaction", "Interaction")
hierarchical_data$Parameter2 <- c("Funct-Dist", "Field fungal composition", "Phylo-Dist", "Soil N", "Soil pH", "Temperature", "Wcont",
                                  "Funct-Dist × Soil N", "Funct-Dist × Wcont", "Phylo-Dist × Temperature")
hierarchical_data$Group <- factor(hierarchical_data$Group, levels = c("Field com", "Plant attributes", "Climate", "Soil properties", "Interaction"))

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
  scale_fill_manual(values = c("#436C88","#49A3A4","#80B09A","#C2887D","#F0C986")) +
  #scale_fill_manual(values = c("#408EA8","#58C5BF","#FFD15C","#FF8942")) +
  theme_classic()+
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_line(linewidth = NA)) +
  labs(y = "Relative effect of estimates (%)") -> Figure_S6b2; Figure_S6b2

(Figure_S6b1|Figure_S6b2) + plot_layout(widths = c(0.7,0.3)) -> Figure_S6b; Figure_S6b

################################## Figure S6c ##################################
pred_mode <- ggeffect(Final_model, terms = c("Soil_N","Funct_Di_log"))
eff_mod_data <- data.frame(pred_mode)
colnames(eff_mod_data)[1] <- "Soil_N"
colnames(eff_mod_data)[6] <- "Funct_Di_log"

# back transform attributes_variable
eff_mod_data["Soil_N"] <- pd_attributes_variable$`scaled:center`["Soil_N"] + 
  pd_attributes_variable$`scaled:scale`["Soil_N"]*eff_mod_data["Soil_N"]

eff_mod_data$Funct_Di_log <- ifelse(eff_mod_data$Funct_Di_log == "-1", "- 1 SD", 
                                    ifelse(eff_mod_data$Funct_Di_log == "0", "Mean", "+ 1 SD"))
eff_mod_data$Funct_Di_log <- factor(eff_mod_data$Funct_Di_log, levels = c("- 1 SD", "Mean", "+ 1 SD"))

ggplot()+
  geom_line(data = eff_mod_data, mapping = aes(Soil_N, predicted, color = factor(Funct_Di_log)), size = 1.25) +
  labs(x = expression("Soil total nitrogen content (%, sqrt)"), 
       y = bquote(atop("Environmental effects", 
                       Ln ~ "(" ~ frac(Fungi-dist["estimated in field survey"], 
                                       Fungi-dist["estimated in greenhouse experiment"]) ~ ")")),
       tag = "(c)") +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = c("#184C3F", "#E4CB8F", "#57320F")) +
  scale_color_manual(values = c("#184C3F", "#E4CB8F", "#57320F"), name = "Funct-Dist") +
  annotate("text", label = expression(italic(p) == 0.021), x = 0.2, y = -0.02, size = 4) + 
  mytheme + theme(legend.position = c(0.80,0.20)) -> Figure_S6c; Figure_S6c


################################## Figure S6d ##################################
pred_mode <- ggeffect(Final_model, terms = c("Tave","Phylo_Di_log"))
eff_mod_data <- data.frame(pred_mode)
colnames(eff_mod_data)[1] = "Tave"
colnames(eff_mod_data)[6] = "Phylo_Di_log"

# back transform attributes_variable
eff_mod_data["Tave"] <- pd_attributes_variable$`scaled:center`["Tave"] + 
  pd_attributes_variable$`scaled:scale`["Tave"]*eff_mod_data["Tave"]

eff_mod_data$Phylo_Di_log <- ifelse(eff_mod_data$Phylo_Di_log == "-1", "- 1 SD", 
                                    ifelse(eff_mod_data$Phylo_Di_log == "0", "Mean", "+ 1 SD"))
eff_mod_data$Phylo_Di_log <- factor(eff_mod_data$Phylo_Di_log, levels = c("- 1 SD", "Mean", "+ 1 SD"))


ggplot()+
  geom_line(data = eff_mod_data, mapping = aes(Tave, predicted, color = factor( Phylo_Di_log)), size = 1.25) +
  labs(x = expression("Annual average temperature (°C)"), 
       y = bquote(atop("Environmental effects", 
                       Ln ~ "(" ~ frac(Fungi-dist["estimated in field survey"], 
                                       Fungi-dist["estimated in greenhouse experiment"]) ~ ")")),
       tag = "(d)") +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = c("#184C3F", "#E4CB8F", "#57320F")) +
  scale_color_manual(values = c("#184C3F", "#E4CB8F", "#57320F"), name = "Phylo-Dist") +
  annotate("text", label = expression(italic(p) == 0.021), x = 15, y = 0.02, size = 4) + 
  mytheme + theme(legend.position = c(0.80,0.20)) -> Figure_S6d; Figure_S6d

