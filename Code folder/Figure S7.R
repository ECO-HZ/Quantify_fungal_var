################################################################################
################################## Figure S7 ###################################
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
mytheme <- theme_classic() + 
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

################################# Figure S7a ###################################
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
  annotate("text", x = 0, y = -0.31, label = "F[1*','*383] == 1397.7", parse = TRUE, size = 4) + 
  annotate("text", x = 0, y = -0.38, label = "italic(R)^2 == 0.785", parse = TRUE, size = 4) + 
  annotate("text", x = 0, y = -0.45, label = "italic(p) < 0.001", parse = TRUE, size = 4) +
  labs(x = "Environmental effects estimated\nbased on common fungal taxa",
       y = "Environmental effects estimated\nbased on complete fungal taxa", 
       tag = "(a)") -> Figure_S7a; Figure_S7a

################################# Figure S7b ###################################
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
var_select <- c("Effect_size","Site_pool","Phylo_Di","Funct_Di","Phylo_Di_log","Funct_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation","PCoA1")
pd_attributes_variable <- attributes(scale(Field_group[var_select]))
total_data <- Field_group
total_data[var_select] <- scale(total_data[var_select])

# performed the global linear regression model
fm1 <- lm(Effect_size ~ PCoA1 + Phylo_Di_log + Funct_Di_log + Soil_ph + Wcont + Soil_N + Precipitation + Tave +
            Phylo_Di_log:Tave + Phylo_Di_log:Precipitation + Phylo_Di_log:Soil_ph + Phylo_Di_log:Wcont + Phylo_Di_log:Soil_N + 
            Funct_Di_log:Tave + Funct_Di_log:Precipitation + Funct_Di_log:Soil_ph + Funct_Di_log:Wcont + Funct_Di_log:Soil_N, data = total_data)

car::vif(fm1)
performance::r2(fm1)

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

#  the proportion of the best-fitting models a factor was selected)
de6 <- model.avg(dd12, subset = delta < 2, fit = TRUE)
de6$msTable[1,]

# calculated the possibility (i.e., the proportion of the best-fitting models a factor was selected) 
# of having an effect for each predictor over the multimodel space composed of all fitting models (ΔAICc < 2)
sw_unclass <- unclass(sw(de6))
weights <- sw_unclass
n_models <- attr(sw_unclass, "n.models")
variables <- names(sw_unclass)
importance_df <- data.frame(
  Parameter = variables,
  Sum_of_weights = as.numeric(weights),
  N_containing_models = as.numeric(n_models))
importance_df <- importance_df[order(-importance_df$Sum_of_weights), ]
importance_df$sel_prop <- importance_df$N_containing_models/6
print(importance_df)

# Identified the best-fitting models for each response variable based on 
# Corrected Akaike Information Criterion (AICc) and weight
Final_model <- get.models(de6,1)[[1]] 
performance::r2(Final_model)
AICc(Final_model)
summary(Final_model)

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

# Obtaining standardized regression coefficients and their 95% CI
MegaModelSummary <- as.data.frame(effectsize::effectsize(Final_model))[-1,]

# add other information (e.g., importance of predictors)
MegaModelSummary <- MegaModelSummary %>% 
  left_join(Table_mod[,c("Parameter","p.adj")], by = "Parameter") %>%
  left_join(importance_df)
head(MegaModelSummary)

# rename predictors
MegaModelSummary_deal <- MegaModelSummary %>%
  mutate(Term_display = case_when(
    Parameter == "PCoA1" ~ "Field fungal composition",
    Parameter == "Funct_Di_log" ~ "Funct-Dist",
    Parameter == "Funct_Di_log:Soil_ph" ~ "Funct-Dist × Soil pH",
    Parameter == "Funct_Di_log:Tave" ~ "Funct-Dist × Temperature",
    Parameter == "Funct_Di_log:Wcont" ~ "Funct-Dist × Wcont",
    Parameter == "Phylo_Di_log" ~ "Phylo-Dist",
    Parameter == "Phylo_Di_log:Precipitation" ~ "Phylo-Dist × Precipitation",
    Parameter == "Phylo_Di_log:Soil_N" ~ "Phylo-Dist × Soil N",
    Parameter == "Phylo_Di_log:Soil_ph" ~ "Phylo-Dist × Soil pH",
    Parameter == "Phylo_Di_log:Tave" ~ "Phylo-Dist × Temperature",
    Parameter == "Precipitation" ~ "Precipitation",
    Parameter == "Soil_N" ~ "Soil N",
    Parameter == "Soil_ph" ~ "Soil pH",
    Parameter == "Tave" ~ "Temperature",
    Parameter == "Wcont" ~ "Wcont",
    Parameter == "Funct_Di_log:Precipitation" ~ "Funct-Dist × Precipitation",
    Parameter == "Funct_Di_log:Soil_N" ~ "Funct-Dist × Soil N",
    Parameter == "Phylo_Di_log:Wcont" ~ "Phylo-Dist × Wcont",
    TRUE ~ Parameter))

# reorder
order = c("PCoA1","Funct_Di_log","Phylo_Di_log","Tave","Precipitation","Soil_N","Soil_ph","Wcont",
          "Funct_Di_log:Tave", "Funct_Di_log:Precipitation","Funct_Di_log:Soil_N","Funct_Di_log:Soil_ph","Funct_Di_log:Wcont",
          "Phylo_Di_log:Tave", "Phylo_Di_log:Precipitation","Phylo_Di_log:Soil_N","Phylo_Di_log:Soil_ph","Phylo_Di_log:Wcont")

MegaModelSummary_deal <- MegaModelSummary_deal %>%
  filter(Parameter %in% order) %>%
  arrange(match(Parameter, order))

MegaModelSummary_deal$Term_display = factor(MegaModelSummary_deal$Term_display, 
                                            levels = rev(MegaModelSummary_deal$Term_display))

# add group information
MegaModelSummary_deal$Group <- c("Field com", rep("Plant attributes", 2), rep("Climate", 1), 
                                 rep("Soil properties", 3), rep("Interaction", 3))
MegaModelSummary_deal$Group <- factor(MegaModelSummary_deal$Group, 
                                      levels = unique(MegaModelSummary_deal$Group))

ggplot(MegaModelSummary_deal, aes(x = Term_display, y = Std_Coefficient, fill = Group))+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0, size = 0.8, color = "black")+
  geom_point(size = 3.5, pch = 21)+
  geom_segment(aes(y = 0, yend = 0, x = 0.5, xend = 10.3), color = "black", linetype = "dashed") + 
  geom_text(aes(y = CI_high, label = paste("italic(p)==", `p.adj`)),
            parse = TRUE, hjust = -0.4, vjust = 0.4, size = 4) + 
  labs(x = '', 
       y = 'Standard regression coefficients', 
       title = "Best model: <i>R</i><sup>2</sup> = 0.262", tag = "(b)") +  
  theme_classic() + coord_flip() +  
  scale_fill_manual(values = c("#BC5546","#2F4590","#6EA3C5","#8E333A","#EFA961")) +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title =  element_text(color = "black", size = 14),
        legend.text = element_text(size = 9, color = "black"),
        plot.title = element_textbox(size = 12, color = "black", fill = "white",     
                                     box.color = "black", width = grid::unit(1, "npc"),padding = margin(5, 5, 5, 5),  
                                     margin = margin(b = 5), halign = 0.5,linetype = "solid"),
        plot.margin = margin(0.5,1.5,0.5,1.5, unit = "cm"),
        legend.position = 'none',
        plot.tag = element_text(size = 14, face = "bold")) +
  scale_x_discrete(expand = expansion(mult = c(0.00, 0.05))) + 
  scale_shape_manual(values = c(16,21)) -> Figure_S7b1; Figure_S7b1


library(viridis)
MegaModelSummary_deal$x_pos <- "Importance"
MegaModelSummary_deal$importance_lab <- sprintf("%.2f", MegaModelSummary_deal$sel_prop)
MegaModelSummary_deal$importance_lab <- round(MegaModelSummary_deal$sel_prop, 2)

ggplot(MegaModelSummary_deal, aes(x = x_pos, y = Term_display, fill = sel_prop)) +
  geom_tile(width = 1, height = 1, color = "white") +
  #scale_fill_viridis(option = "D", direction = -1) + # + 
  scale_fill_gradient2(low = "#0C1E60", mid = "white", high = "#A59590",
                       midpoint = 0, limit = c(-1, 1)) +
  geom_text(aes(label = (importance_lab)), color = "black", size = 3.5) + 
  theme_bw() + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(panel.grid = element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), 
        plot.margin = margin(l = 60, r = 20, t = 10, b = 40),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12, vjust = 1, hjust = 1, angle = 30, color = "black"),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") -> Figure_S7b2; Figure_S7b2

# 9.11 x 10.10
Figure_S7b1+Figure_S7b2 + plot_layout(widths = c(0.9,0.1)) -> Figure_S7b; Figure_S7b

################################## Figure S7c ##################################
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
       tag = "(c)") +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = c("#184C3F", "#E4CB8F", "#57320F")) +
  scale_color_manual(values = c("#184C3F", "#E4CB8F", "#57320F"), name = "Funct-Dist") +
  annotate("text", label = expression(italic(p) == 0.051), x = 0.2, y = -0.02, size = 4) + 
  mytheme + theme(legend.position = c(0.68,0.20)) -> Figure_S7c; Figure_S7c


################################## Figure S6d ##################################
pred_mode <- ggeffect(Final_model, terms = c("Tave","Phylo_Di_log"))
eff_mod_data <- data.frame(pred_mode)
colnames(eff_mod_data)[1] = "Tave"
colnames(eff_mod_data)[2] <- "Effect_size"
colnames(eff_mod_data)[6] = "Phylo_Di_log"

# back transform attributes_variable
eff_mod_data["Tave"] <- pd_attributes_variable$`scaled:center`["Tave"] + 
  pd_attributes_variable$`scaled:scale`["Tave"]*eff_mod_data["Tave"]

eff_mod_data["Effect_size"] <- pd_attributes_variable$`scaled:center`["Effect_size"] + 
  pd_attributes_variable$`scaled:scale`["Effect_size"]*eff_mod_data["Effect_size"]

eff_mod_data$Phylo_Di_log <- ifelse(eff_mod_data$Phylo_Di_log == "-1", "Low Phylo−Dist (- 1 SD)", 
                                    ifelse(eff_mod_data$Phylo_Di_log == "0", "Mean Phylo−Dist", "High Phylo−Dist (+ 1 SD)"))
eff_mod_data$Phylo_Di_log <- factor(eff_mod_data$Phylo_Di_log, levels = c("Low Phylo−Dist (- 1 SD)", "Mean Phylo−Dist", "High Phylo−Dist (+ 1 SD)"))


ggplot()+
  geom_line(data = eff_mod_data, mapping = aes(Tave, Effect_size, color = factor(Phylo_Di_log)), size = 1.25) +
  labs(x = expression("Annual average temperature (°C)"), 
       y = bquote(atop("Environmental effects", 
                       Ln ~ "(" ~ frac(Fungi-dist[" estimated in field"], 
                                       Fungi-dist[" estimated in greenhouse"]) ~ ")")),
       tag = "(d)") +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = c("#184C3F", "#E4CB8F", "#57320F")) +
  scale_color_manual(values = c("#184C3F", "#E4CB8F", "#57320F"), name = "Phylo-Dist") +
  annotate("text", label = expression(italic(p) == 0.001), x = 15, y = 0.02, size = 4) + 
  mytheme + theme(legend.position = c(0.68,0.20)) -> Figure_S7d; Figure_S7d

# 9.11 x 10.10
#ggsave("Figure S7d.pdf", plot = Figure_S7a, width = 18, height = 18, units = "in", dpi = 300)

#Figure_S7a/Figure_S7c -> aaaa
#ggsave("Figure S7aa.pdf", plot = aaaa, width = 5.5, height = 9.90, units = "in", dpi = 300)
