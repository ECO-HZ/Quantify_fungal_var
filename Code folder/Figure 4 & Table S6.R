################################################################################
############################ Figure 4 & Table S6 ###############################
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
library(ggtext)
library(PieGlyph)
library(TITAN2)

# Custom style
mytheme = theme(panel.background = element_rect(fill='white', colour='black'),
                legend.position = "none",
                legend.key = element_blank(),
                #legend.background = element_blank(),   
                legend.box.background = element_blank(),
                panel.grid=element_blank(), 
                legend.title = element_text(size = 11),
                legend.text = element_text(size = 10),
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
Field_group$Years <- as.factor(Field_group$Years)
Field_group$Site <- factor(Field_group$Site, levels = unique(Field_group$Site[order(Field_group$Latitude)]))

# calculate environmental effects
Field_group$Effect_size <- log(Field_group$Fungal_field_Di/Field_group$Fungal_green_Di)

# Data Transformation
Field_group$RS <- sqrt(Field_group$RS)
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Phy_Di_log <- log10(Field_group$Phy_Di)
Field_group$Fun_Di_log <- log10(Field_group$Fun_Di)

# Consider normalizing your data
pd_attributes_variable <- attributes(scale(Field_group[c("Site_pool","Phy_Di","Fun_Di","Phy_Di_log","Fun_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation")]))
total_data <- Field_group
total_data[c("Site_pool","Phy_Di","Fun_Di","Phy_Di_log","Fun_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation")] = 
  scale(total_data[c("Site_pool","Phy_Di","Fun_Di","Phy_Di_log","Fun_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation")])


fm1 <- lm(Effect_size ~ Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Precipitation + Tave +
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


# Identified the best-fitting models for each response variable based on 
# Corrected Akaike Information Criterion (AICc)
Final_model <- get.models(dd12,1)[[1]] 
summary(Final_model)

################################### Table S6 ###################################
Table_S6 <- as.data.frame(car::Anova(Final_model, type = 2))
Table_S6 <- Table_S6[-which(rownames(Table_S6) == "Residuals"),]
Table_S6$p.adj <- p.adjust(Table_S6$`Pr(>F)`, method = "holm")
Table_S6$p.adj <- sprintf("%.3f", Table_S6$p.adj)
Table_S6$`Pr(>F)` <- round(Table_S6$`Pr(>F)`, 3)
Table_S6$`F` <- round(Table_S6$`F value`, 2)
Table_S6$Parameter <- rownames(Table_S6)
Table_S6$VIF <- round(car::vif(Final_model), 2)
Table_S6$df <- paste0("1", ",", "373")
rownames(Table_S6) <- c("Funct-Dist", "Phylo-Dist", "Soil N", "Soil pH", "Temperature", "Wcont", 
                        "Funct-Dist × Soil N", "Funct-Dist × Temperature", "Funct-Dist × Wcont",
                        "Phylo-Dist × Temperature", "Phylo-Dist × Wcont")
print(Table_S6[c(1,2,5,3,4,6,8,7,9,10,11) ,c("df","F","Pr(>F)","p.adj","VIF")])
performance::r2(Final_model)

################################## Figure 4a ###################################
# Obtaining standardized regression coefficients and their 95% CI
Final_model <- lm(Effect_size ~ Fun_Di_log + Phy_Di_log + Soil_N + 
                    Soil_ph + Tave + Wcont + Fun_Di_log:Soil_N + Fun_Di_log:Tave + 
                    Fun_Di_log:Wcont + Phy_Di_log:Tave + Phy_Di_log:Wcont, data = total_data)

MegaModelSummary <- as.data.frame(effectsize::effectsize(Final_model))[-1,]
MegaModelSummary$Parameter2 <- c("Funct-Dist", "Phylo-Dist", "Soil N", "Soil pH", "Temperature", "Wcont", 
                                 "Funct-Dist × Soil N", "Funct-Dist × Temperature", "Funct-Dist × Wcont",
                                 "Phylo-Dist × Temperature", "Phylo-Dist × Wcont")
MegaModelSummary$Parameter2 <- factor(MegaModelSummary$Parameter2, levels = rev(c("Funct-Dist", "Phylo-Dist", "Temperature", "Soil N", "Soil pH", "Wcont",
                                                                                  "Funct-Dist × Temperature", "Funct-Dist × Soil N", "Funct-Dist × Wcont", 
                                                                                  "Phylo-Dist × Temperature", "Phylo-Dist × Wcont")))
MegaModelSummary <- MegaModelSummary %>% left_join(Table_S6[,c("Parameter","p.adj")], by = "Parameter")

MegaModelSummary$Group <- c("Plant attributes", "Plant attributes", "Soil properties", "Soil properties", "Climate", "Soil properties", 
                            "Interaction", "Interaction", "Interaction", "Interaction", "Interaction")
MegaModelSummary$Group <- factor(MegaModelSummary$Group, levels = c("Plant attributes", "Climate", "Soil properties", "Interaction"))

ggplot(MegaModelSummary, aes(x = Parameter2, y = Std_Coefficient, fill = Group))+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0.2, size = 0.8, color = "black")+
  geom_point(size = 3.5, pch = 21)+
  geom_segment(aes(y = 0, yend = 0, x = 0.5, xend = 11.3), color = "black", linetype = "dashed") + 
  geom_text(aes(y = CI_high, label = paste("italic(p)==", `p.adj`)),
            parse = TRUE, hjust = -0.4, vjust = 0.4, size = 4) + 
  #annotate("text", x = 11.7, y = 0.08,
  #         label = expression("Best model: " * italic(R)^2 * "= 0.306"), colour = "black", size = 4) + 
  annotate("text", x = 9.0, y = 0.12,
           label = expression(italic(p) * " < 0.001"), colour = "black", size = 3.5) + 
  labs(x = '', 
       y = 'Standard regression coefficients', 
       title = "Best model: <i>R</i><sup>2</sup> = 0.306", tag = "(a)") +  
  theme_classic() + coord_flip() +  
  scale_fill_manual(values = c("#436C88","#49A3A4","#C2887D","#F0C986")) +
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
  scale_shape_manual(values = c(16,21)) -> p4a; p4a

# Relative contribution of variables
hierarchical_data <- as.data.frame(glmm.hp::glmm.hp(Final_model, type = "R2")$hierarchical.partitioning)
hierarchical_data$Group <- c("Plant attributes", "Plant attributes", "Soil properties", "Soil properties", "Climate", "Soil properties", 
                             "Interaction", "Interaction", "Interaction", "Interaction", "Interaction")
hierarchical_data$Parameter2 <- c("Funct-Dist", "Phylo-Dist", "Soil N", "Soil pH", "Temperature", "Wcont", 
                                  "Funct-Dist × Soil N", "Funct-Dist × Temperature", "Funct-Dist × Wcont",
                                  "Phylo-Dist × Temperature", "Phylo-Dist × Wcont")
hierarchical_data$Group <- factor(hierarchical_data$Group, levels = c("Plant attributes", "Climate", "Soil properties", "Interaction"))

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
  scale_fill_manual(values = c("#436C88","#49A3A4","#C2887D","#F0C986")) +
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
  labs(y = "Relative effect of estimates (%)") -> p4b; p4b

################################## Figure 4b ###################################
eff_mod <- effect("Phy_Di_log:Tave", Final_model, xlevels = 5)
#plot(eff_mod)
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
  geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = c("#184C3F","#6CB2AA", "#E4CB8F", "#BF8E43", "#57320F")) +
  scale_color_manual(values = c("#184C3F","#6CB2AA", "#E4CB8F", "#BF8E43", "#57320F"), name = expression("Phylo-Dist (log"[10]*")")) +
  #annotate("text", label = expression(italic(p) == 0.021), x = 15, y = 0.25, size = 4) + 
  mytheme + theme(legend.position = c(0.8,0.28)) -> P4b; P4b


################################## Figure 4c ###################################
# Soil sample abundance information
Field_otu_raw <- read.xlsx("Field_data_raw_ASVs.xlsx", sheet = "raw_otu", colNames = T, rowNames = T)
Field_otu_raw <- Field_otu_raw[ ,Field_group$Sample_ID]
Field_otu_raw[1:6, 1:6]
colSums(Field_otu_raw)
dim(Field_otu_raw)
#View(as.data.frame(rowSums(Field_otu_raw)))

# ASVs levels
Field_otu_raw_removed <- as.data.frame(t(Field_otu_raw))
Field_otu_raw_removed2 = Field_otu_raw_removed
Field_otu_raw_removed2[Field_otu_raw_removed2>0] <- 1
Field_otu_raw_removed <- Field_otu_raw_removed[, which(colSums(Field_otu_raw_removed2) >= 3)]
colSums(Field_otu_raw_removed)
dim(Field_otu_raw_removed)

set.seed(1234)
colnames(Field_group)
ENV_test <- Field_group[rownames(Field_otu_raw_removed),c("Sample_ID","Effect_size")]
ENV_test$Sample_ID = NULL 
titan_ENV_ASVs <- titan(env = ENV_test, txa = Field_otu_raw_removed, numPerm = 250, nBoot = 500, ncpus = 4)
titan_effect_size_ASVs = data.frame(titan_ENV_ASVs$sppmax) %>% filter(filter != 0) %>%
  arrange(maxgrp, zenv.cp)

# save results
#write.csv(titan_effect_size_ASVs, "titan_effect_size_ASVs.csv")
#load("titan_ENV_ASVs.RData")
#plot_taxa_ridges(titan_ENV_ASVs)

# loading responsive taxa information
responed_ASVs <- read.xlsx("Datasets S1.xlsx", sheet = "responed_ASVs", rowNames = F, colNames = T)
colnames(responed_ASVs)
unique(responed_ASVs$guild)

p1_sum_posit <- responed_ASVs %>% filter(filter == "Z+") %>% group_by(phylum, host_spec) %>%
  summarise(sum_count = sum(n()))

phylum_order <- p1_sum_posit %>%
  group_by(phylum) %>%
  summarise(total = sum(sum_count)) %>%
  arrange(total) %>%
  pull(phylum)

p1_sum_posit$phylum <- factor(p1_sum_posit$phylum, levels = phylum_order)


ggplot(data = p1_sum_posit,aes(x=phylum,y=sum_count,fill=host_spec))+
  geom_bar(stat = "identity",position = "stack", color = "black") + 
  labs(y = "Number of ASVs", x = "Phylum", title = "Positive (Z+) response taxa", tag = "(c)") + 
  theme_minimal() + 
  theme(panel.grid=element_blank(), 
        plot.title = element_textbox(size = 12, color = "black", fill = "white",     
                                     box.color = "black", width = grid::unit(1, "npc"),padding = margin(5, 5, 5, 5),  
                                     margin = margin(b = 5), halign = 0.5,linetype = "solid"),
        plot.margin = margin(0.5,1.5,0.5,1.5, unit = "cm"),
        axis.line.x = element_line(colour = "black"),
        axis.title.y = element_text(colour = "black", size = 14),
        axis.title.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.ticks = element_line(colour = "black"),
        plot.tag = element_text(size = 14, face = "bold"),
        legend.position = "none") + 
  scale_fill_manual(values = c("Core" = "#1C3C63", "Generalist" = "#D3D5D4", "Specialist" = "#93C8C0")) +    
  scale_y_continuous(limits = c(0,330), expand = c(0,0)) + 
  coord_flip() -> p1; p1


p2_sum_posit <- responed_ASVs %>% filter(filter == "Z+") %>% group_by(phylum, guild) %>%
  summarise(sum_count = sum(n()), prediction = 0.1)

p2_sum_posit$phylum <- factor(p2_sum_posit$phylum, levels = phylum_order)

ggplot(data = p2_sum_posit, aes(x = phylum, y = "Guild", fill = guild))+
  # Vertical line of lollipop
  #geom_segment(aes(yend = 0, xend = phylum))+
  # Pies-charts at the centre of the lollipop
  geom_pie_glyph(aes(pie_group = phylum), slices = 'guild', values = 'sum_count',
                 radius = 0.3, colour = 'black')+
  # Axis titles
  labs(y = 'Proportion (%)', x = 'Patient')+
  # Colours for sectors of the pie-chart
  scale_fill_manual(name = "Guild", 
                    values = c('Plant pathogen' = '#C74D26', 'Arbuscular mycorrhiza' = '#2C6344', 
                               'Others' = '#308194', 'Saprotroph' = '#61496D'))+
  theme_minimal() + 
  theme(panel.grid=element_blank(), 
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour = "black", size = 14),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.ticks = element_blank()) + 
  coord_flip() -> p2; p2



## 
p3_sum_negati <- responed_ASVs %>% filter(filter == "Z-") %>% group_by(phylum, host_spec) %>%
  summarise(sum_count = sum(n()))

negati_order <- p3_sum_negati %>%
  group_by(phylum) %>%
  summarise(total = sum(sum_count)) %>%
  arrange(total) %>%
  pull(phylum)

p3_sum_negati$phylum <- factor(p3_sum_negati$phylum, levels = negati_order)


ggplot(data = p3_sum_negati,aes(x=phylum,y=sum_count,fill=host_spec))+
  geom_bar(stat = "identity",position = "stack", color = "black") + 
  labs(y = "Number of ASVs", x = "Phylum", title = "Positive (Z-) response taxa") + 
  theme_minimal() + 
  theme(panel.grid=element_blank(), 
        plot.title = element_textbox(size = 12, color = "black", fill = "white",     
                                     box.color = "black", width = grid::unit(1, "npc"),padding = margin(5, 5, 5, 5),  
                                     margin = margin(b = 5), halign = 0.5,linetype = "solid"),
        plot.margin = margin(0.5,1.5,0.5,1.5, unit = "cm"),
        axis.line.x = element_line(colour = "black"),
        axis.title.y = element_text(colour = "black", size = 14),
        axis.title.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.ticks = element_line(colour = "black"),
        legend.position = c(0.85,0.35)) + 
  scale_fill_manual(values = c("Core" = "#1C3C63", "Generalist" = "#D3D5D4", "Specialist" = "#93C8C0")) +   
  scale_y_continuous(limits = c(0,600), expand = c(0,0)) + 
  coord_flip() -> p3; p3


p4_sum_negati <- responed_ASVs %>% filter(filter == "Z-") %>% group_by(phylum, guild) %>%
  summarise(sum_count = sum(n()), prediction = 0.1)

p4_sum_negati$phylum <- factor(p4_sum_negati$phylum, levels = negati_order)

ggplot(data = p4_sum_negati, aes(x = phylum, y = "Guild", fill = guild))+
  # Vertical line of lollipop
  #geom_segment(aes(yend = 0, xend = phylum))+
  # Pies-charts at the centre of the lollipop
  geom_pie_glyph(aes(pie_group = phylum), slices = 'guild', values = 'sum_count',
                 radius = 0.3, colour = 'black')+
  # Axis titles
  labs(y = 'Proportion (%)', x = 'Patient')+
  # Colours for sectors of the pie-chart
  scale_fill_manual(name = "Guild", 
                    values = c('Plant pathogen' = '#C74D26', 'Arbuscular mycorrhiza' = '#2C6344', 
                               'Others' = '#308194', 'Saprotroph' = '#61496D'))+
  theme_minimal() + 
  theme(panel.grid=element_blank(), 
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour = "black", size = 14),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.ticks = element_blank()) + 
  coord_flip() -> p4; p4


(p1|p2|p3|p4) + plot_layout(widths = c(0.4,0.1,0.4,0.1)) -> P4c; P4c 

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.
