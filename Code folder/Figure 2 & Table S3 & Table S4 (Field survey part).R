################################################################################
############ Figure 2 & Table S3 & Table S4 (Field survey part) ################
################################################################################

# Loading the R packages
library(openxlsx)
library(vegan)
library(Rmisc)
library(AICcPermanova)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(grid)  
library(ggtext)
library(ggstar) 
library(flextable)
library(officer)

# Custom style
mytheme <- theme_bw() + 
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none",
        legend.key = element_blank(),
        panel.grid=element_blank(), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
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

########################## Table S3 (Field survey) #############################
# Soil sample grouping information
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

# Soil sample abundance information
Field_otu_raw <- read.xlsx("Field_data_raw_ASVs.xlsx", sheet = "raw_otu", colNames = T, rowNames = T)
Field_otu_raw <- Field_otu_raw[ ,Field_group$Sample_ID]
Field_otu_raw[1:6, 1:6]
colSums(Field_otu_raw)
#View(as.data.frame(rowSums(Field_otu_raw)))

# Data Transformation
Field_group$RS <- sqrt(Field_group$RS)
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Years <- as.factor(Field_group$Years)
Field_group$Site <- as.factor(Field_group$Site)
Field_group <- Field_group[colnames(Field_otu_raw), ] 
#unique(Field_group$Species)

# Consider normalizing your data
Field_group_scale <- Field_group
shapiro.test(log10(Field_group_scale$Phy_Di))
shapiro.test(log10(Field_group_scale$Fun_Di))
Field_group_scale$Phy_Di_log <- log10(Field_group_scale$Phy_Di)
Field_group_scale$Fun_Di_log <- log10(Field_group_scale$Fun_Di)

#colnames(Field_group_scale)
Field_group_scale[c("Site_pool","Phy_Di","Fun_Di","Phy_Di_log","Fun_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation","Chol","SLA","LDMC","SRL","FRR","RS")] = 
  scale(Field_group_scale[c("Site_pool","Phy_Di","Fun_Di","Phy_Di_log","Fun_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation","Chol","SLA","LDMC","SRL","FRR","RS")])

# Composition of rhizosphere overall fungi (Bray-Curtis distance matrix) in field survey
fungi_relative <- decostand(Field_otu_raw, method = "total", MARGIN = 2)
colSums(fungi_relative)

#BC_dist_raw_abun <- vegdist(t(Field_otu_raw), method = 'bray')
BC_dist_RE_abun <- vegdist(t(fungi_relative), method = 'bray')

######################### Table S3 (field survey part) ######################### 
set.seed(1234)
Table_S3_PERMANOVA <- vegan::adonis2(t(fungi_relative) ~ Family/Species + Site + Years + 
                                       Family:Site + Family:Years + Family:Years:Site + 
                                       Family/Species:Site + Family/Species:Years + Years:Site, by = "term", method = "bray",
                                     data = Field_group_scale, permutations = 9999, parallel = 12)

Table_S3_PERMANOVA_field <- as.data.frame(Table_S3_PERMANOVA)[1:11,]
rownames(Table_S3_PERMANOVA_field)[1:10] <- c("Family", "Site", "Year", "Species", "Family × Site", "Family × Year", "Site × Year", "Family × Site × Year", "Species × Site", "Species × Year")
#Table_S3_PERMANOVA_field$R2_per <- paste0(round(Table_S3_PERMANOVA_field$R2*100,1),"%")
Table_S3_PERMANOVA_field$R2 <- round(Table_S3_PERMANOVA_field$R2,3)
Table_S3_PERMANOVA_field$`F` <- round(Table_S3_PERMANOVA_field$`F`,2)
print(Table_S3_PERMANOVA_field[c(1:7,9,10,8,11), ]) # reorder

########################## Table S4 (Field survey) #############################
total_data <- Field_group_scale
# full model
set.seed(1234)
mod_multifactor <- adonis2(t(fungi_relative) ~ Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation +
                             Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + Phy_Di_log:Soil_N + 
                             Fun_Di_log:Tave + Fun_Di_log:Precipitation + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont + Fun_Di_log:Soil_N, 
                           data = total_data, permutations = 9999, parallel = 10, by = "terms")
mod_multifactor
AICcPermanova::AICc_permanova2(mod_multifactor)
mod_full_results <- as.data.frame(AICcPermanova::AICc_permanova2(mod_multifactor)); mod_full_results$form = "mod0 ~ full"

# Detect Funct-Dist × Precipitation
set.seed(1234)
Simplified_mod1 <- adonis2(t(fungi_relative) ~ Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation +
                             Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + Phy_Di_log:Soil_N + 
                             Fun_Di_log:Tave + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont + Fun_Di_log:Soil_N, 
                           data = total_data, permutations = 9999, parallel = 10, by = "terms")

Simplified_mod1
AICcPermanova::AICc_permanova2(Simplified_mod1)
mod_Simplified1 <- as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod1)); mod_Simplified1$form = "mod1 ~ -Funct-Dist × Precipitation"


# Detect Funct-Dist × Soil N content
set.seed(1234)
Simplified_mod2 <- adonis2(t(fungi_relative) ~ Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation +
                             Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + Phy_Di_log:Soil_N + 
                             Fun_Di_log:Tave + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                           data = total_data, permutations = 9999, parallel = 10, by = "terms")

Simplified_mod2
AICcPermanova::AICc_permanova2(Simplified_mod2)
mod_Simplified2 <- as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod2)); mod_Simplified2$form = "mod2 ~ -Funct-Dist × Soil N content"

# Detect Phylo-Dist × Soil N content 
set.seed(1234)
Simplified_mod3 <- adonis2(t(fungi_relative) ~ Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation  +
                             Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + 
                             Fun_Di_log:Tave + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                           data = total_data, permutations = 9999, parallel = 10, by = "terms")

Simplified_mod3
AICcPermanova::AICc_permanova2(Simplified_mod3)
mod_Simplified3 <- as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod3)); mod_Simplified3$form = "mod3 ~ -Phylo-Dist × Soil N content"

# Detect Phylo-Dist × Wcont  
set.seed(1234)
Simplified_mod4 <- adonis2(t(fungi_relative) ~ Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation  +
                             Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_ph + 
                             Fun_Di_log:Tave + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                           data = total_data, permutations = 9999, parallel = 10, by = "terms")

Simplified_mod4
AICcPermanova::AICc_permanova2(Simplified_mod4)
mod_Simplified4 <- as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod4)); mod_Simplified4$form = "mod4 ~ -Phylo-Dist × Wcont"


# Model comparison (based on based on Corrected Akaike Information Criterion (AICc))
mod_all_results <- rbind(mod_full_results,mod_Simplified1,mod_Simplified2,mod_Simplified3,mod_Simplified4)
min_AICc <- min(mod_all_results$AICc)
mod_all_results$DeltaAICc <- mod_all_results$AICc - min_AICc

### Model Selection-delta_aicc is set to 2
#mod_all_sub <- subset(mod_all_results, DeltaAICc < 2)
mod_all_sub <- mod_all_results
mod_all_sub$ind_Weight <- exp(-0.5 * mod_all_sub$DeltaAICc)
mod_all_sub$AICWeight <- mod_all_sub$ind_Weight/sum(mod_all_sub$ind_Weight)
rownames(mod_all_sub) <- mod_all_sub$form; mod_all_sub$form <- NULL
print(mod_all_sub) # mod4 was the best model

########################## Table S4 (Field survey) #############################
Table_S4_PERMANOVA_field <- as.data.frame(Simplified_mod4)
Table_S4_PERMANOVA_field$R2_per <- paste0(round(Table_S4_PERMANOVA_field$R2*100,1),"%")
Table_S4_PERMANOVA_field$R2 <- round(Table_S4_PERMANOVA_field$R2,4)
Table_S4_PERMANOVA_field$`F` <- round(Table_S4_PERMANOVA_field$`F`,2)
rownames(Table_S4_PERMANOVA_field)[1:13] <- c("Phylo-Dist", "Funct-Dist", "Soil pH", "Wcont",
                                              "Soil N", "Temperature", "Precipitation", 
                                              "Phylo-Dist × Temperature", "Phylo-Dist × Precipitation", "Phylo-Dist × Soil pH",
                                              "Funct-Dist × Temperature", "Funct-Dist × Soil pH", "Funct-Dist × Wcont")
print(Table_S4_PERMANOVA_field[c(1:14), ]) # reorder

########################## Figure 2a (Field survey) ############################
Figure_2a <- Table_S4_PERMANOVA_field[c(1:13), c(1,3,6,4,5)]
Figure_2a$Predictors <- rownames(Figure_2a)
table_data <- Figure_2a[, c("Predictors", "F", "Pr(>F)", "R2")]
table_data$`Pr(>F)` <- round(table_data$`Pr(>F)`, 3)
table_data$`Pr(>F)` <- ifelse(table_data$`Pr(>F)` < 0.001, "<0.001", table_data$`Pr(>F)`)

# Create circle number sequence number
circled_numbers <- c("①", "②", "③", "④", "⑤", "⑥", "⑦", "⑧", "⑨", "⑩",
                     "⑪", "⑫", "⑬", "⑭", "⑮", "⑯", "⑰", "⑱", "⑲", "⑳")

table_data$Predictors <- paste0(circled_numbers[1:nrow(table_data)], table_data$Predictors)

colnames(table_data)[c(3:4)] <- c("p", "R2")

flextable(table_data) %>%
  theme_vanilla() %>%              
  set_table_properties(layout = "autofit") %>% 
  border_remove() %>%               
  align(j = 1, align = "left", part = "all") %>%   
  align(j = 2:4, align = "center", part = "all") %>%   
  italic(j = 3:4, part = "header") %>%     
  hline_top(part = "all", border = fp_border(width = 1.5)) %>% 
  hline_bottom(part = "all", border = fp_border(width = 1.5))
#set_header_labels(R2 = "R²") %>% 
hline(i = 1, part = "header", border = fp_border(width = 1.5)) 

########################## Figure 2b (Field survey) ############################
set.seed(1234)
field_nmds1 <- metaMDS(BC_dist_RE_abun, k = 3, trymax = 200, engine = "monoMDS", autotransform = TRUE, wascores = TRUE)
#field_nmds1 <- metaMDS(Bray_dist_field, k = 2, trymax = 100, autotransform = TRUE, wascores = TRUE)
field_nmds1.stress <- field_nmds1$stress
field_plot_data <- data.frame(field_nmds1$point)
field_plot_data$Sample_ID <- rownames(field_plot_data)
names(field_plot_data)[1:3] <- c('NMDS1', 'NMDS2', 'NMDS3')
field_plot_data <- merge(field_plot_data, Field_group, by = 'Sample_ID', all.x = TRUE)
field_plot_data$Years <- as.factor(field_plot_data$Years)
field_plot_data$Site <- factor(field_plot_data$Site, levels = c("Guangzhou", "Guilin", "Changsha", "Wuhan", "Zhengzhou", "Tai'an"))

# Vector fitting with NMDS
# Create interaction variable
Field_group_scale$Phy_Tave <- Field_group_scale$Phy_Di_log * Field_group_scale$Tave
Field_group_scale$Phy_Prec <- Field_group_scale$Phy_Di_log * Field_group_scale$Precipitation
Field_group_scale$Phy_Soil_ph <- Field_group_scale$Phy_Di_log * Field_group_scale$Soil_ph

Field_group_scale$Fun_Tave <- Field_group_scale$Fun_Di_log * Field_group_scale$Tave
Field_group_scale$Fun_Soil_ph <- Field_group_scale$Fun_Di_log * Field_group_scale$Soil_ph
Field_group_scale$Fun_Wcont <- Field_group_scale$Fun_Di_log * Field_group_scale$Wcont

# Run envfit analysis
NMDS.org.envt.fit <- envfit(data.frame(field_nmds1$point)[,1:2] ~  
                              Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                              Phy_Tave + Phy_Prec + Phy_Soil_ph + 
                              Fun_Tave + Fun_Soil_ph + Fun_Wcont,data = Field_group_scale, na.rm = TRUE)


# Step 1: extract envfit vector coordinates.
arrow_df <- as.data.frame(scores(NMDS.org.envt.fit, display = "vectors"))
arrow_df$Variable <- rownames(arrow_df)

#Step2: Create the coordinates of the arrow end point and scale it.

scale_factor <- 2.5
arrow_df <- arrow_df %>% mutate(x = 0, y = 0,
                                xend = MDS1 * scale_factor,
                                yend = MDS2 * scale_factor)

# Step3: Create a circle number+variable name label.
circled_numbers <- c("①", "②", "③", "④", "⑤", "⑥", "⑦", "⑧", "⑨", "⑩",
                     "⑪", "⑫", "⑬", "⑭", "⑮", "⑯", "⑰", "⑱", "⑲", "⑳")

arrow_df$Label <- circled_numbers[1:nrow(arrow_df)]

# Step 4: Draw NMDS Diagram+Arrow+Label
ggplot(field_plot_data, aes(NMDS1, NMDS2)) +
  geom_star(aes(starshape = Site, fill = Years),
            size = 2.2, color = "black", alpha = 1, show.legend = TRUE) +
  scale_starshape_manual(values = c(13,11,23,15,5,28)) +
  #scale_color_manual(values = c("2018" = "#0D6275", "2020" = "#752D28", "2021" = "#F1A627")) + 
  #scale_fill_manual(values = c("2018" = "#0D6275", "2020" = "#752D28", "2021" = "#F1A627")) +
  #scale_color_manual(values = c("2018" = "#40B0A6", "2020" = "#E1BE6A", "2021" = "#A38E89")) + 
  scale_fill_manual(values = c("2018" = "#40B0A6", "2020" = "#E1BE6A", "2021" = "#A38E89")) + 
  theme_bw() + mytheme +
  theme(legend.position = "right") +
  labs(x = "NMDS1", y = "NMDS2", tag = "(a)") +
  scale_x_continuous(labels = scales::label_comma(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy = 0.01)) +
  geom_segment(data = arrow_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "black", size = 0.6) +
  geom_text(data = arrow_df,
            aes(x = xend, y = yend, label = Variable),
            size = 4, hjust = 0.5, vjust = -0.5) +
  annotate("text", x = -Inf, y = Inf, label = "Stress = 0.19", 
           hjust = -0.1, vjust = 1.5,size = 4.2, color = "black") -> P2a

print(arrow_df)

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.
