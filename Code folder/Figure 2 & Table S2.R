################################################################################
############################## Figure 2 & Table S2 #############################
################################################################################

# Loading the R packages
library(openxlsx)
library(vegan)
library(ggplot2)
library(dplyr)
library(ggtext)
library(patchwork)

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


########################## Table S2 (Field survey) #############################
# Soil sample grouping information
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

# Soil sample abundance information
Field_otu_raw <- read.xlsx("Field_data_raw_ASVs.xlsx", sheet = "raw_otu", colNames = T, rowNames = T)
rownames(Field_otu_raw) <- Field_otu_raw$ASVs_ID
Field_otu_raw <- Field_otu_raw[ ,Field_group$Sample_ID]
Field_otu_raw[1:6, 1:6]
# colSums(Field_otu_raw)
# View(as.data.frame(rowSums(Field_otu_raw)))

# Data Transformation
Field_group <- Field_group[colnames(Field_otu_raw), ] # reorder
Field_group$Years <- as.factor(Field_group$Years)
Field_group$Site <- as.factor(Field_group$Site)
Field_group$RS <- sqrt(Field_group$RS)
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont*100)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Funct_Di_log <- log10(Field_group$Funct_Di)
Field_group$Phylo_Di_log <- log10(Field_group$Phylo_Di)

# Consider normalizing datasets
Field_group_scale <- Field_group

#colnames(Field_group_scale)
var_select <- c("Site_pool","Funct_Di","Phylo_Di","Funct_Di_log","Phylo_Di_log",
                "Soil_ph", "Wcont","Soil_N","Tave","Precipitation",
                "Chol","SLA","LDMC","SRL","FRR","RS")

Field_group_scale[var_select] = scale(Field_group_scale[c(var_select)])

# Composition of rhizosphere overall fungi (Bray-Curtis distance matrix) in field survey
fungi_relative <- decostand(Field_otu_raw, method = "total", MARGIN = 2)
colSums(fungi_relative)
BC_dist_RE_abun <- vegdist(t(fungi_relative), method = 'bray')

######################### Table S2 (Field survey part) #########################
set.seed(1234)
Table_S2_PERMANOVA <- vegan::adonis2(t(fungi_relative) ~ Family/Species + Site + Years +
                                       Family:Site + Family:Years + Family:Years:Site + 
                                       Family/Species:Site + Family/Species:Years + Years:Site, by = "term", method = "bray",
                                     data = Field_group_scale, permutations = 9999, parallel = 6)

Table_S2_PERMANOVA_field <- as.data.frame(Table_S2_PERMANOVA)[1:11,]
rownames(Table_S2_PERMANOVA_field)[1:10] <- c("Family", "Site", "Year", "Species", "Family × Site", 
                                              "Family × Year", "Site × Year", "Family × Site × Year", 
                                              "Species × Site", "Species × Year")
Table_S2_PERMANOVA_field$R2 <- round(Table_S2_PERMANOVA_field$R2,3)
Table_S2_PERMANOVA_field$`F` <- round(Table_S2_PERMANOVA_field$`F`,2)
print(Table_S2_PERMANOVA_field[c(1:7,9,10,8,11), ]) # reorder


################################## Figure 2a ###################################
# Principal Coordinates Analysis (PCoA)
pcoa <- cmdscale(BC_dist_RE_abun , k = 2, eig = TRUE)
plot_data <- data.frame({pcoa$points})[1:2]
plot_data$Sample_ID <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')
eig <- pcoa$eig
plot_data_add <- plot_data %>% 
  left_join(Field_group %>% select(-PCoA1, -PCoA2), by = "Sample_ID")

df = plot_data_add %>% dplyr::group_by(Site, Years) %>% 
  dplyr::summarise(PCoA1_mean = mean(PCoA1), PCoA1_se = sd(PCoA1)/(sqrt(length(PCoA1))),
                   PCoA2_mean = mean(PCoA2), PCoA2_se = sd(PCoA2)/(sqrt(length(PCoA2))))
#Rmisc::summarySE(plot_data_add, measurevar = c("PCoA1"), groupvars = c("Years", "Latitude", "Origin"))

# set colors of site
site_colors <- c("Guangzhou" = "#87898A", "Guilin" = "#C26275", "Changsha" = "#41479F",
                 "Wuhan" = "#32B7B2", "Zhengzhou" = "#75A750", "Tai'an" = "#E69F0D")

ggplot(df, aes(PCoA1_mean, PCoA2_mean))+
  geom_point(plot_data_add, mapping = aes(PCoA1, PCoA2, color = Site, fill = Site, shape = Years), size = 2) + 
  geom_point(df, mapping = aes(PCoA1_mean, PCoA2_mean, fill = Site, shape = Years), size = 4)+
  geom_errorbar(data = df,mapping = aes(ymax = PCoA2_mean+PCoA2_se, ymin=PCoA2_mean-PCoA2_se),width=0,size=0.3,alpha = 1)+#
  geom_errorbarh(data = df,mapping = aes(xmax=PCoA1_mean+PCoA1_se,xmin=PCoA1_mean-PCoA1_se),height=0,size=0.3,alpha = 1) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) +
  labs(x=paste("PCoA1 (", format(100 * eig[1] / sum(eig), digits=3), "%)", sep=""),
       y=paste("PCoA2 (", format(100 * eig[2] / sum(eig), digits=3), "%)", sep=""),
       tag = "(a)", title = "Field survey") +
  scale_shape_manual(values = c(24,21,25)) + 
  theme_bw() + mytheme + theme(legend.position = "right") + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) +
  theme(plot.tag = element_text(size = 14, face = "bold")) + 
  theme(legend.position = "right") -> Figure_2a; Figure_2a

################################## Figure 2b ###################################
permanova_Table_S2_field <- as.data.frame(Table_S2_PERMANOVA_field[c(1:7,9,10,8), ])
permanova_Table_S2_field$R2_per <- sprintf("%.1f", permanova_Table_S2_field$R2*100)
permanova_Table_S2_field$p_label <- ifelse(
  permanova_Table_S2_field$`Pr(>F)` < 0.001,
  "italic(p) < 0.001",
  paste0("italic(p) == ", formatC(permanova_Table_S2_field$`Pr(>F)`, format = "f", digits = 3, flag = "0")))
permanova_Table_S2_field$p_shape <- as.factor(ifelse(permanova_Table_S2_field$`Pr(>F)` >= 0.05, "p ≥ 0.05", "p < 0.05"))
permanova_Table_S2_field$Label <- rownames(permanova_Table_S2_field)
permanova_Table_S2_field$Label <- factor(permanova_Table_S2_field$Label, levels = rownames(permanova_Table_S2_field)) 
permanova_Table_S2_field <- as_tibble(permanova_Table_S2_field)

ggplot(permanova_Table_S2_field, aes(x = Label, y = as.numeric(R2_per))) + 
  geom_bar(stat="identity", aes(fill = p_shape), color = "black") +
  labs(x = NULL, y = "Explained variance (%)", tag = "(b)", fill = "Significance") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,42),
                     labels = scales::label_number(accuracy = 0.1)) +
  theme(panel.background = element_blank(),  
        plot.background = element_blank(),
        legend.position = c(0.8,0.8),
        plot.tag = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(colour='black',size=13),
        axis.text.y = element_text(colour='black',size=11),
        axis.text.x = element_text(colour='black',size=11, angle = 30, vjust = 1, hjust = 1),
        plot.margin = margin(l = 60, r = 20, t = 10, b = 40),
        axis.line.y = element_line(color = "black")) +
  scale_fill_manual(values = c("p ≥ 0.05" = "white", "p < 0.05" = "#3C4D6B")) -> Figure_2b; Figure_2b


###################### Table S2 (Greenhouse exp. part) #########################
# loading sample grouping information in greenhouse exp.
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "green_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)
colnames(Green_group)

# loading sample grouping information in greenhouse exp.
Green_otu_raw <- read.xlsx("Greenhouse_data_raw_ASVs.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
Green_otu_raw <- Green_otu_raw[ ,Green_group$Sample_ID]

Green_fungi_relative <- decostand(Green_otu_raw, method = "total", MARGIN = 2)
colSums(Green_fungi_relative)
Green_fungi_dist <- vegdist(t(Green_fungi_relative), method = 'bray')

set.seed(1234)
Table_S2_PER_green <- vegan::adonis2(t(Green_fungi_relative) ~ Family/Species, Green_group, method = "bray",permutations = 9999)
print(Table_S2_PER_green)

Table_S2_PERMANOVA_Green <- as.data.frame(Table_S2_PER_green[1:3, ])
rownames(Table_S2_PERMANOVA_Green)[1:2] <- c("Family", "Species")
Table_S2_PERMANOVA_Green$R2 <- round(Table_S2_PERMANOVA_Green$R2,3)
Table_S2_PERMANOVA_Green$`F` <- round(Table_S2_PERMANOVA_Green$`F`,2)
print(Table_S2_PERMANOVA_Green[, c(1,3,4,5)]) # reorder

################################### Figure 2c ##################################
# Principal Coordinates Analysis (PCoA)
pcoa <- cmdscale(Green_fungi_dist , k = 2, eig = TRUE)
plot_data <- data.frame({pcoa$points})[1:2]
plot_data$Sample_ID <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')
eig = pcoa$eig
plot_data <- merge(plot_data, Green_group, by = 'Sample_ID', all.x = TRUE)

plot_data$Origin = factor(plot_data$Origin, levels = c("Native","Exotic"))
plot_data$Family = as.factor(plot_data$Family)

###
df = plot_data %>% group_by(Family) %>% 
  summarise(PCoA1_mean = mean(PCoA1), PCoA1_se = sd(PCoA1)/(sqrt(length(PCoA1))),
            PCoA2_mean = mean(PCoA2), PCoA2_se = sd(PCoA2)/(sqrt(length(PCoA2))))

plot_data2 = merge(plot_data, df, by = c("Family"))
unique(df$Family)

family_color = c("Acanthaceae"="#003727","Amaranthaceae"="#005C4F","Asteraceae" ="#047266","Caryophyllaceae" = "#218B82", "Cyperaceae" = "#40A699",
                 "Euphorbiaceae" ="#6FBFB2", "Fabaceae" = "#97CFC9","Lamiaceae" ="#B4E2D8","Malvaceae"="#D2ECE8", "Onagraceae" = "#E3CD89",
                 "Phytolaccaceae" = "#D7B36A","Poaceae" ="#C79046","Polygonaceae"="#AF7423", "Solanaceae" = "#995C12", "Urticaceae" = "#7C4607", "Verbenaceae" = "#5F3401")

ggplot(df, aes(PCoA1_mean, PCoA2_mean))+
  geom_point(plot_data,mapping = aes(PCoA1, PCoA2, fill = Family, color = Family), size = 2, pch= 21) +
  geom_point(df,mapping = aes(PCoA1_mean, PCoA2_mean, color = Family, fill = Family), size = 4, pch = 21, color = "black")+
  geom_errorbar(data = df,mapping = aes(ymax = PCoA2_mean+PCoA2_se, ymin=PCoA2_mean-PCoA2_se),width=0,size=0.3,alpha = 1)+#
  geom_errorbarh(data = df,mapping = aes(xmax=PCoA1_mean+PCoA1_se,xmin=PCoA1_mean-PCoA1_se),height=0,size=0.3,alpha = 1) +
  scale_color_manual(values = family_color) +
  scale_fill_manual(values = family_color) +
  scale_shape_manual(values = c(21,16)) +
  labs(x=paste("PCoA1 (", format(100 * eig[1] / sum(eig), digits=3), "%)", sep=""),
       y=paste("PCoA2 (", format(100 * eig[2] / sum(eig), digits=3), "%)", sep=""),
       tag = "(c)", title = "Greenhouse experiment")+
  theme_bw() + mytheme + theme(legend.position = "right") + 
  theme(plot.tag = element_text(size = 14, face = "bold")) + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) -> Figure_2c; Figure_2c

################################### Figure 2d ##################################
permanova_Table_S2_green <- as.data.frame(Table_S2_PER_green[c(1:2), ])
rownames(permanova_Table_S2_green)[1:2] <- c("Family", "Species")
permanova_Table_S2_green$R2_per <- sprintf("%.1f", permanova_Table_S2_green$R2*100)
permanova_Table_S2_green$R2 <- round(permanova_Table_S2_green$R2,2)
permanova_Table_S2_green$p_label <- ifelse(
  permanova_Table_S2_green$`Pr(>F)` < 0.001,
  "italic(p) < 0.001",
  paste0("italic(p) == ", formatC(permanova_Table_S2_green$`Pr(>F)`, format = "f", digits = 3, flag = "0")))
permanova_Table_S2_green$p_shape <- as.factor(ifelse(permanova_Table_S2_green$`Pr(>F)` >= 0.05, 1, 0))
permanova_Table_S2_green$Label = rownames(permanova_Table_S2_green)
permanova_Table_S2_green$Label = factor(permanova_Table_S2_green$Label, levels = (rownames(permanova_Table_S2_green))) # rev.default
permanova_Table_S2_green <- as_tibble(permanova_Table_S2_green)

ggplot(permanova_Table_S2_green, aes(x = rev(Label), y = as.numeric(R2_per))) + 
  geom_bar(stat="identity", aes(fill = p_shape), color = "black") +
  labs(x = NULL, y = "Explained variance (%)", tag = "b") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,42)) +
  theme(panel.background = element_blank(),  
        plot.background = element_blank(),
        legend.position = "none",
        plot.tag = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(colour='black',size=13),
        axis.text.y = element_text(colour='black',size=11),
        axis.text.x = element_text(colour='black',size=11, angle = 30, vjust = 1, hjust = 1),
        plot.margin = margin(l = 60, r = 20, t = 10, b = 40),
        axis.line.y = element_line(color = "black")) +
  scale_fill_manual(values = c("1" = "white", "0" = "#3C4D6B")) -> Figure_2d; Figure_2d

