################################################################################
#################################### Figure S7 #################################
################################################################################

# Loading the R packages
library(openxlsx)
library(vegan)
library(ggplot2)
library(dplyr)
library(ggtext)
library(patchwork)
library(viridis)

# Custom style
mytheme <- theme_classic() + 
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none",
        legend.key = element_blank(),
        panel.grid=element_blank(), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
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

################################ (Field survey) ################################
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
dim(Field_otu_raw)
Field_otu_raw <- Field_otu_raw[rowSums(Field_otu_raw) > 0, ]
dim(Field_otu_raw)

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

# Composition of rhizosphere overall fungi (Bray-Curtis distance matrix) in field survey
fungi_relative <- decostand(Field_otu_raw, method = "total", MARGIN = 2)
colSums(fungi_relative)
BC_dist_RE_abun <- vegdist(t(fungi_relative), method = 'bray')

################################### Figure S7 ##################################
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
  #geom_point(df, mapping = aes(PCoA1_mean, PCoA2_mean, fill = Site, shape = Years), size = 4)+
  #geom_errorbar(data = df,mapping = aes(ymax = PCoA2_mean+PCoA2_se, ymin=PCoA2_mean-PCoA2_se),width=0,size=0.3,alpha = 1)+#
  #geom_errorbarh(data = df,mapping = aes(xmax=PCoA1_mean+PCoA1_se,xmin=PCoA1_mean-PCoA1_se),height=0,size=0.3,alpha = 1) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) +
  labs(x=paste("PCoA1 (", format(100 * eig[1] / sum(eig), digits=3), "%)", sep=""),
       y=paste("PCoA2 (", format(100 * eig[2] / sum(eig), digits=3), "%)", sep="")) +
  scale_shape_manual(values = c(24,21,25)) + 
  theme_bw() + mytheme + theme(legend.position = "right") + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) +
  theme(plot.tag = element_text(size = 14, face = "bold")) + 
  theme(legend.position = "right") -> Figure_S7; Figure_S7

