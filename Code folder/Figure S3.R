################################################################################
################################# Figure S3 ####################################
################################################################################
# Loading the R packages
library(openxlsx)
library(ggvenn)
library(ggridges)

# Soil sample grouping information in field
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

# Soil sample abundance information in field
Field_otu_raw <- read.xlsx("Field_data_raw_ASVs.xlsx", sheet = "raw_otu", colNames = T, rowNames = T)
rownames(Field_otu_raw) <- Field_otu_raw$ASVs_ID
Field_otu_raw <- Field_otu_raw[ ,Field_group$Sample_ID]
Field_otu_raw[1:6, 1:6]
dim(Field_otu_raw)
Field_otu_raw = Field_otu_raw[rowSums(Field_otu_raw) > 0, ]
#View(as.data.frame(rowSums(Field_otu_raw)))

# loading sample grouping information in greenhouse exp.
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "green_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)
colnames(Green_group)

# loading sample grouping information in greenhouse exp.
Green_otu_raw <- read.xlsx("Greenhouse_data_raw_ASVs.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
rownames(Green_otu_raw) <- Green_otu_raw$ASV_ID
Green_otu_raw <- Green_otu_raw[ ,Green_group$Sample_ID]
dim(Green_otu_raw)
Green_otu_raw <- Green_otu_raw[rowSums(Green_otu_raw) > 0, ]
Green_otu_raw[1:6, 1:6]
sum(Green_otu_raw)

# show venn plot
field_ASV <- data.frame(ASV = rownames(Field_otu_raw), type = "Field")
green_ASV <- data.frame(ASV = rownames(Green_otu_raw), type = "Green")

################################# Figure S3a ###################################
list(`Field survey`=field_ASV$ASV, `Greenhouse experiment` = green_ASV$ASV) %>% 
  ggvenn(show_percentage = F,show_elements = F,label_sep = ",",
         digits = 1, stroke_color = "white",
         fill_color = c("#183C63", "#93C8C0"),
         set_name_color = c("#183C63", "#93C8C0")) -> Figure_S3a; Figure_S3a


################################# Figure S3b ###################################
# Filter common ASV IDs
common_ASVs <- intersect(rownames(Field_otu_raw), rownames(Green_otu_raw))
length(common_ASVs)

# converted the raw ASVs abundance data to relative abundances in greenhouse exp.
Green_otu_rel <- apply(Green_otu_raw, 2, function(x) x / sum(x))
Green_otu_rel <- as.data.frame(Green_otu_rel)
colSums(Green_otu_rel)

Green_long_abun <- data.frame(ASV = row.names(Green_otu_rel), Green_otu_rel) %>%
  tidyr::pivot_longer(cols = -ASV, names_to = "Sample_ID",values_to = "abun") %>%
  left_join(Green_group[,c("Sample_ID", "Species")])
head(Green_long_abun)

Green_mean_rel = Green_long_abun %>% group_by(Species, ASV) %>%
  summarise(mean_abun = mean(abun))
head(Green_mean_rel)

Green_mean_matrix <- Green_mean_rel %>%
  tidyr::pivot_wider(names_from = Species, values_from = mean_abun, values_fill = 0) %>%
  tibble::column_to_rownames("ASV") 
colSums(Green_mean_matrix)

# converted the raw ASVs abundance data to relative abundances in field survey
Field_otu_rel <- apply(Field_otu_raw, 2, function(x) x / sum(x))
Field_otu_rel <- as.data.frame(Field_otu_rel)
colSums(Field_otu_rel)

# The contribution of shared taxa to the variation in community composition was 
# calculated separately for each site and year.

Years <- unique(Field_group$Years)
Site <- unique(Field_group$Site)
diff_BC_merge_all <- NULL

for (i in Years) {
  for (ii in Site) {
    sub_group = subset(Field_group, Years == i & Site == ii)
    
    # Field survey
    Field_otu_rel_sub = Field_otu_rel[,sub_group$Sample_ID]
    
    # Step 1: Calculate compositional variation of the entire community.
    x <- apply(combn(ncol(Field_otu_rel_sub), 2), 2, function(x) sum(abs(Field_otu_rel_sub[,x[1]] - Field_otu_rel_sub[,x[2]]))/2 )
    x_names <- apply(combn(ncol(Field_otu_rel_sub), 2), 2, function(x) paste(colnames(Field_otu_rel_sub)[x], collapse=' - '))
    Field_BC_full <- data.frame(x, x_names)
    
    # Step 2: Calculate compositional variation of the common community.
    Field_common_sub <- Field_otu_rel_sub[common_ASVs, ]
    x <- apply(combn(ncol(Field_common_sub), 2), 2, function(x) sum(abs(Field_common_sub[, x[1]] - Field_common_sub[, x[2]]))/2)
    x_names <- apply(combn(ncol(Field_common_sub), 2), 2, function(x) paste(colnames(Field_common_sub)[x], collapse=' - '))
    Field_BC_common <- data.frame(x, x_names)
    
    # add site, year information
    diff_BC_Field <- left_join(Field_BC_full, Field_BC_common, by=c('x_names'))
    diff_BC_Field$diff_BC <- diff_BC_Field$x.y/diff_BC_Field$x.x
    diff_BC_Field$Years = i; diff_BC_Field$Site = ii; diff_BC_Field$Type = "Field"
    
    # add sample information
    split_names <- strsplit(diff_BC_Field$x_names, " - ")
    diff_BC_Field$Sample_ID1 <- sapply(split_names, function(x) x[1])
    diff_BC_Field$Sample_ID2 <- sapply(split_names, function(x) x[2])
    
    # add species information
    diff_BC_Field$Species_1 <- sub_group$Species[match(diff_BC_Field$Sample_ID1, rownames(sub_group))]
    diff_BC_Field$Species_2 <- sub_group$Species[match(diff_BC_Field$Sample_ID2, rownames(sub_group))]
    
    diff_BC_Field$Species_pair = paste0(diff_BC_Field$Species_1, " - ", diff_BC_Field$Species_2)
    colnames(diff_BC_Field)[2] = "Sample_pair"
    
    
    # Greenhouse experiment
    
    # Step 1: Calculate compositional variation of the entire community.
    Green_otu_rel_sub = Green_mean_matrix[, sub_group$Species]
    
    x <- apply(combn(ncol(Green_otu_rel_sub), 2), 2, function(x) sum(abs(Green_otu_rel_sub[,x[1]] - Green_otu_rel_sub[,x[2]]))/2 )
    x_names <- apply(combn(ncol(Green_otu_rel_sub), 2), 2, function(x) paste(colnames(Green_otu_rel_sub)[x], collapse=' - '))
    Green_BC_full <- data.frame(x, x_names)
    
    # Step 2: Calculate compositional variation of the common community.
    Green_common_matrix <- Green_otu_rel_sub[common_ASVs, ]
    x <- apply(combn(ncol(Green_common_matrix), 2), 2, function(x) sum(abs(Green_common_matrix[, x[1]] - Green_common_matrix[, x[2]]))/2)
    x_names <- apply(combn(ncol(Green_common_matrix), 2), 2, function(x) paste(colnames(Green_common_matrix)[x], collapse=' - '))
    Green_BC_common <- data.frame(x, x_names)
    
    # add site, year information
    diff_BC_Green <- left_join(Green_BC_full, Green_BC_common, by=c('x_names'))
    #head(diff_BC_Green)
    diff_BC_Green$diff_BC <- diff_BC_Green$x.y/diff_BC_Green$x.x
    #mean(diff_BC_Green$diff_BC)
    diff_BC_Green$Years = i; diff_BC_Green$Site = ii; diff_BC_Green$Type = "Greenhouse"
    
    # add sample information
    split_names <- strsplit(diff_BC_Green$x_names, " - ")
    diff_BC_Green$Species_1 <- sapply(split_names, function(x) x[1])
    diff_BC_Green$Species_2 <- sapply(split_names, function(x) x[2])
    
    # add species information
    diff_BC_Green$Sample_ID1 <- sub_group$Sample_ID[match(diff_BC_Green$Species_1, sub_group$Species)]
    diff_BC_Green$Sample_ID2 <- sub_group$Sample_ID[match(diff_BC_Green$Species_2, sub_group$Species)]
    
    diff_BC_Green$Sample_pair = paste0(diff_BC_Green$Sample_ID1, " - ", diff_BC_Green$Sample_ID2)
    colnames(diff_BC_Green)[2] = "Species_pair"
    
    # merge dataset
    # reorder
    diff_BC_Green = diff_BC_Green[,colnames(diff_BC_Field)]
    
    diff_BC_merge = rbind(diff_BC_Field, diff_BC_Green)
    
    # sloop 
    diff_BC_merge_all = rbind(diff_BC_merge_all, diff_BC_merge)
    
  }
}


# set color of site
site_colors <- c("Guangzhou" = "#87898A", "Guilin" = "#C26275", "Changsha" = "#41479F",
                 "Wuhan" = "#32B7B2", "Zhengzhou" = "#75A750", "Tai'an" = "#E69F0D")

diff_BC_merge_all$Site <- factor(diff_BC_merge_all$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))
diff_BC_merge_all$Type <- factor(diff_BC_merge_all$Type, levels = c("Field","Greenhouse"))
diff_BC_merge_all$Years <- as.factor(diff_BC_merge_all$Years)


diff_BC_merge_all %>% 
  group_by(Type) %>% 
  summarise(
    mean_diff_BC = mean(diff_BC*100),
    sd_diff_BC = sd(diff_BC*100),
    se_diff_BC = sd(diff_BC*100) / sqrt(n()))

# plot
ggplot(diff_BC_merge_all, aes(y = Years, x = diff_BC*100, fill = Type, color = Type)) +
  geom_density_ridges(rel_min_height = 0.01, scale = 0.9, alpha = 0.5, linewidth = 0.7,
                      quantile_lines = TRUE, quantile_fun = mean) +
  scale_color_manual(values = c("Field" = "#1C3C63", "Greenhouse" = "#93C8C0")) +
  scale_fill_manual(values = c("Field" = "#1C3C63", "Greenhouse" = "#93C8C0")) +
  scale_x_continuous(labels = scales::label_comma(accuracy = 1), limits = c(0,100), 
                     breaks = seq(0, 100, by = 25), expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(expand = expansion(mult = c(0.36, 0))) + 
  geom_vline(xintercept = 50, linetype = 2) +
  ggh4x::facet_grid2( ~ Site, #switch = "y", 
                      strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = site_colors))) +
  theme_minimal() +
  theme(legend.position = c(0.9,0.15), 
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        axis.line.y = element_line(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        axis.title.x = element_text(colour = 'black', size = 14),
        axis.title.y = element_text(colour = 'black', size = 14),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.text.x = element_text(colour = 'black', size = 12, angle = 90, vjust = -0.01),
        strip.background = element_rect(color=NA, size=0.5, linetype="solid"),
        strip.placement = "outside",
        strip.text.x = element_text(size = 12, colour = "black"),
        panel.spacing = unit(0, "lines")) +
  labs(x = 'Common taxa contributions\nto fungal community dissimilarity', y = NULL) +
  geom_segment(aes(x = 0, xend = 0, y = 1, yend = 3), color = "black") +
  coord_flip() -> Figure_S3b; Figure_S3b

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.
