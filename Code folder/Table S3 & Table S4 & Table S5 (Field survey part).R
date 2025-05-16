################################################################################
############ Table S3 & Table S4 & Table S5 (Field survey part) ################
################################################################################

# Loading the R packages
library(openxlsx)
library(vegan)
library(Rmisc)
library(AICcPermanova)

# Rarefied to minimum sample size

# Soil sample grouping information
#Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
#Field_group$Sample_ID <- rownames(Field_group)
#Rmisc::summarySE(Field_group, measurevar = "Fungal_SR")

# Soil sample abundance information
#Field_otu_row2 = read.xlsx("Field_data_row_ASVs.xlsx", sheet = "row_otu", colNames = T, rowNames = T)
#Field_otu_row = Field_otu_row2[,rownames(Field_group)]
#Field_otu_row[1:6,1:6]
#dim(Field_otu_row)

## checking
#rownames(Field_group) %in% colnames(Field_otu_row)

## Tax INFORMATION
#tax_default <- Field_otu_row2[,(387:394)]

# Resampled by the minimum number of reads per sample
#library(microeco)
#library(mecodev)
#data_default <- microtable$new(sample_table = Field_group,otu_table = Field_otu_row, tax_table = tax_default)
#print(data_default)
#data_default$sample_table
#data_default$tidy_dataset() 
#print(data_default)
#data_default$sample_sums()%>% range

#set.seed(1234)
#data_default$rarefy_samples(sample.size = 9477)
#data_default$sample_sums()%>% range
#print(data_default)
#fungi_Flattening <- data_default$otu_table
#fungi_Flattening[1:5,1:5]

########################## Table S3 (Field survey) #############################

# notes: I have completed the above work, so I directly load the completed file
fungi_Flattening <- read.xlsx("Field_fungi_Flattening.xlsx", sheet = "field_flattening", rowNames = T, colNames = T)
#fungi_Flattening[1:6, 1:6]

Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID = rownames(Field_group)

# Data Transformation
Field_group$RS <- sqrt(Field_group$RS)
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Years <- as.factor(Field_group$Years)
Field_group$Site <- as.factor(Field_group$Site)
Field_group <- Field_group[colnames(fungi_Flattening), ] 
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
species_hel <- as.data.frame(decostand(t(fungi_Flattening), method = 'hellinger'))
species_hel[1:6,1:6]
dim(species_hel)
BC_dist_field <- vegdist(species_hel, method = 'bray')

#field_hel_no <- t(fungi_Flattening)/9477
#Bray_dist_field_no <- vegdist(field_hel_no, method = 'bray')

rownames(species_hel) %in% rownames(Field_group_scale)

# Table S3 (field survey part)
set.seed(1234)
mod1 = vegan::adonis2(BC_dist_field ~ Years + Site + Family/Species + 
                        Family:Years + Family:Site + Family:Years:Site + 
                        Family/Species:Years + Family/Species:Site,
                      data = Field_group_scale, permutations = 999)

Table_S3_PERMANOVA_field <- as.data.frame(mod1)[1:10,]
rownames(Table_S3_PERMANOVA_field)[1:9] <- c("Year", "Site", "Family", "Species", "Year × Family", "Site × Family", "Year × Site × Family", "Years × Species", "Site × Species")
Table_S3_PERMANOVA_field$R2_per <- paste0(round(Table_S3_PERMANOVA_field$R2*100,1),"%")
Table_S3_PERMANOVA_field$R2 <- round(Table_S3_PERMANOVA_field$R2,3)
Table_S3_PERMANOVA_field$`F` <- round(Table_S3_PERMANOVA_field$`F`,2)
Table_S3_PERMANOVA_field$p.adj <- p.adjust(Table_S3_PERMANOVA_field$`Pr(>F)`, method = "holm")
print(Table_S3_PERMANOVA_field[c(1:6,8,9,7,10), c(1,3,6,4,5,7)]) # reorder


########################## Table S4 (Field survey) #############################
set.seed(1234)
total_data <- Field_group_scale
# full model
mod_multifactor <- with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_N + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + 
                                               Fun_Di_log:Tave + Fun_Di_log:Precipitation + Fun_Di_log:Soil_N + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

AICcPermanova::AICc_permanova2(mod_multifactor)
mod_full_results <- as.data.frame(AICcPermanova::AICc_permanova2(mod_multifactor)); mod_full_results$form = "mod0 ~ full"

# Detect Funct-Dist × Precipitation
Simplified_mod1 <- with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_N + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + 
                                               Fun_Di_log:Tave + Fun_Di_log:Soil_N + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod1
AICcPermanova::AICc_permanova2(Simplified_mod1)
mod_Simplified1 <- as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod1)); mod_Simplified1$form = "mod1 ~ -Funct-Dist × Precipitation"


# Detect Funct-Dist × Soil N content
Simplified_mod2 <- with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_N + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + 
                                               Fun_Di_log:Tave + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod2
AICcPermanova::AICc_permanova2(Simplified_mod2)
mod_Simplified2 <- as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod2)); mod_Simplified2$form = "mod2 ~ -Funct-Dist × Soil N content"

# Detect Phylo-Dist × Soil water content 
Simplified_mod3 <- with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_N + Phy_Di_log:Soil_ph +
                                               Fun_Di_log:Tave + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod3
AICcPermanova::AICc_permanova2(Simplified_mod3)
mod_Simplified3 <- as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod3)); mod_Simplified3$form = "mod3 ~ -Phylo-Dist × Soil water content"


# Detect Funct-Dist × Annual average temperature
Simplified_mod4 <- with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_N + Phy_Di_log:Soil_ph +
                                               Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod4
AICcPermanova::AICc_permanova2(Simplified_mod4)
mod_Simplified4 <- as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod4)); mod_Simplified4$form = "mod4 ~ -Funct-Dist × Annual average temperature"


# Detect Phylo_Dist × Soil N content
Simplified_mod5 <- with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_ph +
                                               Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod5
AICcPermanova::AICc_permanova2(Simplified_mod5)
mod_Simplified5 <- as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod5)); mod_Simplified5$form = "mod5 ~ -Phylo_Dist × Soil N content"


# Detect Funct-Dist × Soil pH
Simplified_mod6 <- with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_ph +
                                               Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod6
AICcPermanova::AICc_permanova2(Simplified_mod6)
mod_Simplified6 <- as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod6)); mod_Simplified6$form = "mod6 ~ -Funct-Dsit × Soil pH"


# Detect Phylo-Dist × Soil pH
Simplified_mod7 <- with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation +
                                               Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod7
AICcPermanova::AICc_permanova2(Simplified_mod7)
mod_Simplified7 <- as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod7)); mod_Simplified7$form = "mod7 ~ -Phylo-Dist × Soil pH"

# Model comparison (based on based on Corrected Akaike Information Criterion (AICc))
mod_all_results <- rbind(mod_full_results,mod_Simplified1,mod_Simplified2,mod_Simplified3,mod_Simplified4,
                        mod_Simplified5,mod_Simplified6,mod_Simplified7)
min_AICc <- min(mod_all_results$AICc)
mod_all_results$DeltaAICc <- mod_all_results$AICc - min_AICc

### Model Selection-delta_aicc is set to 2
#mod_all_sub <- subset(mod_all_results, DeltaAICc < 2)
mod_all_sub <- mod_all_results
mod_all_sub$ind_Weight <- exp(-0.5 * mod_all_sub$DeltaAICc)
mod_all_sub$AICWeight <- mod_all_sub$ind_Weight/sum(mod_all_sub$ind_Weight)
rownames(mod_all_sub) <- mod_all_sub$form; mod_all_sub$form <- NULL
print(mod_all_sub) # mod7 was the best model

########################## Table S4 (Field survey) #############################
Table_S4_PERMANOVA_field <- as.data.frame(Simplified_mod7)
Table_S4_PERMANOVA_field$R2_per <- paste0(round(Table_S4_PERMANOVA_field$R2*100,1),"%")
Table_S4_PERMANOVA_field$R2 <- round(Table_S4_PERMANOVA_field$R2,3)
Table_S4_PERMANOVA_field$`F` <- round(Table_S4_PERMANOVA_field$`F`,2)
Table_S4_PERMANOVA_field$p.adj <- p.adjust(Table_S4_PERMANOVA_field$`Pr(>F)`, method = "holm")
rownames(Table_S4_PERMANOVA_field)[1:11] <- c("Site pool", "Phylo-Dist", "Funct-Dist", "Soil pH", "Wcont",
                                        "Soil N", "Temperature", "Precipitation", 
                                        "Phylo-Dist × Temperature", "Phylo-Dist × Precipitation", "Funct-Dist × Wcont")
print(Table_S4_PERMANOVA_field[c(1:12), c(1,3,6,4,5,7)]) # reorder


########################## Table S5 (Field survey) #############################
# Mantel test
# Single trait
Traits <- c("Chol","SLA","LDMC","SRL","FRR","RS")
single_trait_mantel_total <- NULL

for (iii in Traits) {
  ## 
  Traits_test <- Field_group[,c(iii, "Sample_ID")]
  Traits_test$Sample_ID <- NULL
  traits_dis <- vegdist(Traits_test, method = 'euclidean')  
  ## Bray_Curtis distance
  bray_dist <- BC_dist_field
  set.seed(1234)
  mantel_Bray_Curtis <- vegan::mantel(traits_dis, bray_dist, method = 'spearman', permutations = 999, na.rm = TRUE)
  mantel_Bray_Curtis_result <- data.frame(iii,mantel_Bray_Curtis$statistic, mantel_Bray_Curtis$signif)
  colnames(mantel_Bray_Curtis_result) <- c("Traits","Mantel_R","p_value")
  single_trait_mantel_total <- rbind(single_trait_mantel_total,mantel_Bray_Curtis_result)
}


# all traits data
Traits_test <- Field_group[,c("Chol","SLA","LDMC","SRL","FRR","RS","Sample_ID")]
Traits_test$Sample_ID <- NULL
traits_dis <- vegdist(Traits_test, method = 'euclidean')  

# Bray_Curtis
bray_dist <- BC_dist_field
set.seed(1234)
mantel_Bray_Curtis <- vegan::mantel(traits_dis, bray_dist, method = 'spearman', permutations = 999, na.rm = TRUE)
all_trait_mantel_total <- data.frame(mantel_Bray_Curtis$statistic, mantel_Bray_Curtis$signif)
colnames(all_trait_mantel_total) <- c("Mantel_R","p_value")

# plant phylogeny distance
# Bray_Curtis 
bray_dist <- BC_dist_field
##
phylo_dist <- cophenetic(tree)
phylo_dist_data <- reshape2::melt(phylo_dist, na.rm = T)
phylo_dist_data$names <- paste(phylo_dist_data$Var1, phylo_dist_data$Var2, sep = "_")
colnames(phylo_dist_data)[3] <- "Phylo_dist"

# Fungal distance matrix
fungal_dist_data <- reshape2::melt(as.matrix(bray_dist), na.rm = T)
colnames(fungal_dist_data)[1] <- "Sample_ID"
fungal_dist_data <- fungal_dist_data %>% left_join(Field_group[, c("Species","Sample_ID")], by = "Sample_ID")
colnames(fungal_dist_data)[1] <- "Var1"
colnames(fungal_dist_data)[2] <- "Sample_ID"
colnames(fungal_dist_data)[3] <- "Fungal_dist"
fungal_dist_data <- fungal_dist_data %>% left_join(Field_group[, c("Species","Sample_ID")], by = "Sample_ID")
fungal_dist_data$names <- paste(fungal_dist_data$Species.x, fungal_dist_data$Species.y, sep = "_")

#
Total_dist_data <- fungal_dist_data %>% left_join(phylo_dist_data[, c("names","Phylo_dist")], by = "names")
Fungal_dist <- reshape2::dcast(Total_dist_data, Var1  ~ Sample_ID , value.var = "Fungal_dist")
rownames(Fungal_dist) <- Fungal_dist$Var1
Fungal_dist <- Fungal_dist[,-1]
Phylo_dist <- reshape2::dcast(Total_dist_data, Var1  ~ Sample_ID , value.var = "Phylo_dist")
rownames(Phylo_dist) <- Phylo_dist$Var1
Phylo_dist <- Phylo_dist[,-1]

#colnames(Phylo_dist) %in% colnames(Fungal_dist)
#rownames(Phylo_dist) %in% rownames(Fungal_dist)

set.seed(1234)
mantel_Bray_Curtis <- vegan::mantel(as.dist(Phylo_dist), as.dist(Fungal_dist), method = 'spearman', permutations = 999, na.rm = TRUE)
all_phylo_mantel_total <- data.frame(mantel_Bray_Curtis$statistic, mantel_Bray_Curtis$signif)
colnames(all_phylo_mantel_total) <- c("Mantel_R","p_value")

# Merge datasets
all_trait_mantel_total$Traits <- "Trait dissimilarity"
all_phylo_mantel_total$Traits <- "Phylogenetic distance"
Table_S5_mantel_field <- rbind(single_trait_mantel_total, all_trait_mantel_total, all_phylo_mantel_total)
Table_S5_mantel_field$Predictors <- c("Leaf chlorophyll content", "Specific leaf area", "Leaf dry mass content",
                                       "Specific root length", "Fine-to-total root mass ratio", "Root-to-shoot mass ratio",
                                       "Trait dissimilarity", "Phylogenetic distance")
Table_S5_mantel_field$Mantel_R <- round(Table_S5_mantel_field$Mantel_R, 2)
# p.adjust
Table_S5_mantel_field$p.adj <- p.adjust(Table_S5_mantel_field$p_value, method = "holm")
print(Table_S5_mantel_field[,c(4,1:3,5)])

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.
