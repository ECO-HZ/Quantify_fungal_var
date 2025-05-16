################################################################################
############ Table S3 & Table S4 & Table S5 (Greenhouse exp. part) #############
################################################################################

# Loading R packages
library(openxlsx)
library(vegan)
library(dplyr)
library(phytools)
library(funrar)

# Read in the fungal abundance information table
Green_otu <- read.xlsx("Greenhouse_fungi_Flattening.xlsx", sheet = "green_flattening", colNames = T, rowNames = T)
Green_otu[1:6,1:6]
Green_otu <- Green_otu[,-c(1)]
dim(Green_otu)

# Hellinger Transformation
species_hel <- as.data.frame(decostand(t(Green_otu), method = 'hellinger'))
species_hel[1:6,1:6]

# Read in sample grouping information
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "green_group", colNames = T, rowNames = T)
Green_group$sample <- rownames(Green_group)
colnames(Green_group)

# Read plant functional traits (mean values)
traits_mean <- read.xlsx("traits_mean.xlsx", sheet = "traits_mean", colNames = T, rowNames = T)
colnames(traits_mean)
# Transformation
shapiro.test(log10(traits_mean$RS))
shapiro.test(log10(traits_mean$SRL))
traits_mean$RS <- sqrt(traits_mean$RS)
traits_mean$SRL <- log10(traits_mean$SRL)

###################### Table S3 (Greenhouse exp. part) #########################
mod <- lm(Overall_Richness ~ Family/Species , data = Green_group)
Table_S3_lm <- as.data.frame(anova(mod))
Table_S3_lm$p.adj <- round(p.adjust(Table_S3_lm$`Pr(>F)`, method = "holm"), 3)
Table_S3_lm$`Pr(>F)` <- round(Table_S3_lm$`Pr(>F)`, 3)
rownames(Table_S3_lm)[1:2] <- c("Family", "Species")
print(Table_S3_lm)

###################### Table S4 (Greenhouse exp. part) #########################
# Permutational multivariate analysis of variance to explore the effects of plant families and species on fungal composition
set.seed(1234)
bray_dist <- vegdist(species_hel, method = 'bray')
adonis_result <- vegan::adonis2(bray_dist ~ Family/Species, Green_group, permutations = 999)
# p.adjust
Table_S4_PERMANOVA_Green <- as.data.frame(adonis_result)
rownames(Table_S4_PERMANOVA_Green)[1:2] <- c("Family", "Species")
Table_S4_PERMANOVA_Green$R2_per <- paste0(round(Table_S4_PERMANOVA_Green$R2*100,1),"%")
Table_S4_PERMANOVA_Green$R2 <- round(Table_S4_PERMANOVA_Green$R2,3)
Table_S4_PERMANOVA_Green$`F` <- round(Table_S4_PERMANOVA_Green$`F`,2)
Table_S4_PERMANOVA_Green$p.adj <- p.adjust(Table_S4_PERMANOVA_Green$`Pr(>F)`, method = "holm")
print(Table_S4_PERMANOVA_Green[, c(1,3,6,4,5,7)]) # reorder

###################### Table S5 (Greenhouse exp. part) #########################
################################################################################
# Relationship between single traits and trait matrix, phylogenetic distance matrix 
# and community composition (Mantel test)
trait_names <- c("Chol","SLA","LDMC","SRL","FRR","RS")

# single traits
single_trait_mantel <- NULL
Green_group2 <- Green_group %>% left_join(traits_mean, by = c("Species", "Origin"))
colnames(Green_group2)

cor.test(Green_group2$Overall_Richness, Green_group2$SRL)
cor.test(Green_group2$Overall_Richness, Green_group2$LDMC)

for (i in trait_names) {
  all_otu_TEST <- Green_group2[,c(i, "Species")]
  all_otu_TEST$Species <- NULL
  traits_dis <- vegdist(all_otu_TEST, method = 'euclidean')  
  ##Bray-Curtis
  set.seed(1234)
  mantel_Bray_Curtis <- vegan::mantel(traits_dis, bray_dist, method = 'spearman', permutations = 999, na.rm = TRUE)
  mantel_Bray_Curtis_result <- data.frame("Test",mantel_Bray_Curtis$statistic, mantel_Bray_Curtis$signif,i,'Bray_Curtis')
  colnames(mantel_Bray_Curtis_result)  <- c("Test","Mantel_R","P_value","Traits","dist_type")
  single_trait_mantel = rbind(single_trait_mantel,mantel_Bray_Curtis_result)
}
single_trait_mantel

# Trait matrix
Fun_dist <- compute_dist_matrix(traits_mean[,c("Chol","SLA","LDMC","FRR","SRL","RS")], metric = "euclidean", scale = TRUE, center = TRUE)
Fun_dist_data <- reshape2::melt(Fun_dist, na.rm = T)
Fun_dist_data$names <- paste(Fun_dist_data$Var1, Fun_dist_data$Var2, sep = "_")
colnames(Fun_dist_data)[3] <- "Funct_dist"
#head(Fun_dist_data)

# Phylogenetic distance matrix
tree <- read.newick("IQ_tree_plant_2025.NEWICK")
to_drop <- c("Amborella_trichopoda","")
tree <- drop.tip(as.phylo(tree), to_drop) 
plant_dist <- cophenetic(tree)
plant_dist_data <- reshape2::melt(plant_dist, na.rm = T)
plant_dist_data$names <- paste(plant_dist_data$Var1, plant_dist_data$Var2, sep = "_")
colnames(plant_dist_data)[3] <- "Phylo_dist"
#head(plant_dist_data)

Root_dist_data <- reshape2::melt(as.matrix(bray_dist), na.rm = T)
colnames(Root_dist_data)[1] <- "sample"
Root_dist_data = Root_dist_data %>% left_join(Green_group[, c("Species","sample")], by = "sample")
colnames(Root_dist_data)[1] <- "Var1"
colnames(Root_dist_data)[2] <- "sample"
colnames(Root_dist_data)[3] <- "Fungal_dist"
Root_dist_data <- Root_dist_data %>% left_join(Green_group[, c("Species","sample")], by = "sample")
Root_dist_data$names <- paste(Root_dist_data$Species.x, Root_dist_data$Species.y, sep = "_")
#head(Root_dist_data)

#
Total_dist_data <- Root_dist_data %>% left_join(plant_dist_data[, c("names","Phylo_dist")], by = "names") %>% 
  left_join(Fun_dist_data[, c("names","Funct_dist")], by = "names")

head(Total_dist_data)

Fungal_dist <- reshape2::dcast(Total_dist_data, Var1  ~ sample , value.var = "Fungal_dist")
rownames(Fungal_dist) <- Fungal_dist$Var1
Fungal_dist <- Fungal_dist[,-1]
##
Phylo_dist <- reshape2::dcast(Total_dist_data, Var1  ~ sample , value.var = "Phylo_dist")
rownames(Phylo_dist) <- Phylo_dist$Var1
Phylo_dist <- Phylo_dist[,-1]
##
Funct_dist <- reshape2::dcast(Total_dist_data, Var1  ~ sample , value.var = "Funct_dist")
rownames(Funct_dist) <- Funct_dist$Var1
Funct_dist = Funct_dist[,-1]

colnames(Phylo_dist) %in% colnames(Fungal_dist)
rownames(Phylo_dist) %in% rownames(Fungal_dist)

rownames(Funct_dist) %in% rownames(Fungal_dist)
rownames(Funct_dist) %in% rownames(Fungal_dist)

set.seed(1234)
mantel_fun <- vegan::mantel(as.dist(Funct_dist), as.dist(Fungal_dist), method = 'spearman', permutations = 999, na.rm = TRUE)
mantel_result_fun <- data.frame("Test",mantel_fun$statistic, mantel_fun$signif,"Trait dissimilarity",'Phylogenetic distance')
colnames(mantel_result_fun) <- c("Test","Mantel_R","P_value","Traits","dist_type")

set.seed(1234)
mantel_phy <- vegan::mantel(as.dist(Phylo_dist), as.dist(Fungal_dist), method = 'spearman', permutations = 999, na.rm = TRUE)
mantel_result_phy <- data.frame("Test",mantel_phy$statistic, mantel_phy$signif,"Trait dissimilarity",'Phylogenetic distance')
colnames(mantel_result_phy) <- c("Test","Mantel_R","P_value","Traits","dist_type")

# Merge databases
Table_S5_mantel_Green <- rbind(single_trait_mantel, mantel_result_fun, mantel_result_phy)
Table_S5_mantel_Green$Predictors <- c("Leaf chlorophyll content", "Specific leaf area", "Leaf dry mass content",
                                      "Specific root length", "Fine-to-total root mass ratio", "Root-to-shoot mass ratio",
                                      "Trait dissimilarity", "Phylogenetic distance")

Table_S5_mantel_Green$Mantel_R <- round(Table_S5_mantel_Green$Mantel_R, 2)
# p.adjust
Table_S5_mantel_Green$p.adj <- p.adjust(Table_S5_mantel_Green$P_value, method = "holm")
print(Table_S5_mantel_Green[,c(6,4,2,3,7)])


