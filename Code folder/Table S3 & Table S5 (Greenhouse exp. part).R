################################################################################
############ Table S3 & Table S4 & Table S5 (Greenhouse exp. part) #############
################################################################################

# Loading R packages
library(openxlsx)
library(vegan)
library(dplyr)
library(phytools)
library(funrar)

# loading sample grouping information
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "green_group", colNames = T, rowNames = T)
Green_group$sample <- rownames(Green_group)
colnames(Green_group)

# loading sample grouping information
Green_otu_raw <- read.xlsx("Greenhouse_data_raw_ASVs.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
Green_otu_raw <- Green_otu_raw[ ,Green_group$sample]

# Read plant functional traits (mean values)
traits_mean <- read.xlsx("traits_mean.xlsx", sheet = "traits_mean", colNames = T, rowNames = T)
colnames(traits_mean)
# Transformation
shapiro.test(log10(traits_mean$RS))
shapiro.test(log10(traits_mean$SRL))
traits_mean$RS <- sqrt(traits_mean$RS)
traits_mean$SRL <- log10(traits_mean$SRL)

###################### Table S3 (Greenhouse exp. part) #########################
# Richness of overall fungi
Green_richness <- as.data.frame(specnumber(t(Green_otu_raw)))
Green_richness$sample <- rownames(Green_richness); colnames(Green_richness)[1] <- "Green_SR"

# Effect of family and species on fungal richness
Green_group <- Green_group %>% left_join(Green_richness)
mod <- lm(Green_SR ~ Family/Species , data = Green_group)
#mod <- aov(Green_SR ~ Family , data = Green_group)
#summary(mod)
Table_Fig_S3 <- as.data.frame(car::Anova(mod, type = 2))
Table_Fig_S3 <- Table_Fig_S3[-which(rownames(Table_Fig_S3) == "Residuals"),]
Table_Fig_S3$p.adj <- p.adjust(Table_Fig_S3$`Pr(>F)`, method = "holm")
Table_Fig_S3$p.adj <- sprintf("%.3f", Table_Fig_S3$p.adj)
Table_Fig_S3$`Pr(>F)` <- round(Table_Fig_S3$`Pr(>F)`, 3)
Table_Fig_S3$`F` <- round(Table_Fig_S3$`F value`, 2)
rownames(Table_Fig_S3) <- c("Family", "Species")
Table_Fig_S3$df <- c("15,108", "38,108")
print(Table_Fig_S3[ ,c("df","F","Pr(>F)","p.adj")])

###################### Table S4 (Greenhouse exp. part) #########################
# Permutational multivariate analysis of variance to explore the effects of plant families and species on fungal composition
# Composition of rhizosphere overall fungi (Bray-Curtis distance matrix) in Greenhouse experiment
Green_fungi_relative <- decostand(Green_otu_raw, method = "total", MARGIN = 2)
colSums(Green_fungi_relative)
bray_dist <- vegdist(t(Green_fungi_relative), method = 'bray')

set.seed(1234)
adonis_result <- vegan::adonis2(bray_dist ~ Family/Species, Green_group, method = "bray",permutations = 9999)
# p.adjust
Table_S4_PERMANOVA_Green <- as.data.frame(adonis_result)
rownames(Table_S4_PERMANOVA_Green)[1:2] <- c("Family", "Species")
Table_S4_PERMANOVA_Green$R2_per <- paste0(round(Table_S4_PERMANOVA_Green$R2*100,1),"%")
Table_S4_PERMANOVA_Green$R2 <- round(Table_S4_PERMANOVA_Green$R2,3)
Table_S4_PERMANOVA_Green$`F` <- round(Table_S4_PERMANOVA_Green$`F`,2)
#Table_S4_PERMANOVA_Green$p.adj <- p.adjust(Table_S4_PERMANOVA_Green$`Pr(>F)`, method = "holm")
print(Table_S4_PERMANOVA_Green[, c(1,3,6,4,5)]) # reorder

###################### Table S5 (Greenhouse exp. part) #########################
################################################################################
# Relationship between single traits and trait matrix, phylogenetic distance matrix 
# and community composition (Mantel test)
Green_group2 <- Green_group %>% left_join(traits_mean, by = c("Species", "Origin"))
rownames(Green_group2) <- Green_group2$Sample_ID

cor.test(Green_group2$Green_SR, Green_group2$LDMC)

# single traits
trait_names <- c("Chol","SLA","LDMC","SRL","FRR","RS")
single_trait_mantel <- NULL
for (i in trait_names) {
  all_otu_TEST <- Green_group2[,c(i, "Species")]
  all_otu_TEST$Species <- NULL
  traits_dis <- compute_dist_matrix(all_otu_TEST, metric = "euclidean", scale = TRUE, center = TRUE)
  ##Bray-Curtis
  set.seed(1234)
  mantel_Bray_Curtis <- vegan::mantel(as.dist(traits_dis), bray_dist, method = 'spearman', permutations = 999, na.rm = TRUE)
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
