################################################################################
#################### Distinctiveness index calculation #########################
################################################################################

# We sequentially calculated distinctiveness indices for plant functional traits, phylogenetic relationships, and fungal community composition based on a common distinctiveness metric.

# Loading the R packages
library(openxlsx)
library(dplyr)
library(betapart)
library(phytools)
library(treeio)
library(funrar)

# Function to standardize (only if needed) ----
standr <- function(x){(x-min(x))/(max(x)-min(x))} 

################################# Field survey ################################# 
# Soil sample grouping information
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

# Soil samples-abundance table
Field_otu_raw <- read.xlsx("Field_data_raw_ASVs.xlsx", sheet = "raw_otu", colNames = T, rowNames = T)
Field_otu_raw <- Field_otu_row[ ,Field_group$Sample_ID]
Field_otu_raw[1:6, 1:6]
#colSums(Field_otu_raw)

# Simpson distance matrix
Field_otu_01 <- t(Field_otu_raw)
Field_otu_01[Field_otu_01 > 0] = 1 
fd <- beta.pair(Field_otu_01, index.family = "sorensen")
Sim_dist_field_mean <- fd$beta.sim

############################# Greenhouse experiment ############################
# Soil sample grouping information in greenhouse exp.
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "green_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)

# Soil samples-abundance table in greenhouse exp.
Green_otu_raw <- read.xlsx("Greenhouse_data_raw_ASVs.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
Green_otu_raw <- Green_otu_raw[ ,Green_group$Sample_ID]
Green_otu_raw[1:6,1:6]
#colSums(Green_otu_raw)

# Simpson distance matrix
Green_otu_01 <- t(Green_otu_raw)
Green_otu_01[Green_otu_01 > 0] <- 1 
fd <- beta.pair(Green_otu_01, index.family = "sorensen")
Sim_dist_green <- fd$beta.sim

# note:
# To match the rhizosphere fungal data with the corresponding species from the 
# field survey, we calculated the mean pairwise dissimilarity values among 
# species (averaged across three replicate samples).

Green_dist <- as.matrix(Sim_dist_green)
Green_dist_data <- reshape2::melt(Green_dist, varnames = c("Sample_ID_A", "Sample_ID_B"),
                                  value.name = "dist", na.rm = T)
colnames(Green_dist_data)[1] <- "Sample_ID"
Green_dist_data <- Green_dist_data %>% left_join(Green_group[,c("Sample_ID","Species")], by = "Sample_ID")
colnames(Green_dist_data)[c(1,2)] <- c("Sample_ID2","Sample_ID")
Green_dist_data <- Green_dist_data %>% left_join(Green_group[,c("Sample_ID","Species")], by = "Sample_ID")
aa = Rmisc::summarySE(Green_dist_data, measurevar = c("dist"), groupvars = c("Species.x", "Species.y"))
Sim_dist_green_mean <- reshape2::dcast(aa, Species.x  ~ Species.y , value.var = "dist")
rownames(Sim_dist_green_mean) <- Sim_dist_green_mean$Species.x
Sim_dist_green_mean <- Sim_dist_green_mean[,-1]
diag(Sim_dist_green_mean) = 0

############################# Functional traits ################################
# loading functional traits databases estimated in the greenhouse experiment
traits_mean = read.xlsx("traits_mean.xlsx", sheet = "traits_mean", colNames = T, rowNames = T)
traits_mean$Species = rownames(traits_mean)
colnames(traits_mean)
shapiro.test(sqrt(traits_mean$RS))
shapiro.test(log10(traits_mean$SRL))
traits_mean$RS = sqrt(traits_mean$RS)
traits_mean$SRL = log10(traits_mean$SRL)

########################## Plant phylogenetic tree #############################
tree <- read.newick("IQ_tree_plant_2025.NEWICK")
# Remove outergroup
to_drop <- c("Amborella_trichopoda","")
tree <- drop.tip(as.phylo(tree), to_drop) 

# Input Simpson distance matrix of fungi
Field_Fungal_dist <- as.matrix(Sim_dist_field_mean)
Green_Fungal_dist <- Sim_dist_green_mean

# Distinctiveness index calculation
Years_names <- as.vector(unique(Field_group$Years))
Latitude_names <- as.vector(unique(Field_group$Latitude))
Final_total_distinct <- NULL

for (Y in Years_names) {
  sample_group1 = subset(Field_group, Years == Y)
  for (L in Latitude_names) {
    sample_group2 = subset(sample_group1, Latitude == L)
    result <- sample_group2 %>% group_by(Species) %>% summarize(count = n()) %>% filter(count > 1)
    to_remove = c(result$Species)
    # If there are multiple samples of a species at a certain site
    if (length(to_remove) > 0) {
      for (AA in 1:length(to_remove)) {
        variable_name <- paste0("group", AA)
        # Extract data and assign it to variables
        data <- c((subset(sample_group2, Species == to_remove[AA]))$Sample_ID)
        assign(variable_name, data)
      }
      # According to the length of to_remove, the number of packets with corresponding length is created.
      selected_groups <- lapply(paste0("group", 1:length(to_remove)), get)
      # Use the expand.grid () function to generate all combinations.
      combinations <- as.data.frame(do.call(expand.grid, selected_groups))
      combinations[] <- lapply(combinations, as.character)
      dist.matr_field <- list()
      for (BB in 1:nrow(combinations)) {
        # Sample combination name after removing duplicate species data.
        rest_sample <- (sample_group2[!sample_group2$Species %in% to_remove, ])$Sample_ID
        sample_id <- unlist(c(as.vector(combinations[BB,]), rest_sample))
        # Evaluation of the distinctiveness of fungal composition based on simpson distance
        Fungal_mat_total <- as.matrix(Field_Fungal_dist) 
        standard.fungal_dist <- Fungal_mat_total[sample_id, sample_id]
        # distances are standardized to fit between 1 and 0 (if other distance metric is used)
        standard.fungal_dist <- standr(dist_mat) 
        dist.matr_field <- c(dist.matr_field,list(standard.fungal_dist))
        print(paste("Years:", Y, "Latitude:", L, "Repeat:", BB))
        dist_matrix_field <- as.data.frame(mean_matrix(dist.matr_field))
        # Mean of the computed distances for each species, to obtain the
        # average distance to the other species 
        standard.fungal_Di <- colSums(dist_matrix_field)/(nrow(dist_matrix_field)-1) #Mean of the computed distances for each species, to obtain the
        #average distance to the other species based on all trait combinations
        Fungal_field_Di <- as.data.frame(standard.fungal_Di)
        Fungal_field_Di$Sample_ID <- rownames(Fungal_field_Di); colnames(Fungal_field_Di)[1]<-c("Fungal_field_Di")
        Fungal_field_Di <- Fungal_field_Di %>% left_join(sample_group2[,c("Sample_ID","Species")], by = "Sample_ID")
        Fungal_field_Di$Sample_ID <- NULL
        Fungal_field_Di <- sample_group2[,c("Sample_ID","Species")] %>% left_join(Fungal_field_Di, by = "Species")
      }
    }
    else {
      # If there is only one sample of a species at a certain site
      sample_id <- sample_group2$Sample_ID
      Fungal_mat_total <- Field_Fungal_dist
      dist_matrix_field <- Fungal_mat_total[sample_id, sample_id]
      # distances are standardized to fit between 1 and 0 (if other distance metric is used)
      dist_matrix_field <- standr(dist_matrix_field)
      # Mean of the computed distances for each species, to obtain the
      # average distance to the other species 
      standard.fungal_Di <- colSums(dist_matrix_field)/(nrow(dist_matrix_field)-1) 
      Fungal_field_Di <- as.data.frame(standard.fungal_Di)
      Fungal_field_Di$Sample_ID <- rownames(Fungal_field_Di); colnames(Fungal_field_Di)[1] <- c("Fungal_field_Di")
      Fungal_field_Di <- Fungal_field_Di %>% left_join(sample_group2[,c("Sample_ID","Species")], by = "Sample_ID")
      Fungal_field_Di$Sample_ID = NULL
      Fungal_field_Di <- sample_group2[,c("Sample_ID","Species")] %>% left_join(Fungal_field_Di, by = "Species")
    }
    
    # Greenhouse experiment
    green_sample_id <- unique(sample_group2$Species)
    Fungal_Dist_green <- Green_Fungal_dist[green_sample_id, green_sample_id] 
    Fungal_Dist_green <- standr(Fungal_Dist_green)
    standard.fungal_Di <- colSums(Fungal_Dist_green)/(nrow(Fungal_Dist_green)-1)
    Fungal_green_Di <- as.data.frame(standard.fungal_Di)
    Fungal_green_Di$Species <- rownames(Fungal_green_Di); colnames(Fungal_green_Di)[1]<-c("Fungal_green_Di")
    
    # Functional Distinctiveness
    Fun_dist <- compute_dist_matrix((traits_mean[, c("Chol","SLA","LDMC","FRR","SRL","RS")]), metric = "euclidean", scale = TRUE, center = TRUE)
    Fun_dist_matrix <- Fun_dist[unique(sample_group2$Species), unique(sample_group2$Species)]
    Fun_dist_matrix <- standr(Fun_dist_matrix)
    Fun_Di <- as.data.frame(colSums(Fun_dist_matrix)/(nrow(Fun_dist_matrix)-1)) 
    Fun_Di$Species <- rownames(Fun_Di); colnames(Fun_Di)[1] <- c("Fun_Di")
    
    # Phylogenetic Distinctiveness
    Phylo_dist_mat <- as.data.frame(cophenetic(tree))
    Phylo_dist_mat <- Phylo_dist_mat[unique(sample_group2$Species), unique(sample_group2$Species)]
    Phylo_dist_mat <- standr(Phylo_dist_mat)
    Phy_Di <- as.data.frame(colSums(Phylo_dist_mat)/(nrow(Phylo_dist_mat)-1))
    Phy_Di$Species <- rownames(Phy_Di); colnames(Phy_Di)[1] <- c("Phy_Di")
    
    ### Merge data sets
    Final_Di <- left_join(Fungal_field_Di, Fungal_green_Di, by = "Species") %>%
      left_join(Phy_Di, by = "Species") %>%
      left_join(Fun_Di, by = "Species")
    Final_total_distinct = rbind(Final_total_distinct,Final_Di)
  }
}


head(Final_total_distinct)
#write.xlsx(Final_total_distinct, "Field_Distinctiveness.xlsx")
