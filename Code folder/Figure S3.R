################################################################################
################################## Figure S3 ###################################
################################################################################
# Robustness evaluation of fungal composition distinctiveness index (based on greenhouse exp.)

# Loading the R packages
library(betapart)
library(openxlsx)
library(dplyr)
library(Rmisc)
library(reshape2)
library(emmeans)
library(multcomp)
library(ggplot2)

# Custom style
mytheme <- theme(panel.background = element_rect(fill = 'transparent', colour = 'black'),
                 legend.position = 'none',
                 panel.grid = element_blank(), 
                 legend.title = element_blank(),
                 legend.text = element_text(size = 8),
                 legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
                 axis.ticks = element_line(color = 'black'),
                 axis.line = element_line(colour = 'black'), 
                 axis.title.x = element_text(colour = 'black', size = 13),
                 axis.title.y = element_text(colour = 'black', size = 13),
                 axis.text = element_text(colour = 'black', size = 11),
                 plot.tag = element_text(size = 14, face = "bold")) 

# Function to standardize (only if needed) ----
standr = function(x){(x-min(x))/(max(x)-min(x))} 

# Obtain the Simpson pair-wise dissimilarity (Species Level)
# Add group information
Group_otu <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "green_group", colNames = T, rowNames = F)

# Loading the soil fungal abundance samples from greenhouse experiments
Green_otu <- read.xlsx("Greenhouse_data_raw_ASVs.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
Green_otu <- Green_otu[ ,Group_otu$Sample_ID]
colSums(Green_otu)
Green_otu_01 <- t(Green_otu)
Green_otu_01[Green_otu_01>0] = 1 
#View(Green_otu_01)
fd <- beta.pair(Green_otu_01, index.family = "sorensen")
Green_total_dis5 <- fd$beta.sim
Green_dist <- as.matrix(Green_total_dis5)

# Transform sample-level dissimilarity matrix into species-level dissimilarity matrix
Green_dist_data <- reshape2::melt(Green_dist, varnames = c("Sample_ID_A", "Sample_ID_B"), value.name = "dist", na.rm = T)
colnames(Green_dist_data)[1] <- "Sample_ID"
Green_dist_data <- Green_dist_data %>% left_join(Group_otu[,c("Sample_ID","Species")], by = "Sample_ID")
colnames(Green_dist_data)[c(1,2)] <- c("Sample_ID2","Sample_ID")
Green_dist_data <- Green_dist_data %>% left_join(Group_otu[,c("Sample_ID","Species")], by = "Sample_ID")
aa <- Rmisc::summarySE(Green_dist_data, measurevar = c("dist"), groupvars = c("Species.x", "Species.y"))
Fungal_dist <- reshape2::dcast(aa, Species.x  ~ Species.y , value.var = "dist")
rownames(Fungal_dist) <- Fungal_dist$Species.x
Fungal_dist <- Fungal_dist[,-1]
diag(Fungal_dist) <- 0

# standardization
Fungal_dist <- standr(Fungal_dist)

# Sloop
Species <- unique(Group_otu$Species)
M <- 2
final_sp2 <- NULL
final_sp4 <- NULL
final_sp8 <- NULL
final_sp14 <- NULL
final_sp22 <- NULL
final_sp32 <- NULL
final_sp44 <- NULL

set.seed(123456)
for (S in Species) {
  # Target species
  remain_sp <- unique((subset(Group_otu, Species != S))$Species)
  for (N in 1:500) {
    ### 2
    sample_sp2 <- sample(remain_sp, M, replace = FALSE)
    sample_sp_dist2 <- Fungal_dist[c(S,sample_sp2), c(S,sample_sp2)]
    sp_dist2_DI <- data.frame(colSums(sample_sp_dist2)/(nrow(sample_sp_dist2)-1));colnames(sp_dist2_DI) = "Fungal_Di"
    sp_dist2_DI$taxon <- rownames(sp_dist2_DI)
    sp_dist2_DI$num <- M
    sp_dist2_DI <- sp_dist2_DI[S,]
    final_sp2 <- rbind(final_sp2, sp_dist2_DI)
    ### 4
    sample_sp4 <- sample(remain_sp, M+2)
    sample_sp_dist4 <- Fungal_dist[c(S,sample_sp4), c(S,sample_sp4)]
    sp_dist4_DI <- data.frame(colSums(sample_sp_dist4)/(nrow(sample_sp_dist4)-1));colnames(sp_dist4_DI) = "Fungal_Di"
    sp_dist4_DI$taxon <- rownames(sp_dist4_DI)
    sp_dist4_DI$num <- M+2
    sp_dist4_DI <- sp_dist4_DI[S,]
    final_sp4 <- rbind(final_sp4, sp_dist4_DI)
    ### 8
    sample_sp8 <- sample(remain_sp, M+6)
    sample_sp_dist8 <- Fungal_dist[c(S,sample_sp8), c(S,sample_sp8)]
    sp_dist8_DI <- data.frame(colSums(sample_sp_dist8)/(nrow(sample_sp_dist8)-1));colnames(sp_dist8_DI) = "Fungal_Di"
    sp_dist8_DI$taxon <- rownames(sp_dist8_DI)
    sp_dist8_DI$num <- M+6
    sp_dist8_DI <- sp_dist8_DI[S,]
    final_sp8 <- rbind(final_sp8, sp_dist8_DI)
    ### 14
    sample_sp14 <- sample(remain_sp, M+12)
    sample_sp_dist14 <- Fungal_dist[c(S,sample_sp14), c(S,sample_sp14)]
    sp_dist14_DI <- data.frame(colSums(sample_sp_dist14)/(nrow(sample_sp_dist14)-1));colnames(sp_dist14_DI) = "Fungal_Di"
    sp_dist14_DI$taxon <- rownames(sp_dist14_DI)
    sp_dist14_DI$num <- M+12
    sp_dist14_DI <- sp_dist14_DI[S,]
    final_sp14 <- rbind(final_sp14, sp_dist14_DI)
    ### 22
    sample_sp22 <- sample(remain_sp, M+20)
    sample_sp_dist22 <- Fungal_dist[c(S,sample_sp22), c(S,sample_sp22)]
    sp_dist22_DI <- data.frame(colSums(sample_sp_dist22)/(nrow(sample_sp_dist22)-1));colnames(sp_dist22_DI) = "Fungal_Di"
    sp_dist22_DI$taxon <- rownames(sp_dist22_DI)
    sp_dist22_DI$num <- M+20
    sp_dist22_DI <- sp_dist22_DI[S,]
    final_sp22 <- rbind(final_sp22, sp_dist22_DI)
    ### 32
    sample_sp32 <- sample(remain_sp, M+30)
    sample_sp_dist32 <- Fungal_dist[c(S,sample_sp32), c(S,sample_sp32)]
    sp_dist32_DI <- data.frame(colSums(sample_sp_dist32)/(nrow(sample_sp_dist32)-1));colnames(sp_dist32_DI) = "Fungal_Di"
    sp_dist32_DI$taxon <- rownames(sp_dist32_DI)
    sp_dist32_DI$num <- M+30
    sp_dist32_DI <- sp_dist32_DI[S,]
    final_sp32 <- rbind(final_sp32, sp_dist32_DI)
    ### 44
    sample_sp44 <- sample(remain_sp, M+42)
    sample_sp_dist44 <- Fungal_dist[c(S,sample_sp44), c(S,sample_sp44)]
    sp_dist44_DI <- data.frame(colSums(sample_sp_dist44)/(nrow(sample_sp_dist44)-1));colnames(sp_dist44_DI) = "Fungal_Di"
    sp_dist44_DI$taxon <- rownames(sp_dist44_DI)
    sp_dist44_DI$num <- M+42
    sp_dist44_DI <- sp_dist44_DI[S,]
    final_sp44 <- rbind(final_sp44, sp_dist44_DI)
    print(paste("Species:", S, "Repeat:", N))
  }
}

total_di <- rbind(final_sp2,final_sp4,final_sp8,final_sp14,final_sp22,final_sp32,final_sp44)  
#head(total_di)
total_di$Species <- gsub("_", " ", total_di$taxon)
#unique(total_di$taxon)

# rename of Latin names
first_char <- substr(total_di$Species, 1, 1)
sub_str <- gsub(".*_", "", total_di$taxon)
Latin_name <- paste(first_char, sub_str, sep = ". ")
total_di$Latin_name <- Latin_name

#plot
ggplot() +
  #geom_jitter(total_di, mapping = aes(x = as.factor(num) , y = (Fungal_Di)),
  #position = position_jitter(0.3), size=0.3, color = "#00000022") +  #fill = "#00000022", stroke = 2
  labs(x = 'Add number of species',
       y = 'Fungal compositional distinctiveness')+
  geom_boxplot(data = total_di, mapping = aes(x = as.factor(num) , y = (Fungal_Di)), outlier.size = 0.8, outlier.shape = 21, size = 0.5) + 
  facet_wrap(~ Latin_name,ncol = 9, nrow = 6) + 
  mytheme + 
  theme(strip.text = element_text(size = 9, face = "italic"),
        #strip.background = element_rect(colour = "black"),
        axis.text = element_text(colour = 'black', size = 8)) ->  fig_S3

# one-way Anova
taxon <- unique(total_di$taxon)
final_data <- NULL

for (i in taxon) {
  select_data <- subset(total_di, taxon == i)
  mod <- aov(Fungal_Di ~ factor(num), data = select_data)
  results <- summary(mod)
  summary_data <- data.frame(Species = i, F_value = results[[1]][1,4], p_value = results[[1]][1,5])
  final_data <- rbind(final_data, summary_data)
}

print(subset(final_data, p_value <= 0.05))

# 
sp_select <- subset(final_data, p_value <= 0.05)$Species

for (sp_sel in sp_select) {
  total_di$num <- as.factor(total_di$num)
  select_data <- subset(total_di, taxon == sp_sel)
  mod <- aov(Fungal_Di ~ (num), data = select_data)
  print(sp_sel);print(summary(mod))
  LSSeasonNrate <- emmeans(mod, ~ num, adjust = "Tukey")
  Multiple_comp <- cld(LSSeasonNrate, alpha = 0.05, Letters = letters, adjust = "none", sort = F, rev = TRUE) 
  print(Multiple_comp)
}

