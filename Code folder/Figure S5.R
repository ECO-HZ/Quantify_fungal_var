################################################################################
################################## Figure S5 ###################################
################################################################################

# Loading the R packages
library(openxlsx)
library(corrplot)
library(Hmisc)

# Soil sample grouping information
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

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

cor.test(Field_group$Phylo_Di_log, Field_group$Precipitation)

Field_group_cor = Field_group[,c("Funct_Di_log","Phylo_Di_log",
                                 "Tave","Precipitation","Soil_N","Soil_ph","Wcont")]
colnames(Field_group_cor) = c("Funct-Dist","Phylo-Dist","Temperature","Precipitation","Soil N","Soil pH", "Wcont")

corr_matrix <- rcorr(as.matrix(Field_group_cor), type = 'spearman')
corr_matrix$r   
corr_matrix$P  

p.mat = corr_matrix$P
diag(p.mat) = 0

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(corr_matrix$r, p.mat = p.mat, sig.level = 0.05, insig = 'blank', method = 'number',type = 'lower',
         diag = F, col=col(200), tl.cex = 0.8,tl.col = "black", number.cex = 0.8, order = "original",tl.srt = 45)

corrplot(corr_matrix$r, p.mat = p.mat, sig.level = 0.05, insig = 'blank', method = 'square',
         add = TRUE, type = 'lower', diag = F, col=col(200), tl.pos = 'n', cl.pos = 'n',outline = F, order = "original",
         addCoef.col = "black", number.cex = 0.8)
