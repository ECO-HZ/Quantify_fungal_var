################################################################################
################################## Figure S4 ###################################
################################################################################

# Loading R packages
library(openxlsx)
library(vegan)

############################## (Field survey part) #############################
# Soil sample grouping information
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

# Soil sample abundance information
Field_otu_raw <- read.xlsx("Field_data_raw_ASVs.xlsx", sheet = "raw_otu", colNames = T, rowNames = T)
Field_otu_raw <- Field_otu_raw[ ,Field_group$Sample_ID]

par(las = 1, mar = c(5, 5, 4, 2))
Field_min_depth <- min(rowSums(t(Field_otu_raw)))
rarecurve(as.data.frame(t(Field_otu_raw)), step = 100, 
          col = "black", cex = 0.6, label = FALSE, xlab = "", ylab = "", main = "")
abline(v = Field_min_depth, col = "red", lty = 2, lwd = 2)
title(main = "Field survey", xlab = "Number of sequences", ylab = "Number of ASVs")


###########################  (Greenhouse exp. part) ############################
Green_group <- read.xlsx("Greenhouse_data_group.xlsx", sheet = "green_group", colNames = T, rowNames = T)
Green_group$Sample_ID <- rownames(Green_group)

Green_otu_raw <- read.xlsx("Greenhouse_data_row_ASVs.xlsx", sheet = "raw_ASVs", colNames = T, rowNames = T)
Green_otu_raw <- Green_otu_raw[ ,Green_group$Sample_ID]
dim(Green_otu_raw)

par(las = 1, mar = c(5, 5, 4, 2))
Green_min_depth <- min(rowSums(t(Green_otu_raw)))
rarecurve(as.data.frame(t(Green_otu_raw)), step = 100, 
          col = "black", cex = 0.6, label = FALSE, xlab = "", ylab = "", main = "")
abline(v = Green_min_depth, col = "red", lty = 2, lwd = 2)
title(main = "Greenhouse experiment", xlab = "Number of sequences", ylab = "Number of ASVs")

