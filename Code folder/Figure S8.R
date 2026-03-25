################################################################################
################################# Figure S8 ####################################
################################################################################

# Loading the R packages
library(openxlsx)
library(car)
library(MuMIn)
library(ggplot2)
library(glmm.hp)
library(ggeffects)
library(patchwork)
library(ggtext)
library(dplyr)

# Custom style
mytheme = theme(panel.background = element_rect(fill='white', colour='black'),
                legend.position = "none",
                legend.key = element_blank(),
                #legend.background = element_blank(),   
                legend.box.background = element_blank(),
                panel.grid=element_blank(), 
                legend.title = element_text(size = 11),
                legend.text = element_text(size = 10),
                legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
                axis.ticks = element_line(color='black'),
                axis.line = element_line(colour = "black"), 
                axis.title.x = element_text(colour='black', size=13),
                axis.title.y = element_text(colour='black', size=13),
                axis.text = element_text(colour='black',size=11),
                plot.tag = element_text(size = 14, face = "bold")) 

# Loading field survey data
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)
Field_group$Years <- as.factor(Field_group$Years)
Field_group$Site <- factor(Field_group$Site, levels = unique(Field_group$Site[order(Field_group$Latitude)]))

# calculate environmental effects
Field_group$Effect_size <- log(Field_group$Fungal_Di_field_all/Field_group$Fungal_Di_green_all)

# Data Transformation
Field_group$RS <- sqrt(Field_group$RS)
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont*100)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Funct_Di_log <- log10(Field_group$Funct_Di)
Field_group$Phylo_Di_log <- log10(Field_group$Phylo_Di)

# Consider normalizing your data
var_select <- c("Effect_size","Site_pool","Phylo_Di","Funct_Di","Phylo_Di_log","Funct_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation","PCoA1")
pd_attributes_variable <- attributes(scale(Field_group[var_select]))
total_data <- Field_group
total_data[var_select] <- scale(total_data[var_select])

# performed the global linear regression model
fm1 <- lm(Effect_size ~ PCoA1 + Phylo_Di_log + Funct_Di_log + Soil_ph + Wcont + Soil_N + Precipitation + Tave +
            Phylo_Di_log:Tave + Phylo_Di_log:Precipitation + Phylo_Di_log:Soil_ph + Phylo_Di_log:Wcont + Phylo_Di_log:Soil_N + 
            Funct_Di_log:Tave + Funct_Di_log:Precipitation + Funct_Di_log:Soil_ph + Funct_Di_log:Wcont + Funct_Di_log:Soil_N, data = total_data)
summary(fm1)
car::vif(fm1)

options(na.action = "na.fail")
dd12 <- dredge(fm1, subset = ~ PCoA1 &
                 dc(Funct_Di_log, Tave, Funct_Di_log:Tave) & 
                 dc(Funct_Di_log, Precipitation, Funct_Di_log:Precipitation) & 
                 dc(Funct_Di_log, Soil_ph, Funct_Di_log:Soil_ph) & 
                 dc(Funct_Di_log, Wcont, Funct_Di_log:Wcont) & 
                 dc(Funct_Di_log, Soil_N, Funct_Di_log:Soil_N) &
                 ##
                 dc(Phylo_Di_log, Tave, Phylo_Di_log:Tave) & 
                 dc(Phylo_Di_log, Precipitation, Phylo_Di_log:Precipitation) & 
                 dc(Phylo_Di_log, Soil_ph, Phylo_Di_log:Soil_ph) & 
                 dc(Phylo_Di_log, Wcont, Phylo_Di_log:Wcont) & 
                 dc(Phylo_Di_log, Soil_N, Phylo_Di_log:Soil_N), trace = 2, rank = "AICc")

de6 <- model.avg(dd12, subset = delta < 2, fit = TRUE)
mod_list <- as.data.frame(de6$msTable)
mod_list$var_list <- rownames(mod_list); rownames(mod_list) <- NULL
mod_list$model_code <- 1:nrow(mod_list)
head(mod_list)
df_data = mod_list 

all_vars <- unique(unlist(strsplit(df_data$var_list, "\\+")))
all_vars <- sort(as.numeric(all_vars))

selection_matrix <- matrix(0, nrow = nrow(df_data), ncol = length(all_vars))
rownames(selection_matrix) <- paste0("Model_", df_data$model_code)
colnames(selection_matrix) <- all_vars

for(i in 1:nrow(df_data)) {
  selected_vars <- as.numeric(unlist(strsplit(df_data$var_list[i], "\\+")))
  selection_matrix[i, which(all_vars %in% selected_vars)] <- 1
}

plot_data <- melt(selection_matrix)
colnames(plot_data) <- c("Model", "Variable", "Selected")
head(plot_data)

plot_data <- plot_data %>%
  dplyr::left_join(df_data %>% 
                     dplyr::select(model_code, df, AICc, delta, weight) %>%
                     dplyr::mutate(Model = paste0("Model_", model_code)),
            by = "Model")

# adding the informations of code of predictors
print(de6$formula)

var_number_mapping <- data.frame(
  Variable = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18),
  Term = c(
    "Funct_Di_log",           
    "PCoA1",                  
    "Phylo_Di_log",           
    "Precipitation",         
    "Soil_N",                 
    "Soil_ph",               
    "Tave",                   
    "Wcont",                 
    "Funct_Di_log:Precipitation", 
    "Funct_Di_log:Soil_N",        
    "Funct_Di_log:Soil_ph",        
    "Funct_Di_log:Tave",           
    "Funct_Di_log:Wcont",          
    "Phylo_Di_log:Precipitation",  
    "Phylo_Di_log:Soil_N",         
    "Phylo_Di_log:Soil_ph",       
    "Phylo_Di_log:Tave",          
    "Phylo_Di_log:Wcont"))

print(var_number_mapping)

# rename predictors
var_number_mapping <- var_number_mapping %>%
  mutate(Term_display = case_when(
    Term == "PCoA1" ~ "Field fungal composition",
    Term == "Funct_Di_log" ~ "Funct-Dist",
    Term == "Funct_Di_log:Soil_ph" ~ "Funct-Dist × Soil pH",
    Term == "Funct_Di_log:Tave" ~ "Funct-Dist × Temperature",
    Term == "Funct_Di_log:Wcont" ~ "Funct-Dist × Wcont",
    Term == "Phylo_Di_log" ~ "Phylo-Dist",
    Term == "Phylo_Di_log:Precipitation" ~ "Phylo-Dist × Precipitation",
    Term == "Phylo_Di_log:Soil_N" ~ "Phylo-Dist × Soil N",
    Term == "Phylo_Di_log:Soil_ph" ~ "Phylo-Dist × Soil pH",
    Term == "Phylo_Di_log:Tave" ~ "Phylo-Dist × Temperature",
    Term == "Precipitation" ~ "Precipitation",
    Term == "Soil_N" ~ "Soil N",
    Term == "Soil_ph" ~ "Soil pH",
    Term == "Tave" ~ "Temperature",
    Term == "Wcont" ~ "Wcont",
    Term == "Funct_Di_log:Precipitation" ~ "Funct-Dist × Precipitation",
    Term == "Funct_Di_log:Soil_N" ~ "Funct-Dist × Soil N",
    Term == "Phylo_Di_log:Wcont" ~ "Phylo-Dist × Wcont",
    TRUE ~ Term))

plot_data <- plot_data %>% left_join(var_number_mapping)

# reorder
order = c("PCoA1","Funct_Di_log","Phylo_Di_log","Tave","Precipitation","Soil_N","Soil_ph","Wcont",
          "Funct_Di_log:Tave", "Funct_Di_log:Precipitation","Funct_Di_log:Soil_N","Funct_Di_log:Soil_ph","Funct_Di_log:Wcont",
          "Phylo_Di_log:Tave", "Phylo_Di_log:Precipitation","Phylo_Di_log:Soil_N","Phylo_Di_log:Soil_ph","Phylo_Di_log:Wcont")

plot_data <- plot_data %>% arrange(match(Term, order))
plot_data$Term_display <- factor(plot_data$Term_display, levels = unique(plot_data$Term_display))
head(plot_data)

ggplot(plot_data, aes(x = Term_display, y = reorder(Model, -AICc), fill = factor(Selected))) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "#A59590"),
                    labels = c("0" = "Not selected", "1" = "Selected")) +
  labs(title = NULL, x = NULL, y = "Model list (ΔAICc < 2)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 11, color='black'),
        axis.text.y = element_blank(),
        axis.title.y = element_text(colour='black', size=13),
        axis.ticks.x = element_line(color='black'),
        panel.grid = element_blank(),
        plot.margin = margin(t = 30, r = 40, b = 30, l = 80, unit = "pt"),
        legend.position = "none") -> Figure_S8a; Figure_S8a

# Model parameters
mod_list_filter = mod_list[,c("model_code", "AICc", "delta", "weight")]
head(mod_list_filter)

mod_list_long <- mod_list_filter %>%
  tidyr::pivot_longer(cols = c(AICc, delta, weight),
                      names_to = "Parameter",
                      values_to = "Value") %>%
  mutate(Value = round(Value, 2))
head(mod_list_long)

ggplot(mod_list_long, aes(x = factor(Parameter), y = factor(reorder(model_code, -model_code)))) +
  geom_tile(color = "white", size = 0.5, fill = "white") +
  geom_text(aes(label = Value), size = 3.5, color = "black") +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 11, color='black'),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(color='black'),
        panel.grid = element_blank(),
        plot.margin = margin(t = 30, r = 40, b = 30, l = 80, unit = "pt"),
        legend.position = "none") -> Figure_S8b; Figure_S8b

#
Figure_S8a|Figure_S8b + plot_layout(widths = c(0.9,0.1))
