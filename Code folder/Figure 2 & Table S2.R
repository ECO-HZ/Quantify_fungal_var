################################################################################
############################## Figure 2 & Table S2 #############################
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


########################## Table S2 (Field survey) #############################
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

# Consider normalizing datasets
Field_group_scale <- Field_group

#colnames(Field_group_scale)
var_select <- c("Site_pool","Funct_Di","Phylo_Di","Funct_Di_log","Phylo_Di_log",
                "Soil_ph", "Wcont","Soil_N","Tave","Precipitation",
                "Chol","SLA","LDMC","SRL","FRR","RS","Field_SR")

Field_group_scale[var_select] = scale(Field_group_scale[c(var_select)])

# Composition of rhizosphere overall fungi (Bray-Curtis distance matrix) in field survey
fungi_relative <- decostand(Field_otu_raw, method = "total", MARGIN = 2)
colSums(fungi_relative)
BC_dist_RE_abun <- vegdist(t(fungi_relative), method = 'bray')

# performed distance-based redundancy analysis (db-RDA)
Field_db_rda <- dbrda(BC_dist_RE_abun ~ Funct_Di_log + Phylo_Di_log + Tave + Precipitation + Soil_N + Soil_ph + Wcont +  
                  Funct_Di_log:Tave + Funct_Di_log:Precipitation + Funct_Di_log:Soil_N + Funct_Di_log:Soil_ph + Funct_Di_log:Wcont + 
                  Phylo_Di_log:Tave + Phylo_Di_log:Precipitation + Phylo_Di_log:Soil_N + Phylo_Di_log:Soil_ph + Phylo_Di_log:Wcont + 
                  Condition(Field_SR), 
                data = Field_group_scale, permutations = 999)
set.seed(123456)
db_rda_results <- anova(Field_db_rda, by = "term", permutations = 999)
print(db_rda_results)

# check results
db_rda_summary <- summary(Field_db_rda)

# Importance of components (proportion of total variance explained)
db_rda_summary$cont$importance 
# Importance of constrained components / Constrained axes importance
db_rda_summary$concont$importance 

site_scores <- scores(Field_db_rda, display = "sites") 
head(site_scores)
soil_scores <- scores(Field_db_rda, display = "bp")
head(soil_scores)

df_sites <- as.data.frame(site_scores)
df_soil <- as.data.frame(soil_scores)
multiplier <- 2.5
df_soil <- df_soil * multiplier

df_sites$Sample_ID = rownames(df_sites)
df_sites = df_sites %>% left_join(Field_group)
colnames(df_sites)

circled_numbers <- c("①", "②", "③", "④", "⑤", "⑥", "⑦", "⑧", "⑨", "⑩",
                     "⑪", "⑫", "⑬", "⑭", "⑮", "⑯", "⑰", "⑱", "⑲", "⑳")

df_soil$Label <- circled_numbers[1:nrow(df_soil)]

ggplot() +
  geom_point(data = df_sites, aes(x = dbRDA1, y = dbRDA2, color = Funct_Di_log), 
             alpha = 1, size = 3,  pch = 16) +                                
  geom_segment(data = df_soil,aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", alpha = 1, linewidth = 0.6) +
  #geom_text(data = df_soil,aes(x = dbRDA1 * 1.1, y = dbRDA2 * 1.1, label = rownames(df_soil)),
  #          color = "#B75E5F", size = 3, check_overlap = TRUE) +
  geom_text(data = df_soil,aes(x = dbRDA1 * 1.1, y = dbRDA2 * 1.1, label = Label),
            color = "black", size = 8, check_overlap = TRUE) +
  labs(x=paste0("dbRDA1 (", round(Field_db_rda$CCA$eig[1]/sum(Field_db_rda$CCA$eig)*100,2), "%)"), 
       y=paste0("dbRDA2 (", round(Field_db_rda$CCA$eig[2]/sum(Field_db_rda$CCA$eig)*100,2), "%)"),title = "") +
  #scale_fill_gradient(low = "white", high = "#7F4C2E",  # 白色到 viridis 的深紫色
  #                    na.value = "transparent") + 
  scale_color_viridis(option = "D", direction = -1) + # 
  mytheme + theme(legend.position = "right") + 
  labs(title = NULL, fill = "Funct-Dist (log10)") + 
  #scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  #scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) +
  theme(plot.tag = element_text(size = 14, face = "bold")) +
  guides(color = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,barwidth = 10, barheight = 1)) +  
  theme(legend.position = "top") -> Figure2_part1; Figure2_part1


ggplot() +
  geom_point(data = df_sites, aes(x = dbRDA1, y = dbRDA2, color = Phylo_Di_log), 
             alpha = 1, size = 3,  pch = 16) +                                
  geom_segment(data = df_soil,aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", alpha = 1, linewidth = 0.6) +
  geom_text(data = df_soil,aes(x = dbRDA1 * 1.1, y = dbRDA2 * 1.1, label = Label),
            color = "black", size = 8, check_overlap = TRUE) +
  labs(x=paste0("dbRDA1(", round(Field_db_rda$CCA$eig[1]/sum(Field_db_rda$CCA$eig)*100,2), "%)"), 
       y=paste0("dbRDA2(", round(Field_db_rda$CCA$eig[2]/sum(Field_db_rda$CCA$eig)*100,2), "%)"),title = "") +
  #scale_fill_gradient(low = "white", high = "#7F4C2E",  # 白色到 viridis 的深紫色
  #                    na.value = "transparent") + 
  scale_color_viridis(option = "E", direction = -1) + # 
  mytheme + theme(legend.position = "right") + 
  labs(title = NULL, fill = "Phylo-Dist (log10)") + 
  #scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  #scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) +
  theme(plot.tag = element_text(size = 14, face = "bold")) +
  guides(color = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,barwidth = 10, barheight = 1)) +  
  theme(legend.position = "top") -> Figure2_part2; Figure2_part2                   

Figure2_part1|Figure2_part2 -> Figure2_top

#ggsave("Figure2_top.pdf", plot = Figure2_top, width = 16, height = 16, units = "in", dpi = 300)

# Sensitivity analyses on distance-based redundancy analysis (db-RDA)
vars <- c("Funct_Di_log", "Phylo_Di_log", "Tave", "Precipitation", "Soil_N", "Soil_ph", "Wcont")
orders <- combinat::permn(vars) 

#results_all <- NULL
i = 1

for(i in 1:length(orders)) {
  
  main_effects <- paste(orders[[i]], collapse = " + ")
  formula_str <- paste0("BC_dist_RE_abun ~", main_effects, 
  "+ Funct_Di_log:Tave + Funct_Di_log:Precipitation + Funct_Di_log:Soil_N + Funct_Di_log:Soil_ph + Funct_Di_log:Wcont + Phylo_Di_log:Tave + Phylo_Di_log:Precipitation + Phylo_Di_log:Soil_N + Phylo_Di_log:Soil_ph + Phylo_Di_log:Wcont + Condition(Field_SR)")
  
  # performed db-RDA
  set.seed(123456)
  Sensitivity_db_rda <- dbrda(as.formula(formula_str), data = Field_group_scale)
  Sensitivity_db_rda_results <- anova(Sensitivity_db_rda, by = "term", permutations = 999, parallel = 6)
  #print(Sensitivity_db_rda_results)

  # add sloop information
  result_df <- as.data.frame(Sensitivity_db_rda_results)
  result_df$Term <- rownames(Sensitivity_db_rda_results)
  result_df$R2 <- result_df$SumOfSqs/sum(result_df$SumOfSqs)  
  result_df$Order_ID <- i
  
  # 合并到总结果
  results_all <- rbind(results_all, result_df)
  print(paste("orders: ", i))
}

# loading data 
results_all = data.table::fread("db-RDA_allorder_results.csv", header = T)
head(results_all)
results_all = results_all[results_all$Term != "Residual", ]

name_mapping <- c(
  "Funct_Di_log:Tave" = "Funct_Di_log:Tave",
  "Tave:Funct_Di_log" = "Funct_Di_log:Tave",
  
  "Funct_Di_log:Precipitation" = "Funct_Di_log:Precipitation",
  "Precipitation:Funct_Di_log" = "Funct_Di_log:Precipitation",
  
  "Funct_Di_log:Soil_N" = "Funct_Di_log:Soil_N",
  "Soil_N:Funct_Di_log" = "Funct_Di_log:Soil_N",
  
  "Funct_Di_log:Soil_ph" = "Funct_Di_log:Soil_ph",
  "Soil_ph:Funct_Di_log" = "Funct_Di_log:Soil_ph",
  
  "Funct_Di_log:Wcont" = "Funct_Di_log:Wcont",
  "Wcont:Funct_Di_log" = "Funct_Di_log:Wcont",
  
  "Phylo_Di_log:Tave" = "Phylo_Di_log:Tave",
  "Tave:Phylo_Di_log" = "Phylo_Di_log:Tave",
  
  "Phylo_Di_log:Precipitation" = "Phylo_Di_log:Precipitation",
  "Precipitation:Phylo_Di_log" = "Phylo_Di_log:Precipitation",
  
  "Phylo_Di_log:Soil_N" = "Phylo_Di_log:Soil_N",
  "Soil_N:Phylo_Di_log" = "Phylo_Di_log:Soil_N",
  
  "Phylo_Di_log:Soil_ph" = "Phylo_Di_log:Soil_ph",
  "Soil_ph:Phylo_Di_log" = "Phylo_Di_log:Soil_ph",
  
  "Phylo_Di_log:Wcont" = "Phylo_Di_log:Wcont",
  "Wcont:Phylo_Di_log" = "Phylo_Di_log:Wcont"
)

results_all_deal <- results_all %>%
  distinct(Order_ID, Term, .keep_all = TRUE) %>%
  mutate(Term_unified = ifelse(Term %in% names(name_mapping), name_mapping[Term], Term))

check_unification <- results_all_deal %>%
  filter(grepl(":", Term)) %>%
  group_by(Term_unified) %>%
  summarise(original_names = paste(unique(Term), collapse = " | "),
            n_orders = n_distinct(Order_ID),
            .groups = "drop")

print(check_unification)

significance_results <- results_all_deal %>%
  group_by(Term_unified) %>%
  summarise(
    total_orders = n_distinct(Order_ID),
    sig_count = sum(`Pr(>F)` < 0.05, na.rm = TRUE),
    nonsig_count = sum(`Pr(>F)` >= 0.05, na.rm = TRUE), 
    sig_probability = sig_count / total_orders * 100,
    nonsig_probability = nonsig_count / total_orders * 100, 
    mean_p = mean(`Pr(>F)`, na.rm = TRUE),
    sd_p = sd(`Pr(>F)`, na.rm = TRUE),
    min_p = min(`Pr(>F)`, na.rm = TRUE),
    max_p = max(`Pr(>F)`, na.rm = TRUE)
  ) %>%
  arrange(desc(sig_probability))

significance_results_deal <- significance_results %>%
  mutate(Term_display = case_when(
    Term_unified == "Funct_Di_log" ~ "Funct-Dist",
    Term_unified == "Funct_Di_log:Soil_ph" ~ "Funct-Dist × Soil pH",
    Term_unified == "Funct_Di_log:Tave" ~ "Funct-Dist × Temperature",
    Term_unified == "Funct_Di_log:Wcont" ~ "Funct-Dist × Wcont",
    Term_unified == "Phylo_Di_log" ~ "Phylo-Dist",
    Term_unified == "Phylo_Di_log:Precipitation" ~ "Phylo-Dist × Precipitation",
    Term_unified == "Phylo_Di_log:Soil_N" ~ "Phylo-Dist × Soil N",
    Term_unified == "Phylo_Di_log:Soil_ph" ~ "Phylo-Dist × Soil pH",
    Term_unified == "Phylo_Di_log:Tave" ~ "Phylo-Dist × Temperature",
    Term_unified == "Precipitation" ~ "Precipitation",
    Term_unified == "Soil_N" ~ "Soil N",
    Term_unified == "Soil_ph" ~ "Soil pH",
    Term_unified == "Tave" ~ "Temperature",
    Term_unified == "Wcont" ~ "Wcont",
    Term_unified == "Funct_Di_log:Precipitation" ~ "Funct-Dist × Precipitation",
    Term_unified == "Funct_Di_log:Soil_N" ~ "Funct-Dist × Soil N",
    Term_unified == "Phylo_Di_log:Wcont" ~ "Phylo-Dist × Wcont",
    Term_unified == "Residuals" ~ "Residuals",
    Term_unified == "Total" ~ "Total",
    TRUE ~ Term_unified 
  ))

#
Term_order = c("Funct-Dist","Phylo-Dist","Temperature","Precipitation","Soil N","Soil pH","Wcont",
               "Funct-Dist × Temperature", "Funct-Dist × Precipitation","Funct-Dist × Soil N","Funct-Dist × Soil pH","Funct-Dist × Wcont",
               "Phylo-Dist × Temperature", "Phylo-Dist × Precipitation","Phylo-Dist × Soil N","Phylo-Dist × Soil pH","Phylo-Dist × Wcont")

significance_results_deal <- significance_results_deal %>%
  mutate(Term_display = factor(Term_display, levels = Term_order)) %>%
  arrange(Term_display)

plot_data_percent <- significance_results_deal %>%
  filter(!Term_unified %in% c("Residuals", "Total")) %>%
  select(Term_display, sig_probability, nonsig_probability) %>%
  tidyr::pivot_longer(cols = c(sig_probability, nonsig_probability), 
                      names_to = "Significance", 
                      values_to = "Percentage") %>%
  mutate(
    Significance = case_when(
      Significance == "sig_probability" ~ "p < 0.05)",
      Significance == "nonsig_probability" ~ "p ≥ 0.05)"
    ),
    Significance = factor(Significance, levels = c("p < 0.05)", 
                                                   "p ≥ 0.05)"))
  )

ggplot(plot_data_percent, aes(x = Term_display, y = Percentage, fill = Significance)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = c("p < 0.05)" = "#214A7B", 
                               "p ≥ 0.05)" = "#949FBB")) +
  #geom_hline(yintercept = 50, linetype = "dashed", color = "grey50", alpha = 0.5) + 
  labs(x = NULL, y = "Proportion of\ntotal models (%)", fill = "Significance") +
  theme(panel.background = element_blank(),  
        plot.background = element_blank(),
        legend.position = "right",
        plot.tag = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(colour='black',size=13),
        axis.text.y = element_text(colour='black',size=11),
        axis.text.x = element_text(colour='black',size=11, angle = 30, vjust = 1, hjust = 1),
        plot.margin = margin(l = 60, r = 20, t = 10, b = 40),
        axis.line.y = element_line(color = "black")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = expansion(mult = c(0, 0))) -> Figure2_part3; Figure2_part3


r2_statistics <- results_all_deal %>%
  group_by(Term_unified) %>%
  summarise(
    n_orders = n_distinct(Order_ID),
    n_total = n(),
    mean_R2 = mean(R2, na.rm = TRUE),
    median_R2 = median(R2, na.rm = TRUE),
    sd_R2 = sd(R2, na.rm = TRUE),
    se_R2 = sd_R2 / sqrt(n()), 
    ci_lower = mean_R2 - 1.96 * se_R2,
    ci_upper = mean_R2 + 1.96 * se_R2,
    n_unique_R2 = n_distinct(R2),
    .groups = "drop")

print(r2_statistics)

r2_statistics$R2_per <- sprintf("%.1f", r2_statistics$mean_R2*100)

r2_statistics_deal <- r2_statistics %>%
  mutate(Term_display = case_when(
    Term_unified == "Funct_Di_log" ~ "Funct-Dist",
    Term_unified == "Funct_Di_log:Soil_ph" ~ "Funct-Dist × Soil pH",
    Term_unified == "Funct_Di_log:Tave" ~ "Funct-Dist × Temperature",
    Term_unified == "Funct_Di_log:Wcont" ~ "Funct-Dist × Wcont",
    Term_unified == "Phylo_Di_log" ~ "Phylo-Dist",
    Term_unified == "Phylo_Di_log:Precipitation" ~ "Phylo-Dist × Precipitation",
    Term_unified == "Phylo_Di_log:Soil_N" ~ "Phylo-Dist × Soil N",
    Term_unified == "Phylo_Di_log:Soil_ph" ~ "Phylo-Dist × Soil pH",
    Term_unified == "Phylo_Di_log:Tave" ~ "Phylo-Dist × Temperature",
    Term_unified == "Precipitation" ~ "Precipitation",
    Term_unified == "Soil_N" ~ "Soil N",
    Term_unified == "Soil_ph" ~ "Soil pH",
    Term_unified == "Tave" ~ "Temperature",
    Term_unified == "Wcont" ~ "Wcont",
    Term_unified == "Funct_Di_log:Precipitation" ~ "Funct-Dist × Precipitation",
    Term_unified == "Funct_Di_log:Soil_N" ~ "Funct-Dist × Soil N",
    Term_unified == "Phylo_Di_log:Wcont" ~ "Phylo-Dist × Wcont",
    Term_unified == "Residuals" ~ "Residuals",
    Term_unified == "Total" ~ "Total",
    TRUE ~ Term_unified  
  ))


r2_statistics_deal <- r2_statistics_deal %>%
  filter(!Term_unified %in% c("Residuals", "Total")) %>%
  mutate(Term_display = factor(Term_display, levels = Term_order)) %>%
  arrange(Term_display) %>% as.data.frame()

colnames(r2_statistics_deal)

r2_statistics_deal$y_name = "Explained\nvariance (%)"
str(r2_statistics_deal)
ggplot(r2_statistics_deal, aes(x = Term_display, y = y_name)) +
  geom_point(shape = 21, fill = "#3C4D6B", color = "#3C4D6B",alpha = 1, size = 9) +
  geom_text(aes(label = R2_per), 
            color = "white", fontface = "bold", size = 3) +
  labs(x = NULL, y = NULL) +
  theme(panel.background = element_blank(),  
        plot.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank()) -> Figure2_part4; Figure2_part4

(Figure2_part4/Figure2_part3) + plot_layout(heights = c(0.1,0.9))


