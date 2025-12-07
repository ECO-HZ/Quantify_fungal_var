################################################################################
################################## Figure 3 ####################################
################################################################################

# Loading the R packages
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggstar)
library(ggtext)
library(BestFitM)
#library(ggtrendline)
library(patchwork)

# Custom style
mytheme <- theme_bw() + 
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none",
        legend.key = element_blank(),
        panel.grid=element_blank(), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
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

# Loading field survey data
total_data <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
total_data$Sample_ID <- rownames(total_data)

# Data Transformation
total_data$Years <- as.factor(total_data$Years)
total_data$Site <- factor(total_data$Site, levels = unique(total_data$Site[order(total_data$Latitude)]))

total_data$RS <- sqrt(total_data$RS)
total_data$SRL <- log10(total_data$SRL)
total_data$Wcont <- sqrt(total_data$Wcont*100)
total_data$Soil_N <- sqrt(total_data$Soil_N)
total_data$Funct_Di_log <- log10(total_data$Funct_Di)
total_data$Phylo_Di_log <- log10(total_data$Phylo_Di)

cor.test(total_data$Funct_Di_log, total_data$Phylo_Di_log, method = "spearman")
cor.test(total_data$Funct_Di_log, total_data$Fungal_Di_field_all, method = "spearman")
cor.test(total_data$Phylo_Di_log, total_data$Fungal_Di_field_all, method = "spearman")

mod = lm(Fungal_Di_field_all ~ Fungal_Di_green_all, data = total_data)
anova(mod)
AIC(mod)
AICc(mod)
summary(mod)

# Mean ± 1SE
total_data_mean <- total_data %>% dplyr::group_by(Site, Years, Latitude) %>%
  dplyr::summarise(Field_Di = mean(Fungal_Di_field_all, na.rm = TRUE),
                   Field_Di_se = sd(Fungal_Di_field_all, na.rm = TRUE) / sqrt(n()),
                   Green_Di = mean(Fungal_Di_green_all, na.rm = TRUE),
                   Green_Di_se = sd(Fungal_Di_green_all, na.rm = TRUE) / sqrt(n()),
                   plant_richness = n(),
                   Funct_DI = mean(Funct_Di_log, na.rm = TRUE),
                   Phylo_DI = mean(Phylo_Di_log, na.rm = TRUE),
                   site_pool = mean(Site_pool)) %>%
  as.data.frame()

# set colors of site
site_colors <- c("Guangzhou" = "#87898A", "Guilin" = "#C26275", "Changsha" = "#41479F",
                 "Wuhan" = "#32B7B2", "Zhengzhou" = "#75A750", "Tai'an" = "#E69F0D")

ggplot()+
  geom_point(total_data, mapping = aes(Fungal_Di_green_all, Fungal_Di_field_all, color = Site, fill = Site, shape = Years), size = 2) + 
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) +
  labs(y = "Composition distinctiveness\nestimated in field survey", 
       x = "Composition distinctiveness\nestimated in greenhouse experiment", tag = "(a)") +
  scale_shape_manual(values = c(24,21,25)) + 
  geom_abline(intercept=0,slope=1 , linetype = 1, color = "#96383E", size = 0.8)+
  mytheme + 
  geom_smooth(data = total_data, aes(x =Fungal_Di_green_all , y = Fungal_Di_field_all), 
              method = "lm", formula = y ~ x, color = "black", se = F, linetype = 2, size = 0.8) + 
  theme_bw() + mytheme + 
  theme(legend.position = "right") + 
  annotate("segment", x = 0.56, xend = 0.59, y = 0.93, yend = 0.93, size = 1, color = "#96383E") + 
  annotate("text", x = 0.62, y = 0.93, label = "1:1 line", size = 4) +
  annotate("segment", x = 0.56, xend = 0.59, y = 0.91, yend = 0.91, size = 1, color = "black") + 
  annotate("text", x = 0.625, y = 0.91, label = "Model fit", size = 4) + 
  annotate("text", x = 0.68, y = 0.610, label = "F[1*','*383] < 0.01", parse = TRUE, size = 4) + 
  annotate("text", x = 0.68, y = 0.590, label = "italic(R)^2 < 0.01", parse = TRUE, size = 4) + 
  annotate("text", x = 0.68, y = 0.560, label = "italic(p) == 0.945", parse = TRUE, size = 4) +
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.56,0.95)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.56,0.95)) -> Figure_3a; Figure_3a


# calculate environmental effects
total_data$Effect_size <- log(total_data$Fungal_Di_field_all/total_data$Fungal_Di_green_all)

total_data2 = total_data %>% dplyr::group_by(Site, Years, Latitude) %>%
  mutate(plant_richness = n())
cor.test(total_data2$plant_richness, total_data2$Effect_size)

# Mean ± 1SE
Effect_size_mean <- total_data %>% dplyr::group_by(Site, Years, Latitude) %>%
  dplyr::summarise(Effect_size_mean = mean(Effect_size, na.rm = TRUE),
                   Effect_size_se = sd(Effect_size, na.rm = TRUE) / sqrt(n()),
                   Phylo_DI = mean(Phylo_Di_log, na.rm = TRUE),
                   Funct_DI = mean(Funct_Di_log, na.rm = TRUE),
                   plant_richness = n(),
                   site_pool = mean(Site_pool)) %>%
  as.data.frame()

cor.test(Effect_size_mean$Funct_DI, Effect_size_mean$Effect_size_mean, method = "spearman")
cor.test(Effect_size_mean$Phylo_DI, Effect_size_mean$Effect_size_mean, method = "spearman")
cor.test(Effect_size_mean$site_pool, Effect_size_mean$Effect_size_mean, method = "spearman")
cor.test(Effect_size_mean$plant_richness, Effect_size_mean$Effect_size_mean, method = "spearman")


# T.test
Years <- unique(total_data$Years)
Latitude <- unique(total_data$Latitude)
final_result_t <- NULL

for (Y in Years) {
  for (L in Latitude) {
    group_di <- subset(total_data, Years == Y & Latitude == L)
    X <-  group_di$Effect_size
    # 
    result_t <- t.test(X, mu = 0)
    result_t_test <- data.frame(Years = Y, Latitude = L, t_value = result_t$statistic,
                                p_value = result_t$p.value)
    final_result_t <- rbind(final_result_t, result_t_test)
    rownames(final_result_t) = NULL
  }
}

final_result_t$p_value <- round(final_result_t$p_value, 3)
print(subset(final_result_t, p_value >= 0.05))


Effect_size_mean = Effect_size_mean %>% left_join(final_result_t)
Effect_size_mean$sig <- ifelse(Effect_size_mean$p_value > 0.05, 0, 1)

Effect_size_mean$Years = as.factor(Effect_size_mean$Years)

# plot
ggplot() + 
  geom_errorbar(data = Effect_size_mean, 
                mapping = aes(x = Site, ymax = Effect_size_mean+Effect_size_se, 
                              ymin=Effect_size_mean-Effect_size_se, color = Years, group = Years),
                width=0.3,alpha = 1, color = "black", position = position_dodge(width = 0.3))+
  geom_point(data = subset(Effect_size_mean, Years == 2018 & Latitude != 23.1), 
             mapping = aes(x = Site, y = Effect_size_mean, shape = factor(Years), fill = Site),
             position = position_nudge(x = -0.1), size = 2.8, show.legend = FALSE) +
  geom_point(data = subset(Effect_size_mean, Years == 2018 & Latitude == 23.1), 
             mapping = aes(x = Site, y = Effect_size_mean, shape = factor(Years), fill = Site),
             position = position_nudge(x = -0.1), size = 2.8, show.legend = FALSE, fill = "white", color = "#87898A") +
  geom_point(data = subset(Effect_size_mean, Years == 2020 & Latitude != 23.1), 
             mapping = aes(x = Site, y = Effect_size_mean, shape = factor(Years), fill = Site),
             position = position_nudge(x = -0.0), size = 2.8, show.legend = FALSE) +
  geom_point(data = subset(Effect_size_mean, Years == 2020 & Latitude == 23.1), 
             mapping = aes(x = Site, y = Effect_size_mean, shape = factor(Years), fill = Site),
             position = position_nudge(x = -0.0), size = 2.8, show.legend = FALSE, fill = "white", color = "#87898A") +
  geom_point(data = subset(Effect_size_mean, Years == 2021 & Latitude != 27.9 & Latitude != 30.5), 
             mapping = aes(x = Site, y = Effect_size_mean, shape = factor(Years), fill = Site),
             position = position_nudge(x = 0.1), size = 2.8, show.legend = FALSE) +
  geom_point(data = subset(Effect_size_mean, Years == 2021 & Latitude == 27.9), 
             mapping = aes(x = Site, y = Effect_size_mean, shape = factor(Years), fill = Site),
             position = position_nudge(x = 0.1), size = 2.8, show.legend = FALSE, fill = "white", color = "#41479F") +
  geom_point(data = subset(Effect_size_mean, Years == 2021 & Latitude == 30.5), 
             mapping = aes(x = Site, y = Effect_size_mean, shape = factor(Years), fill = Site),
             position = position_nudge(x = 0.1), size = 2.8, show.legend = FALSE, fill = "white", color = "#32B7B2") +
  scale_shape_manual(values = c("2018" = 24, "2020" = 21, "2021" = 25)) +
  #scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) +
  geom_hline(yintercept = 0, linetype = 2) +
  mytheme + theme(legend.position = "right") +
  theme(legend.position = c(0.35,0.25),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + 
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) + 
  labs(#x = "Latitude (North degrees)", 
    x = NULL,
    y = bquote(atop("Environmental effects", 
                    Ln ~ "(" ~ frac(Fungi-dist["estimated in field survey"], 
                                    Fungi-dist["estimated in greenhouse experiment"]) ~ ")")),
    tag = "(b)") -> Figure_3b; Figure_3b
