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
total_data$Years <- as.factor(total_data$Years)
total_data$Site <- factor(total_data$Site, levels = unique(total_data$Site[order(total_data$Latitude)]))

# Data Transformation
total_data$RS <- sqrt(total_data$RS)
total_data$SRL <- log10(total_data$SRL)
total_data$Wcont <- sqrt(total_data$Wcont)
total_data$Soil_N <- sqrt(total_data$Soil_N)
total_data$Phy_Di_log <- log10(total_data$Phy_Di)
total_data$Fun_Di_log <- log10(total_data$Fun_Di)

# Pearson correlation analyses
cor.test(subset(total_data, Years == "2018")$Fungal_green_Di, subset(total_data, Years == "2018")$Fungal_field_Di)
summary(lm(subset(total_data, Years == "2018")$Fungal_field_Di ~ subset(total_data, Years == "2018")$Fungal_green_Di))
# r = 0.08, p = 0.403

cor.test(subset(total_data, Years == "2020")$Fungal_green_Di, subset(total_data, Years == "2020")$Fungal_field_Di)
summary(lm(subset(total_data, Years == "2020")$Fungal_field_Di ~ subset(total_data, Years == "2020")$Fungal_green_Di))
# r = -0.17, p = 0.054

cor.test(subset(total_data, Years == "2021")$Fungal_green_Di, subset(total_data, Years == "2021")$Fungal_field_Di)
summary(lm(subset(total_data, Years == "2021")$Fungal_field_Di ~ subset(total_data, Years == "2021")$Fungal_green_Di))
# r = -0.08, p = 0.368

# Mean ± 1SE
total_data_mean <- total_data %>% dplyr::group_by(Site, Years, Latitude) %>%
  dplyr::summarise(Field_Di = mean(Fungal_field_Di, na.rm = TRUE),
                   Field_Di_se = sd(Fungal_field_Di, na.rm = TRUE) / sqrt(n()),
                   Green_Di = mean(Fungal_green_Di, na.rm = TRUE),
                   Green_Di_se = sd(Fungal_green_Di, na.rm = TRUE) / sqrt(n()),
                   plant_richness = n()) %>%
  as.data.frame()


cor.test(total_data_mean$Field_Di, total_data_mean$plant_richness)
plot(total_data_mean$plant_richness, total_data_mean$Field_Di)
# Merge data
total_data2 <-  merge(total_data, total_data_mean, by = c("Years", "Site"))

################################## Figure 3a ###################################
ggplot() + 
  geom_star(data = subset(total_data, Years == "2018"), aes(x =Fungal_green_Di , y = Fungal_field_Di, color = Years, fill = Years, starshape = Site), 
            size = 1.2, alpha = 0.6, show.legend = F) + 
  scale_starshape_manual(values = c(13,11,23,15,5,28)) +
  geom_segment(subset(total_data2, Years == "2018"),mapping = aes(x = Green_Di, y = Field_Di, xend = Fungal_green_Di, yend = Fungal_field_Di),
               linetype = 1, linewidth = 0.2, alpha = 0.2, color = "#40B0A6")+
  geom_errorbar(data = subset(total_data_mean, Years == "2018"),mapping = aes(x = Green_Di,ymax = Field_Di+Field_Di_se, ymin=Field_Di-Field_Di_se),width=0.008,alpha = 1, color = "black")+
  geom_errorbarh(data = subset(total_data_mean, Years == "2018"),mapping = aes(y = Field_Di,xmax=Green_Di+Green_Di_se,xmin=Green_Di-Green_Di_se),height=0.008,alpha = 1, color = "black")+
  geom_star(data = subset(total_data_mean, Years == "2018"),mapping = aes(x = Green_Di, y = Field_Di, fill = Years, starshape = Site),size=2.5, color = "black")+
  geom_abline(intercept=0,slope=1 , linetype = 2, color = "black")+
  mytheme + 
  theme(plot.title = element_textbox(size = 14, color = "black", fill = "#40B0A6",
                                     box.color = "grey50",padding = margin(5, 5, 5, 5), margin = margin(b = 0),       
                                     halign = 0.5, width = grid::unit(1, "npc"))) + 
  ggpmisc::stat_poly_eq(data = subset(total_data, Years == "2018"), 
                        aes(x =Fungal_green_Di , y = Fungal_field_Di, 
                            label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.y.npc = "top") + 
  #annotate("text", x = 0.68, y = 0.60, label = bquote(atop(italic(r) == 0.08, italic(p) == 0.403)), size = 4) + 
  labs(y = "Composition distinctiveness\nestimated in field survey", x = NULL, tag = "(a)", title = "2018") + 
  scale_color_manual(values = c("2018" = "#40B0A6", "2020" = "#E1BE6A", "2021" = "#A38E89")) + 
  scale_fill_manual(values = c("2018" = "#40B0A6", "2020" = "#E1BE6A", "2021" = "#A38E89")) + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.56,0.94)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.56,0.94)) -> p2a;p2a


ggplot() + 
  geom_star(data = subset(total_data, Years == "2020"), aes(x =Fungal_green_Di , y = Fungal_field_Di, color = Years, fill = Years, starshape = Site), 
            size = 1.2, alpha = 0.6, show.legend = F) + 
  scale_starshape_manual(values = c(13,11,23,15,5,28)) +
  geom_segment(subset(total_data2, Years == "2020"),mapping = aes(x = Green_Di, y = Field_Di, xend = Fungal_green_Di, yend = Fungal_field_Di),
               linetype = 1, linewidth = 0.2, alpha = 0.2, color = "#E1BE6A")+
  geom_errorbar(data = subset(total_data_mean, Years == "2020"),mapping = aes(x = Green_Di,ymax = Field_Di+Field_Di_se, ymin=Field_Di-Field_Di_se),width=0.008,alpha = 1, color = "black")+
  geom_errorbarh(data = subset(total_data_mean, Years == "2020"),mapping = aes(y = Field_Di,xmax=Green_Di+Green_Di_se,xmin=Green_Di-Green_Di_se),height=0.008,alpha = 1, color = "black")+
  geom_star(data = subset(total_data_mean, Years == "2020"),mapping = aes(x = Green_Di, y = Field_Di, fill = Years, starshape = Site),size=2.5, color = "black")+
  geom_abline(intercept=0,slope=1 , linetype = 2, color = "black")+
  mytheme + 
  ggpmisc::stat_poly_eq(data = subset(total_data, Years == "2020"), 
                        aes(x =Fungal_green_Di , y = Fungal_field_Di, 
                            label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.y.npc = "top") + 
  theme(axis.text.y = element_blank(),
        plot.title = element_textbox(size = 14, color = "black", fill = "#E1BE6A",
          box.color = "grey50",padding = margin(5, 5, 5, 5), margin = margin(b = 0),       
          halign = 0.5, width = grid::unit(1, "npc"))) + 
  labs(y = NULL,
       x = "Composition distinctiveness estimated in greenhouse experiment", title = "2020") + 
  #annotate("text", x = 0.68, y = 0.60, label = bquote(atop(italic(r) == -0.17, italic(p) == 0.054)), size = 4) + 
  scale_color_manual(values = c("2018" = "#40B0A6", "2020" = "#E1BE6A", "2021" = "#A38E89")) + 
  scale_fill_manual(values = c("2018" = "#40B0A6", "2020" = "#E1BE6A", "2021" = "#A38E89")) + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.56,0.95)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.56,0.95)) -> p2b;p2b


ggplot() + 
  geom_star(data = subset(total_data, Years == "2021"), aes(x =Fungal_green_Di , y = Fungal_field_Di, color = Years, fill = Years, starshape = Site), 
            size = 1.2, alpha = 0.6, show.legend = F) + 
  scale_starshape_manual(values = c(13,11,23,15,5,28)) +
  geom_segment(subset(total_data2, Years == "2021"), mapping = aes(x = Green_Di, y = Field_Di, xend = Fungal_green_Di, yend = Fungal_field_Di),
               linetype = 1, linewidth = 0.2, alpha = 0.2, color = "#A38E89")+
  geom_errorbar(data = subset(total_data_mean, Years == "2021"), mapping = aes(x = Green_Di, ymax = Field_Di+Field_Di_se, ymin=Field_Di-Field_Di_se), width=0.008, alpha = 1, color = "black")+
  geom_errorbarh(data = subset(total_data_mean, Years == "2021"), mapping = aes(y = Field_Di, xmax = Green_Di+Green_Di_se, xmin=Green_Di-Green_Di_se), height=0.008, alpha = 1, color = "black")+
  geom_star(data = subset(total_data_mean, Years == "2021"), mapping = aes(x = Green_Di, y = Field_Di, fill = Years, starshape = Site), size=2.5, color = "black")+
  geom_abline(intercept = 0, slope = 1 , linetype = 2, color = "black")+
  mytheme + 
  theme(axis.text.y = element_blank(),
        plot.title = element_textbox(size = 14, color = "black", fill = "#A38E89",
                                     box.color = "grey50",padding = margin(5, 5, 5, 5), margin = margin(b = 0),       
                                     halign = 0.5, width = grid::unit(1, "npc"))) + 
  ggpmisc::stat_poly_eq(data = subset(total_data, Years == "2021"), 
                        aes(x =Fungal_green_Di , y = Fungal_field_Di, 
                            label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.y.npc = "top") + 
  labs(x = NULL, y = NULL, title = "2021")+ 
  #annotate("text", x = 0.68, y = 0.60, label = bquote(atop(italic(r) == -0.08, italic(p) == 0.368)), size = 4) + 
  scale_color_manual(values = c("2018" = "#40B0A6", "2020" = "#E1BE6A", "2021" = "#A38E89")) + 
  scale_fill_manual(values = c("2018" = "#40B0A6", "2020" = "#E1BE6A", "2021" = "#A38E89")) + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.56,0.95)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.56,0.95)) -> p2c;p2c

(p2a|p2b|p2c) -> Fig_3a; Fig_3a


################################## Figure 3b ###################################
# Environmental effects
total_data$Effect_size <- log(total_data$Fungal_field_Di/total_data$Fungal_green_Di)

total_data$Effect_change <- (total_data$Fungal_field_Di-total_data$Fungal_green_Di)/total_data$Fungal_green_Di

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
print(subset(final_result_t, p_value > 0.05))


# Mean ± 1SE
head(total_data)
mean_size <- total_data %>% dplyr::group_by(Latitude, Years) %>%
  dplyr::summarise(Effect_mean = mean(Effect_size, na.rm = TRUE),
                   Effect_se = sd(Effect_size, na.rm = TRUE) / sqrt(n()),
                   Effect_Funct_DI = mean(Fun_Di_log, na.rm = TRUE),
                   Effect_Phylo_DI = mean(Phy_Di_log, na.rm = TRUE),
                   plant_richness = n()) %>%
  dplyr::left_join(final_result_t, by = c("Latitude", "Years"))

cor.test(mean_size$Effect_mean, mean_size$Effect_Funct_DI)
cor.test(mean_size$Effect_mean, mean_size$Effect_Phylo_DI)
cor.test(mean_size$Latitude, mean_size$plant_richness)

colnames(total_data)
mean_change <- total_data %>% dplyr::group_by(Latitude, Years) %>%
  dplyr::summarise(Change_mean = mean(Effect_change, na.rm = TRUE),
                   Change_se = sd(Effect_change, na.rm = TRUE) / sqrt(n())) %>%
  dplyr::left_join(final_result_t, by = c("Latitude", "Years"))


mean_size$Years <- as.factor(mean_size$Years)
mean_size$p_value <- round(mean_size$p_value, 3)
mean_size$sig <- ifelse(mean_size$p_value > 0.05, 0, 1)
mean_size$sig <- as.factor(mean_size$sig)
total_data$Years <- as.factor(total_data$Years)

bestFitM2(data= mean_size, x= "Latitude", y = "Effect_mean")
FitM(data= mean_size, x= "Latitude", y = "Effect_mean", model = "line2P")

# plot
ggplot() + 
  geom_smooth(data = mean_size, mapping = aes(x = Latitude, y = Effect_mean),
              method = "lm", se = T, color = "black", fill = "grey80") +
  geom_errorbar(data = subset(mean_size, Years == 2018),mapping = aes(x = Latitude-0.10,ymax = Effect_mean+Effect_se, ymin=Effect_mean-Effect_se, color = Years),width=0.2,alpha = 1, color = "black")+
  geom_errorbar(data = subset(mean_size, Years == 2020),mapping = aes(x = Latitude,ymax = Effect_mean+Effect_se, ymin=Effect_mean-Effect_se, color = Years),width=0.2,alpha = 1, color = "black")+
  geom_errorbar(data = subset(mean_size, Years == 2021),mapping = aes(x = Latitude+0.10,ymax = Effect_mean+Effect_se, ymin=Effect_mean-Effect_se, color = Years), width = 0.2, alpha = 1, color = "black")+
  geom_star(data = subset(mean_size, Years == 2018 & Latitude != 23.1), mapping = aes(x = Latitude-0.10, y = Effect_mean, starshape = factor(Latitude), fill = Years), size = 2.8, show.legend = F) + 
  geom_star(data = subset(mean_size, Years == 2018 & Latitude == 23.1), mapping = aes(x = Latitude-0.10, y = Effect_mean),size = 2.8, show.legend = F, color = "#40B0A6", fill = "white", starshape = 13) + 
  geom_star(data = subset(mean_size, Years == 2020 & Latitude != 23.1), mapping = aes(x = Latitude, y = Effect_mean, starshape = factor(Latitude), fill = Years),size = 2.8, show.legend = F) + 
  geom_star(data = subset(mean_size, Years == 2020 & Latitude == 23.1), mapping = aes(x = Latitude, y = Effect_mean),size = 2.8, show.legend = F, color = "#E1BE6A", fill = "white", starshape = 13) + 
  geom_star(data = subset(mean_size, Years == 2021 & Latitude != 27.9 & Latitude != 30.5), mapping = aes(x = Latitude+0.10, y = Effect_mean, starshape = factor(Latitude), fill = Years),size = 2.8, show.legend = F) + 
  geom_star(data = subset(mean_size, Years == 2021 & Latitude == 27.9), mapping = aes(x = Latitude+0.10, y = Effect_mean),size = 2.8, show.legend = F, color = "#A38E89", fill = "white", starshape = 23) + 
  geom_star(data = subset(mean_size, Years == 2021 & Latitude == 30.5), mapping = aes(x = Latitude+0.10, y = Effect_mean),size = 2.8, show.legend = F, color = "#A38E89", fill = "white", starshape = 15) + 
  ggpmisc::stat_poly_eq(data = mean_size, aes(x =Latitude , y = Effect_mean, 
                            label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.x = 0.1, label.y = 0.2) + 
  scale_starshape_manual(values = c("23.1" = 13, "25.2" = 11, "27.9" = 23, "30.5" = 15, "34.6" = 5, "36.2" = 28)) +
  #scale_color_manual(values = c("2018" = "#2D6FA6", "2020" = "#CA5516", "2021" = "#40A47B")) + 
  #scale_fill_manual(values = c("2018" = "#2D6FA6", "2020" = "#CA5516", "2021" = "#40A47B")) + 
  #scale_color_manual(values = c("2018" = "#40B0A6", "2020" = "#E1BE6A", "2021" = "#A38E89")) + 
  scale_fill_manual(values = c("2018" = "#40B0A6", "2020" = "#E1BE6A", "2021" = "#A38E89")) + 
  geom_hline(yintercept = 0, linetype = 2) +
  mytheme + theme(legend.position = "right") +
  annotate("text", x = 37.7, y = 0.185, label = "Heterogenization", angle = 270, size = 4, hjust = 0) +
  annotate("text", x = 37.7, y = -0.225, label = "Homogenization", angle = 270, size = 4, hjust = 1) +
  theme(legend.position = c(0.35,0.25)) + 
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) + 
  labs(x = "Latitude (North degrees)", 
       y = bquote(atop("Environmental effects", 
                       Ln ~ "(" ~ frac(Fungi-dist["estimated in field survey"], 
                                       Fungi-dist["estimated in greenhouse experiment"]) ~ ")")),
       tag = "(b)") +
  annotate("segment",x = 37.3, xend = 37.3, y = 0, yend = 0.2,  
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"),color = "#A5CCC9", size = 2) +
  annotate("segment",x = 37.3, xend = 37.3, y = 0, yend = -0.3,  
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"),color = "#D48286", size = 2) -> Fig_3b; Fig_3b
  
# combination
Fig_3b <- Fig_3b + plot_layout(widths = c(0.8,0.2))

Fig_3a/(Fig_3b) + plot_layout(heights = c(0.46,0.54))

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.

