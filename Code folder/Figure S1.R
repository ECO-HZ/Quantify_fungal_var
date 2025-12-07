################################################################################
################################# Figure S1 ####################################
################################################################################
library(ggplot2)
library(openxlsx)
library(emmeans)

# Custom style
mytheme = theme(
  legend.position = "none",
  panel.grid=element_blank(), 
  strip.background = element_rect(color="black", fill="white", size=0.5, linetype="solid"),
  strip.text.x = element_text(size = 9, color = "black"), # face = "bold.italic"
  legend.title = element_blank(),
  legend.key = element_blank(),
  legend.text = element_text(size = 10),
  legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
  axis.ticks = element_line(color='black'),
  axis.line = element_line(colour = "black"), 
  axis.title.x = element_text(colour='black', size=14),
  axis.title.y = element_text(colour='black', size=14),
  axis.text = element_text(colour='black',size=12),
  plot.title = element_textbox(
    size = 14, color = "black", fill = "grey90",
    box.color = "grey50",padding = margin(5, 5, 5, 5), margin = margin(b = 0),       
    halign = 0.5, width = grid::unit(1, "npc")), #r = unit(3, "pt")     
  plot.tag = element_text(size = 16, face = "bold")) 

# Loading the grouping metadata of soil samples
Field_group <- read.xlsx("Field_data_group.xlsx", sheet = "Field_group", rowNames = T, colNames = T)
Field_group$Sample_ID <- rownames(Field_group)

# Data Transformation
Field_group$SRL <- log10(Field_group$SRL)
Field_group$Wcont <- sqrt(Field_group$Wcont*100)
Field_group$Soil_N <- sqrt(Field_group$Soil_N)
Field_group$Years <- as.factor(Field_group$Years)
Field_group$Site <- factor(Field_group$Site, levels = c("Guangzhou","Guilin","Changsha","Wuhan","Zhengzhou","Tai'an"))
Field_group$Origin <- factor(Field_group$Origin, levels = c("Native","Exotic"))

## annual average temperature 
data_Tave <- unique(Field_group[,c("Years" ,"Site", "Latitude", "Tave")])

mod1 = lm(Tave ~ Years * Latitude, data_Tave)
anova(mod1)
mod1_emtrends <- emtrends(mod1, pairwise ~ Years, var = "Latitude")
test(mod1_emtrends, adjust = "BH")

ggplot() +
  geom_point(data_Tave, mapping = aes(x = Latitude, y = Tave, fill = Years),
             position = position_dodge(width = 0.4), size = 2.5, pch = 21) + 
  geom_smooth(data_Tave, mapping = aes(x = Latitude, y = Tave, color = Years),
              method = "lm", formula = y ~ x, se = F, size = 0.8, linetype = 1, alpha = 1) + 
  ggpmisc::stat_poly_eq(data_Tave, mapping = aes(x = Latitude, y = Tave, color = Years, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
                        formula = y ~ x, parse = TRUE, size = 4) + 
  geom_smooth(data = data_Tave, mapping = aes(x = Latitude, y = Tave), color = "black", 
              method = "lm", formula = y ~ x, se = F, linetype = 1, size = 1.2) + 
  ggpmisc::stat_poly_eq(data = data_Tave,  mapping = aes(x = Latitude, y = Tave, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.y.npc = "center") + 
  theme_bw() + mytheme + 
  theme(legend.position = c(0.85,0.85)) + 
  scale_color_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_fill_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  labs(x = NULL, y = expression("Annual average temperature ( " * degree * "C)"), tag = "(a)") -> p1; p1


## annual Precipitation
data_Precipitation <- unique(Field_group[,c("Years" ,"Site", "Latitude", "Precipitation")])
# anova(lm(Precipitation ~ Latitude, subset(data_Precipitation, Years == "2020")))
mod2 = lm(Precipitation ~ Years * Latitude, data_Precipitation)
anova(mod2)
mod2_emtrends <- emtrends(mod2, pairwise ~ Years, var = "Latitude")
test(mod2_emtrends, adjust = "BH")

ggplot() +
  geom_point(data_Precipitation, mapping = aes(x = Latitude, y = Precipitation, fill = Years),
             position = position_dodge(width = 0.4), size = 2.5, pch = 21) + 
  geom_smooth(data_Precipitation, mapping = aes(x = Latitude, y = Precipitation, color = Years, linetype = Years),
              method = "lm", formula = y ~ x, se = F, size = 0.8, alpha = 1) + 
  ggpmisc::stat_poly_eq(data_Precipitation, mapping = aes(x = Latitude, y = Precipitation, color = Years, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
                        formula = y ~ x, parse = TRUE, size = 4) + 
  geom_smooth(data = data_Precipitation, mapping = aes(x = Latitude, y = Precipitation), color = "black", 
              method = "lm", formula = y ~ x, se = F, linetype = 1, size = 1.2) + 
  ggpmisc::stat_poly_eq(data = data_Precipitation,  mapping = aes(x = Latitude, y = Precipitation, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.y.npc = "center") + 
  theme_bw() + mytheme +
  scale_color_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_fill_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_linetype_manual(values = c(1,2,2)) + 
  #scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  #scale_x_continuous(breaks=c(23.1,25.2,27.9,30.5,34.6,36.2)) + 
  labs(x = NULL, y = expression("Annual precipitation (mm)"), tag = "(b)") -> p2; p2


## Soil pH
data_Soil_ph <- Field_group[,c("Origin", "Years" ,"Site", "Latitude", "Soil_ph")]

mod3 = lm(Soil_ph ~ Years * Latitude, data = data_Soil_ph)
anova(mod3)
mod3_emtrends <- emtrends(mod3, pairwise ~ Years, var = "Latitude")
test(mod3_emtrends, adjust = "BH")

ggplot() +
  geom_point(data_Soil_ph, mapping = aes(x = Latitude, y = Soil_ph, fill = Years, group = Years),
             position = position_dodge(width = 0.4), size = 2.5, pch = 21) +  
  geom_smooth(data_Soil_ph, mapping = aes(x = Latitude, y = Soil_ph, color = Years, linetype = Years),
              method = "lm", formula = y ~ x, se = F, size = 0.8, alpha = 1) + 
  ggpmisc::stat_poly_eq(data_Soil_ph, mapping = aes(x = Latitude, y = Soil_ph, color = Years, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
                        formula = y ~ x, parse = TRUE, size = 4) +
  geom_smooth(data = data_Soil_ph, mapping = aes(x = Latitude, y = Soil_ph), color = "black", 
              method = "lm", formula = y ~ x, se = F, linetype = 1, size = 1.2) + 
  ggpmisc::stat_poly_eq(data = data_Soil_ph,  mapping = aes(x = Latitude, y = Soil_ph, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.y.npc = "center") + 
  theme_bw() + mytheme + 
  scale_color_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_fill_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_linetype_manual(values = c(2,1,1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = NULL, y = "Soil pH", tag = "(c)") -> p3; p3


## Soil total nitrogen content
data_Soil_N <- Field_group[,c("Origin", "Years" ,"Site", "Latitude", "Soil_N")]

mod4 = lm(Soil_N ~ Years * Latitude, data = data_Soil_N)
shapiro.test(residuals(mod4))
mod4_emtrends <- emtrends(mod4, pairwise ~ Years, var = "Latitude")
test(mod4_emtrends, adjust = "BH")

ggplot() +
  geom_point(data_Soil_N, mapping = aes(x = Latitude, y = Soil_N, fill = Years, group = Years),
             position = position_dodge(width = 0.4), size = 2.5, pch = 21) +  
  geom_smooth(data_Soil_N, mapping = aes(x = Latitude, y = Soil_N, color = Years),
              method = "lm", formula = y ~ x, se = F, size = 0.8, linetype = 1, alpha = 1) + 
  ggpmisc::stat_poly_eq(data_Soil_N, mapping = aes(x = Latitude, y = Soil_N, color = Years, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
                        formula = y ~ x, parse = TRUE, size = 4) +
  geom_smooth(data = data_Soil_N, mapping = aes(x = Latitude, y = Soil_N), color = "black", 
              method = "lm", formula = y ~ x, se = F, linetype = 1, size = 1.2) + 
  ggpmisc::stat_poly_eq(data = data_Soil_N,  mapping = aes(x = Latitude, y = Soil_N, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.y.npc = "center") + 
  theme_bw() + mytheme + 
  scale_color_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_fill_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  labs(x = NULL, y = "Soil total nitrogen content (%, sqrt)", tag = "(d)") -> p4; p4


## Soil water content
data_Wcont <- Field_group[,c("Origin", "Years" ,"Site", "Latitude", "Wcont")]

mod5 = lm(Wcont ~ Origin * Years * Latitude, data = data_Wcont)
shapiro.test(residuals(mod5))
mod5_emtrends <- emtrends(mod5, pairwise ~ Years, var = "Latitude")
test(mod5_emtrends, adjust = "BH")

ggplot() +
  geom_point(data_Wcont, mapping = aes(x = Latitude, y = Wcont, fill = Years, group = Years),
             position = position_dodge(width = 0.4), size = 2.5, pch = 21) +  
  #geom_errorbar(data_Wcont_mean, mapping = aes(x = Latitude, y = mean_Wcont, ymin = mean_Wcont - se_Wcont, ymax = mean_Wcont + se_Wcont, color = Years, group = Years), 
  #              width = 0.5, size = 0.5, position = position_dodge(width = 0)) +
  geom_smooth(data_Wcont, mapping = aes(x = Latitude, y = Wcont, color = Years, linetype = Years),
              method = "lm", formula = y ~ x, se = F, size = 0.8) + 
  ggpmisc::stat_poly_eq(data_Wcont, mapping = aes(x = Latitude, y = Wcont, color = Years, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
                        formula = y ~ x, parse = TRUE, size = 4) +
  ## all dataset
  geom_smooth(data = data_Wcont, mapping = aes(x = Latitude, y = Wcont), color = "black", 
              method = "lm", formula = y ~ x, se = F, linetype = 1, size = 1.2) + 
  ggpmisc::stat_poly_eq(data = data_Wcont,  mapping = aes(x = Latitude, y = Wcont, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
                        formula = y ~ x, parse = TRUE, size = 4, color = "black", label.y.npc = "center") + 
  theme_bw() + mytheme + 
  scale_color_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_fill_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_linetype_manual(values = c(1,1,2)) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  #scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  #scale_x_continuous(breaks=c(23.1,25.2,27.9,30.5,34.6,36.2)) + 
  labs(x = NULL, y = "Soil water content (%, sqrt)", tag = "(e)") -> p5; p5

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.
