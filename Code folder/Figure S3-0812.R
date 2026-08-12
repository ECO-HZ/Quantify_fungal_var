################################################################################
################################# Figure S3 ####################################
################################################################################
library(ggplot2)
library(openxlsx)
library(emmeans)
library(ggeffects)

# Custom style
mytheme = theme_classic()+ 
  theme(
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

mod1 = lmer(Tave ~ Years * Latitude + (1|Site), data_Tave)
anova(mod1, ddf = "Kenward-Roger")
shapiro.test(residuals(mod1))

pred_mod1 <- ggeffect(mod1, terms = c("Latitude"))
eff_mod1_data <- data.frame(pred_mod1)
colnames(eff_mod1_data)[1] <- "Latitude"
colnames(eff_mod1_data)[2] <- "Tave"

ggplot() +
  geom_point(data_Tave, mapping = aes(x = Latitude, y = Tave, fill = Years),
             position = position_dodge(width = 0.4), size = 2.5, pch = 21) + 
  geom_line(data = eff_mod1_data, mapping = aes(Latitude, Tave), size = 1.2) +
  mytheme + 
  theme(legend.position = c(0.85,0.85)) + 
  scale_color_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_fill_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  labs(x = NULL, y = expression("Annual average temperature ( " * degree * "C)"), tag = "(a)") -> p1; p1


## annual Precipitation
data_Precipitation <- unique(Field_group[,c("Years" ,"Site", "Latitude", "Precipitation")])
# anova(lm(Precipitation ~ Latitude, subset(data_Precipitation, Years == "2020")))
mod2 = lmer(Precipitation ~ Years * Latitude + (1|Site), data_Precipitation)
anova(mod2, ddf = "Kenward-Roger")

pred_mod2 <- ggeffect(mod2, terms = c("Latitude"))
eff_mod2_data <- data.frame(pred_mod2)
colnames(eff_mod2_data)[1] <- "Latitude"
colnames(eff_mod2_data)[2] <- "Precipitation"

ggplot() +
  geom_point(data_Precipitation, mapping = aes(x = Latitude, y = Precipitation, fill = Years),
             position = position_dodge(width = 0.4), size = 2.5, pch = 21) + 
  geom_line(data = eff_mod2_data, mapping = aes(Latitude, Precipitation), size = 1.2) +
  theme_bw() + mytheme +
  scale_color_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_fill_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_linetype_manual(values = c(1,2,2)) + 
  #scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  #scale_x_continuous(breaks=c(23.1,25.2,27.9,30.5,34.6,36.2)) + 
  labs(x = NULL, y = expression("Annual precipitation (mm)"), tag = "(b)") -> p2; p2


## Soil pH
data_Soil_ph <- Field_group[,c("Origin", "Years" ,"Site", "Latitude", "Soil_ph", "Species")]

mod3 = lmer(Soil_ph ~ Years * Latitude + (1|Site) + (1|Site:Years), data = data_Soil_ph)
anova(mod3, ddf = "Kenward-Roger")

pred_mod3 <- ggeffect(mod3, terms = c("Latitude"))
eff_mod3_data <- data.frame(pred_mod3)
colnames(eff_mod3_data)[1] <- "Latitude"
colnames(eff_mod3_data)[2] <- "Soil_ph"

ggplot() +
  geom_point(data_Soil_ph, mapping = aes(x = Latitude, y = Soil_ph, fill = Years, group = Years),
             position = position_dodge(width = 0.4), size = 2.5, pch = 21) +  
  geom_line(data = eff_mod3_data, mapping = aes(Latitude, Soil_ph), size = 1.2, linetype = 2) +
  theme_bw() + mytheme + 
  scale_color_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_fill_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_linetype_manual(values = c(2,1,1)) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = NULL, y = "Soil pH", tag = "(c)") -> p3; p3


## Soil total nitrogen content
data_Soil_N <- Field_group[,c("Origin", "Years" ,"Site", "Latitude", "Soil_N", "Species")]

mod4 = lmer(Soil_N ~ Years * Latitude + (1|Site) + (1|Site:Years), data = data_Soil_N)
anova(mod4, ddf = "Kenward-Roger")

pred_mod4 <- ggeffect(mod4, terms = c("Latitude"))
eff_mod4_data <- data.frame(pred_mod4)
colnames(eff_mod4_data)[1] <- "Latitude"
colnames(eff_mod4_data)[2] <- "Soil_N"

ggplot() +
  geom_point(data_Soil_N, mapping = aes(x = Latitude, y = Soil_N, fill = Years, group = Years),
             position = position_dodge(width = 0.4), size = 2.5, pch = 21) +  
  geom_line(data = eff_mod4_data, mapping = aes(Latitude, Soil_N), size = 1.2, linetype = 2) +
  theme_bw() + mytheme + 
  scale_color_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_fill_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  #scale_y_continuous(expand = expansion(mult = c(0.1, 0.3))) + 
  labs(x = NULL, y = "Soil total nitrogen content (%, sqrt)", tag = "(d)") -> p4; p4


## Soil water content
data_Wcont <- Field_group[,c("Origin", "Years" ,"Site", "Latitude", "Wcont", "Species")]

mod5 = lmer(Wcont ~ Years * Latitude + (1|Site) + (1|Site:Years), data = data_Wcont)
shapiro.test(residuals(mod5))
anova(mod5, ddf = "Kenward-Roger")

pred_mod5 <- ggeffect(mod5, terms = c("Latitude"))
eff_mod5_data <- data.frame(pred_mod5)
colnames(eff_mod5_data)[1] <- "Latitude"
colnames(eff_mod5_data)[2] <- "Wcont"

ggplot() +
  geom_point(data_Wcont, mapping = aes(x = Latitude, y = Wcont, fill = Years, group = Years),
             position = position_dodge(width = 0.4), size = 2.5, pch = 21) +  
  geom_line(data = eff_mod5_data, mapping = aes(Latitude, Wcont), size = 1.2, linetype = 2) +
  theme_bw() + mytheme + 
  scale_color_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) + 
  scale_fill_manual(values = c("2018" = "#898EA1", "2020" = "#CF9742", "2021" = "#3A7C72")) +
  scale_linetype_manual(values = c(1,1,2)) + 
  #scale_y_continuous(expand = expansion(mult = c(0.1, 0.3))) + 
  #scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  #scale_x_continuous(breaks=c(23.1,25.2,27.9,30.5,34.6,36.2)) + 
  labs(x = NULL, y = "Soil water content (%, sqrt)", tag = "(e)") -> p5; p5

# (p1/p4|p2/p5|p3/p5)

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.
