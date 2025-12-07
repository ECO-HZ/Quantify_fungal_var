################################################################################
################################## Figure S2 ###################################
################################################################################

# Loading the R packages
library(openxlsx)
library(dplyr)
library(tidygeocoder)
library(sf)
library(ggplot2)
library(ggspatial)
library(gt)
library(tidyverse)
library(glue)
library(ggtree)
library(phytools)

# Loading seed collection location information
seed_data <- read.xlsx("Site_information_seed.xlsx", sheet = "site_map", rowNames = F, colNames = T)

# Obtain latitude and longitude information of the collection location
options(amap_key = '989a66c96e4750756c11390fa782d8bc')
seed_site_data <- seed_data %>% 
  geocode("Address", method = 'arcgis', lat = latitude, 
          long = longitude, full_results = TRUE) %>% as.data.frame()
#head(seed_site_data)
colnames(seed_site_data)[c(4,5)] <- c("Latitude", "Longitude")
seed_site_data <- seed_site_data %>% left_join(seed_data) 

# Loading China map
china_map <- sf::st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json") 

################################# Figure S2a ###################################
ggplot(china_map)+
  geom_sf(data = china_map, fill = "grey95",size = 1) + 
  xlim(105, 122)+ ylim(18, 42)+ 
  ggnewscale::new_scale_fill() + 
  annotate("text", x = 106, y = 32, label = "China", size = 4, fontface = "bold") + 
  geom_point(data = seed_site_data, aes(x = Longitude, y = Latitude),
             size = 3, alpha = 1, shape = 21, color = "black", fill = "#32B0A5") + 
  annotation_scale(location = "br", style = "ticks", line_width = 1.2, pad_y = unit(0.5, "cm"), text_cex = 0.9) + 
  annotation_north_arrow(location = "tl", which_north = "true", height = unit(1, "cm"), width = unit(1, "cm"),
                         pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"), style = north_arrow_fancy_orienteering) +
  theme_bw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(size = 0.5, colour = "NA"),
        axis.line.y = element_line(size = 0.5, colour = "NA"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 11),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10),
        plot.tag = element_text(size = 14, face = "bold"))+
  labs(x = " ", y = " ", tag = "(a)") +
  ggrepel::geom_text_repel(mapping = aes(x = Longitude,y = Latitude,label = City), data = seed_site_data, size = 2.8,
                           segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 25))


################################# Figure S2b ###################################
tree <- read.newick("IQ_tree_plant_2025.NEWICK")
to_drop <- c("Amborella_trichopoda","")
tree <- drop.tip(as.phylo(tree), to_drop) 
plot(tree)
ggtree(tree, size = 0.4, color = "black", branch.length = "branch.length", ladderize = F) + 
  geom_tiplab(aes(label = sub("_", " ", label)), size = 2, offset = 0.01,
              fontface = "italic", linetype = "dotted", align = TRUE) + 
  xlim(0,8)-> p.tree; p.tree

ggtree::rotate(p.tree, node = 56)

# Add collection point information
seed_site_infor = read.xlsx("Site_information_seed.xlsx", sheet = "sp_site", rowNames = F, colNames = T)
rownames(seed_site_infor) = seed_site_infor$Species
seed_site_infor = seed_site_infor[rev(tree$tip.label), ]
seed_site_infor$Species = NULL
gt(seed_site_infor)

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.


