The following files allow one to reproduce analyses in the manuscript entitled "Alien invasion effect on soil biota: the interplay between alien species and environmental variation".

DATA & FILE OVERVIEW

***In Data folder***

The experimental data are stored in Figshare [![DOI](https://zenodo.org/badge/DOI/10.6084/m9.figshare.27880494.svg)](https://doi.org/10.6084/m9.figshare.27880494.v1).
Before the manuscript is officially published, experimental and analytical data must remain confidential. 
If needed, please contact the first or corresponding author in advance to obtain the relevant experimental data. 
All data will be made available upon acceptance of the manuscript.

*List of experimental data files (.xlsx)*

    * 1. Field_data_row_ASVs.xlsx  
    * 2. Field_data_group.xlsx  
    * 3. Field_fungi_Flattening.xlsx
    * 4. Greenhouse_data_row_ASVs.xlsx  
    * 5. Greenhouse_data_group.xlsx  
    * 6. Greenhouse_fungi_Flattening.xlsx  
    * 7. Site_information_seed.xlsx  
    * 8. traits_mean.xlsx  
    
*List of phylogenetic tree data files (.newick)*  

    * 1. IQ_tree_plant_2025.newick  (Phylogenetic tree of studying plant species)

***In Code folder***

The names of R-scripts correspond to the statistical analysis and visualization of the corresponding figures or tables in this manuscript.

*List of R-scripts*

    * 1. Figure_1.tiff (Draw in ppt)
    * 2. Figure_2.R  
    * 3. Figure_3 & Table S6.R  
    * 4. Figure_S1.R  
    * 5. Figure_S2.R  
    * 6. Figure_S3.R 
    * 7. Figure_S4.R  
    * 8. Table S3 & Table S4 & Table S5 (Field survey part).R  
    * 9. Table S3 & Table S4 & Table S5 (Greenhouse exp. part).R  
    
**Data-specific onformation for:** ***Field_data_row_ASVs.xlsx***

    * Abundance table of raw sequencing data of rhizosphere fungi from field survey (not rarefied to minimum sample size).

**Data-specific onformation for:** ***Field_data_group.xlsx***

    Variable list	         Description
    * Sample_ID	         Sample id of  plant rhizosphere soil 
    * Latitude	         Latitude of sampling point
    * Longitude	         Longitude of sampling point
    * Year	                 Sampling year
    * Site	                 Name of sampling site
    * Chinese_name	         Chinese name of study species
    * Species	         Latin name of study species
    * Genus	               Genus name of research species
    * Family	            Family name of research species
    * Origin	            Geographical origin of plants (native vs. exotic)
    * Site_pool	            The total richness of fungi for each site and year
    * Soil_ph	               Soil pH
    * Wcont	                  Soil water content
    * Soil_N	               Soil total nitrogen content
    * Tave	                  Annual average temperature (℃)
    * Precipitation	         Annual precipitation (mm)
    * Fungal_SR	            Fungi richness
    * Chol	                  Leaf chlorophyll (SPAD)
    * SLA	                  Specific leaf area (cm2 g-1)
    * LDMC	                  Leaf dry matter content (g g-1)
    * SRL	                  Specific root length (cm2 g-1)
    * FRR	                  Fine-to-total root mass (g g-1)
    * RS	                     Root-to-shoot mass ratio (g g-1)
    * Fungal_field_Di	 Fungal compositional distinctiveness estimated in the field
    * Fungal_green_Di	 Fungal compositional distinctiveness estimated in the greenhouse experiment
    * Fun_Di	 Species functional distinctiveness for each year and site
    * Phy_Di	 Species phylogenetic distinctiveness for each year and site

**Data-specific onformation for:** ***Field_fungi_Flattening.xlsx***

    * Abundance table of raw sequencing data of rhizosphere fungi from field survey (rarefied to minimum sample size)
      
**Data-specific onformation for:** ***Greenhouse_ASVs_row_data.xlsx***

    * Abundance table of raw sequencing data of rhizosphere fungi from greenhouse experiment (not rarefied to minimum sample size)

**Data-specific onformation for:** ***Greenhouse_data_group.xlsx***

    Variable list	 Description
    * Sample_ID	 Sample id of  plant rhizosphere soil 
    * Repeats	 Repeat number of soil samples
    * Chinese_name	 Chinese name of study species
    * Species	 Latin name of study species
    * Genus	 Genus name of research species
    * Family	 Family name of research species
    * Origin	 Geographical origin of plants (native vs. exotic)
    * Overall_Richness	 Fungi richess

**Data-specific onformation for:** ***Greenhouse_fungi_Flattening.xlsx***

    * Abundance table of raw sequencing data of rhizosphere fungi from greenhouse experiment (rarefied to minimum sample size)
    
**Data-specific onformation for:** ***Site_information_seed.xlsx***

    Sheets: sp_site	Information on the location of seed collection of studying species
    Variable list	 Description
    * Species_powo	 Latin species (powo)
    * Species	 Latin species
    * Site1	 Collection site 1
    * Site2	 Collection site 2
    * Site3	 Collection site 3
    
    Sheets: site_map	City information of studying species seed collection point
    Variable list	 Description
    * Pinyin	 City name
    * Address	 City name
    * City	 Plant richness of native species (number of species added)

**Data-specific onformation for:** ***traits_mean.xlsx***

    Variable list	 Description
    * Chol	 Leaf chlorophyll (SPAD)
    * SLA	 Specific leaf area (cm2 g-1)
    * LDMC	 Leaf dry matter content (g g-1)
    * SRL	 Specific root length (cm2 g-1)
    * FRR	 Fine-to-total root mass (g g-1)
    * RS	 Root-to-shoot mass ratio (g g-1)
    * Origin	 Geographical origin of plants (native vs. exotic)
    * taxon	 Latin name of study species
    * Species	 Latin name of study species
