The following files allow one to reproduce analyses in the manuscript entitled "Environmental context modulates rhizosphere fungal distinctiveness across plant phylogeny and functional traits".

DATA & FILE OVERVIEW

***In Data folder***

The experimental data are stored in Figshare [![DOI](https://zenodo.org/badge/DOI/10.6084/m9.figshare.27880494.svg)](https://doi.org/10.6084/m9.figshare.27880494.v5).
Before the manuscript is officially published, experimental and analytical data must remain confidential. 
If needed, please contact the first or corresponding author in advance to obtain the relevant experimental data. 
All data will be made available upon acceptance of the manuscript.

*List of experimental data files (.xlsx)*

    * 1. Field_data_raw_ASVs.xlsx  
    * 2. Field_data_group.xlsx  
    * 3. Greenhouse_data_raw_ASVs.xlsx  
    * 4. Greenhouse_data_group.xlsx  
    * 5. Site_information_seed.xlsx  
    * 6. traits_mean.xlsx
    
*List of phylogenetic tree data files (.newick)*  

    * 1. IQ_tree_plant_2025.newick (Phylogenetic tree of studying plant species)

***In Code folder***

The names of R-scripts correspond to the statistical analysis and visualization of the corresponding figures or tables in this manuscript.

*List of R-scripts*

    * 1. Figure 1.pdf
    * 2. Figure 2 & Table S2.R  
    * 3. Figure 3.R  
    * 4. Figure 4 & Table S4.R  
    * 5. Figure S1.R  
    * 6. Figure S2.R 
    * 7. Figure S3.R  
    * 8. Figure S4.R  
    * 9. Figure S5.R  
    * 10. Figure S6.R  
    * 11. Figure S7 & Table S3.R  
    Ohter codes
    * 01_Calculation of asymptotic exponential Shannon diversity.R  
      All the samples, regardless of number of reads, were well characterized with saturating accumulation curves such that ASV diversity measured in the samples 
      was nearly identical to asymptotic exponential Shannon diversity estimated in iNEXT (R2 > 0.99 for both the field survey and greenhouse experiment, based on 
      general linear regression). Therefore, we used the raw ASVs abundance data for all the following analyses.
      
    * 02_Distinctiveness index calculation.R  
      We sequentially calculated distinctiveness indices for plant functional traits, phylogenetic relationships, and fungal community composition based on a common distinctiveness metric.
    
**Data-specific onformation for:** ***Field_data_raw_ASVs.xlsx***

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
    * Genus	                 Genus name of research species
    * Family	         Family name of research species
    * Origin	         Geographical origin of plants (native vs. exotic)
    * Site_pool	         The total richness of fungi for each site and year
    * Soil_ph	         Soil pH
    * Wcont	                 Soil water content
    * Soil_N	         Soil total nitrogen content
    * Tave	                 Annual average temperature (℃)
    * Precipitation	         Annual precipitation (mm)
    * Chol	                 Leaf chlorophyll (SPAD)
    * SLA	                 Specific leaf area (cm2 g-1)
    * LDMC	                 Leaf dry matter content (g g-1)
    * SRL	                 Specific root length (cm2 g-1)
    * FRR	                 Fine-to-total root mass (g g-1)
    * RS	                 Root-to-shoot mass ratio (g g-1)
    * Fungal_field_Di	 Fungal compositional distinctiveness estimated in the field
    * Fungal_green_Di	 Fungal compositional distinctiveness estimated in the greenhouse experiment
    * Funct_Di	         Species functional distinctiveness for each year and site
    * Phylo_Di	         Species phylogenetic distinctiveness for each year and site
    * PCoA1	          The first PCoA axis of rhizosphere fungal communities of all field soil samples
    * PCoA2	         The second PCoA axis of rhizosphere fungal communities of all field soil samples
    * Fungal_Di_field_com	          Fungal compositional distinctiveness estimated in the field based on shared ASVs
    * Fungal_Di_green_com	          Fungal compositional distinctiveness estimated in the greenhouse experiment based on shared ASVs
  
**Data-specific onformation for:** ***Greenhouse_data_raw_ASVs.xlsx***

    * Abundance table of raw sequencing data of rhizosphere fungi from greenhouse experiment (not rarefied to minimum sample size)

**Data-specific onformation for:** ***Greenhouse_data_group.xlsx***

    Variable list	      Description
    * Sample_ID	      Sample id of  plant rhizosphere soil 
    * Repeats	      Repeat number of soil samples
    * Chinese_name	      Chinese name of study species
    * Species	      Latin name of study species
    * Genus	              Genus name of research species
    * Family	      Family name of research species
    * Origin	      Geographical origin of plants (native vs. exotic)

**Data-specific onformation for:** ***Site_information_seed.xlsx***

    Sheets: sp_site-Information on the location of seed collection of studying species
    Variable list	 Description
    * Species_powo	 Latin species (powo)
    * Species	 Latin species
    * Site1	         Collection site 1
    * Site2	         Collection site 2
    * Site3	         Collection site 3
    
    Sheets: site_map-City information of studying species seed collection point
    Variable list	 Description
    * Pinyin	 City name
    * Address	 City name
    * City	         Plant richness of native species (number of species added)

**Data-specific onformation for:** ***traits_mean.xlsx***

    Variable list	 Description
    * Chol	         Leaf chlorophyll (SPAD)
    * SLA	         Specific leaf area (cm2 g-1)
    * LDMC	         Leaf dry matter content (g g-1)
    * SRL	         Specific root length (cm2 g-1)
    * FRR	         Fine-to-total root mass (g g-1)
    * RS	         Root-to-shoot mass ratio (g g-1)
    * Origin	 Geographical origin of plants (native vs. exotic)
    * taxon	         Latin name of study species
    * Species	 Latin name of study species
