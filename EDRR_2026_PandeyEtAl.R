

# title: "Optimizing trap deployment for the detection of non-native scolytines"
# author: Manoj Pandey


# ---- Libraries ----
library(dplyr) 
library(ggplot2)
library(forcats)
library(car) 
library(sf)
library(sp)
library(ggspatial)
library(RColorBrewer)
library(FedData)
library(raster)
library(terra)
library(base)
library(tidyr)#
library(tmap)
library(devtools)
library(prism)
library(reshape2)
library(rmarkdown)
library(knitr)
library(DT)
library(data.table)
library(purrr)
library(stringr)
library(vegan)
library(landscapemetrics)
library(exactextractr)
library(parallel)
library(data.table)
library(tigris)
library(tidycensus)
library(units)
library(lme4)
library(usmap)
library(maps)
library(corrplot)
library(GGally)
library(pbapply)
library(ggmap)
library(USA.state.boundaries)
library(furrr)
library(progressr)
library(qs)
library(ggtext)
library(MuMIn)
library(broom.mixed)
library(glmmTMB)
library(DHARMa)
library(emdi)
library(performance)
library(lubridate)
library(viridis)
library(scales)
library(ggpubr)
library(patchwork)
library(tibble)

writeLines(toBibtex(citation("DHARMa")), "DHARMa.bib")


 
# ---- Cleaning the EDRR dataset ----

  
  edrr<-read.csv(file="Data/Processed/edrr.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings="", strip.white = T) 
# Importing the Dataset , xyleborus intrusus was corrected to Xyleborus intrusus, and hypothenemus crudiae was corrected to Hypothenemus crudiae before importing

# Correcting species names
  edrr <- edrr %>%
    mutate(PestScientificName= dplyr::recode(PestScientificName,
                            "Crypturgus alutaceous" = "Crypturgus alutaceus", 
                            "Dryoxylon onoharaensum" = "Dryoxylon onoharaense",
                            "Ips tridens englemanni" = "Ips tridens engelmanni",
                            "Pityophthorus annectans" = "Pityophthorus annectens",
                            "Pseudopityophthorus pruinosis"= "Pseudopityophthorus pruinosus",
                            "Scierus annectans"= "Scierus annectens",
                            "Chaetophloeus pencillatus"= "Chaetophloeus penicillatus",
                           
    ))
  
  str(edrr) #checking the structure of the dataset
  anyNA(edrr)# checking if there are any NAs in the dataset
  colSums(is.na(edrr)) # count of total NAs in each column in the dataset
  edrr_clean<-edrr[-c(29943,29944,92279), -c(1,2, 10,15, 20, 22:31, 33)] #Removing the rows and columns with missing data
  unique(edrr_clean$PestScientificName)
  edrr_clean<-edrr_clean[edrr_clean$PestScientificName != "%", ]  #Removing the observations with missing Pest Scientific Names and have % instead

  edrr_species<-as.data.frame(unique(edrr_clean$PestScientificName))
  print(edrr_species)
  edrr_species <- edrr_species %>% rename(Species = `unique(edrr_clean$PestScientificName)`)
  edrr_species <- edrr_species[order(edrr_species$Species),]
  edrr_species <- as.data.frame(edrr_species)
  unique(edrr$StateCode)

  
  # write.csv(edrr_species,"Data/Saved/edrr_species_list.csv")  ### List of all species from the EDRR data
  
 
  # EDRR dataset also has Cerambycidae, so adding taxa info (Cerambycidae/Scolytinae), information (complete/missing) to filter out all taxa other than Scolytines
  # Taxa information was manually added/assigned to the edrr_species_list data, the new dataset with taxa information is edrr_species 
 
  

  edrr_species<-read.csv(file="Data/Saved/edrr_species.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings="", strip.white = T, row.names = 1) # edrr_species is the dataframe with taxa information
  table(edrr_species$Taxa)
  edrr_species_Cerambycidae<-edrr_species[edrr_species$Taxa=="Cerambycidae", ]
  edrr_species_Scolytinae<-edrr_species[edrr_species$Taxa=="Scolytinae", ]
  
  
  edrr_taxa_info <- left_join(edrr_clean, edrr_species, by="PestScientificName")

  # Checking
  edrr_taxa_info_aggregate<- edrr_taxa_info %>% group_by(PestScientificName) %>% summarise(Taxa=first(Taxa),Info=first(Info))
  print(edrr_taxa_info_aggregate, n=408)
  # 
  
# Note: SiteNumber is not unique to a location
# FIPS are unique to each site (1-3 traps per site)
# FIPS and TrapSite both are unique to each site

# Creating unique TrapSite corresponding to each unique geometry  # There are multiple TrapIDs for different years for same locations
  edrr_taxa_info$Coordinates<-paste(edrr_taxa_info$SiteLatitude, edrr_taxa_info$SiteLongitude)

  unique_sites<-edrr_taxa_info %>%
  distinct(Coordinates) %>%
  mutate(TrapSite=paste0("Site", row_number()))
  
 
  edrr_taxa_info<- edrr_taxa_info %>%
    left_join(unique_sites, by="Coordinates")


  edrr_Scolytinae<-edrr_taxa_info[edrr_taxa_info$Taxa=="Scolytinae", ]
  edrr_Cerambycidae<-edrr_taxa_info[edrr_taxa_info$Taxa=="Cerambycidae", ]
  length(unique(edrr_Scolytinae$PestScientificName))  # 344
  length(unique(edrr_Cerambycidae$PestScientificName))  # 63


  write.csv(edrr_taxa_info,"Data/Saved/edrr_taxa_info.csv", row.names = FALSE)
  #
  write.csv(edrr_Scolytinae,"Data/Saved/edrr_Scolytinae.csv", row.names = FALSE)

#### edrr_taxa_info is the cleaned EDRR dataset with taxa information
#### edrr_Scolytinae is the edrr_taxa_info with Scolytines only

  

  
# ---- Clean edrr_Scolytinae data ----
  
  edrr_Scolytinae<-read.csv(file="Data/Saved/edrr_Scolytinae.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings="", strip.white = T)

  colnames(edrr_Scolytinae)

  # edrr_Scolytinae <- edrr_Scolytinae %>%
  #   rename(
  #     Species = PestScientificName,
  #     SiteLoc = SiteLocation,
  #     Latitude = SiteLatitude,
  #     Longitude = SiteLongitude,
  #     CollDate = CollectionDate,
  #     CollID = CollectionID,
  #     )
  # colnames(edrr_Scolytinae)

  # Headers edited to keep it to 10 characters max so that its easier when working with the spatial data and extracting info
  # names(edrr_Scolytinae) <- substr(names(edrr_Scolytinae), 1, 10)



  # write.csv(edrr_Scolytinae,"Data/Saved/edrr_Scolytinae.csv", row.names = FALSE)

  # 
# ---- EXTRACTING PREDICTORS ---- 
  edrr_Scolytinae<-read.csv(file="Data/Saved/edrr_Scolytinae.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings="", strip.white = T) 
  
  
  length(unique(edrr_Scolytinae$Species))  # 344
  
# ---- NLCD Land cover ----
 
  ## Link to download Annual NLCD -->  https://www.mrlc.gov/data?f%5B0%5D=project_tax_term_term_parents_tax_term_name%3AAnnual%20NLCD  
    # Project --> Annual NLCD
       # Products --> Land Cover
          # Regions --> CONUS
             # Time Period --> 2016
  
  
   edrr_Scolytinae_sf<-st_as_sf(edrr_Scolytinae, coords = c("Longitude", "Latitude"),
                               crs="EPSG:4326")
  

  # Annual NLCD Land Cover 2016 (Coterminous United States)
   nlcd_raster <- rast("Data/Raw/Annual_NLCD_LndCov_2016_CU_C1V2.tif")
  
  #plot(nlcd_raster)
  edrr_Scolytinae_sf <- st_transform(edrr_Scolytinae_sf, st_crs(nlcd_raster))
  
  
  #Extracting land cover values at EDRR point locations
  edrr_Scolytinae_sf$landcover <- terra::extract(nlcd_raster, vect(edrr_Scolytinae_sf))[, 2]  
  # 
  # str(edrr_Scolytinae_sf$landcover)
  #  
  edrr_Scolytinae_sf$landcover<-as.factor(edrr_Scolytinae_sf$landcover)
  
  unique(edrr_Scolytinae_sf$landcover)  # NO Class 12/Perennial Ice/Snow 
  levels(edrr_Scolytinae_sf$landcover)
  #  
  #  
  #  # NLCD land cover class names
  # nlcd_classes <- c(
  #   "11" = "Open Water",
  #   "12" = "Perennial Ice/Snow",
  #   "21" = "Developed, Open Space",
  #   "22" = "Developed, Low Intensity",
  #   "23" = "Developed, Medium Intensity",
  #   "24" = "Developed, High Intensity",
  #   "31" = "Barren Land",
  #   "41" = "Deciduous Forest",
  #   "42" = "Evergreen Forest",
  #   "43" = "Mixed Forest",
  #   "52" = "Shrub/Scrub",
  #   "71" = "Grassland/Herbaceous",
  #   "81" = "Pasture/Hay",
  #   "82" = "Cultivated Crops",
  #   "90" = "Woody Wetlands",
  #   "95" = "Emergent Herbaceous Wetlands"
  # )
  # 
  # # Numeric land cover codes to descriptive names
  # edrr_Scolytinae_sf$LC_class <- nlcd_classes[as.character(edrr_Scolytinae_sf$landcover)]
  
  #   
  # 
  unique(edrr_Scolytinae_sf$landcover)
  sum(is.na(edrr_Scolytinae_sf$landcover))

  
  
  
  
  
  


# ---- Ecoregions ----
  
  ## Link to download Ecoregions data -->  https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states
     ## Level III Ecoregions Downloads --> GIS data
        ## US Level III Ecoregions shapefile with state boundaries 
  
 ecoregions<-st_read("Data/Raw/us_eco_l3_state_boundaries.shp")
 unique(ecoregions$NA_L3NAME)
 unique(ecoregions$US_L3NAME)
 unique(ecoregions$US_L3CODE)
  
 
  edrr_Scolytinae_sf<- st_transform(edrr_Scolytinae_sf, st_crs(ecoregions))
  
  edrr_Scolytinae_sf<-st_join(edrr_Scolytinae_sf, ecoregions) # default is st_intersects()
  
  # unique(edrr_Scolytinae_sf$US_L3NAME)
  # length(unique(edrr_Scolytinae_sf$US_L3NAME))  # 66 including NA

  
  edrr_Scolytinae$landcover <- edrr_Scolytinae_sf$landcover
  # edrr_Scolytinae$LC_class <- edrr_Scolytinae_sf$LC_class
  edrr_Scolytinae$Ecoregion3 <- edrr_Scolytinae_sf$US_L3NAME
 
  
  # Removing observations from Puerto Rico  and Alaska # NLCD data has the data for continental US only  # Therefore no land cover class is allocated
  
   edrr_Scolytinae<-edrr_Scolytinae %>% filter(!is.na(edrr_Scolytinae$landcover))
  
  write.csv(edrr_Scolytinae, "Data/Saved/edrr_Scolytinae_landcover_ecoregions_CONUS.csv", row.names = FALSE)
  
  # edrr_Scolytinae_landcover_ecoregions_CONUS.csv is the EDRR Scolytinae dataset with landcover and ecoregion information
  
  
 
 
   length(unique(edrr_Scolytinae$Species))  # 320
   
  

  
  
  
  
  
# ---- Landscape metrics ---- 
  library(landscapemetrics)
  # landscapemetrics::list_lsm()
  all_metrics<-list_lsm() 
  
  edrr<-read.csv(file="Data/Saved/edrr_Scolytinae_landcover_ecoregions_CONUS.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings="", strip.white = T) 
 
  edrr_sf <- st_as_sf(edrr, coords = c("Longitude", "Latitude"), crs = "EPSG:4326")
 
  nlcd_raster <- rast("Data/Raw/Annual_NLCD_LndCov_2016_CU_C1V2.tif")
  # plot(nlcd_raster)
  
 
  edrr_sf<-st_transform(edrr_sf, st_crs(nlcd_raster))
  # st_crs(nlcd_raster)
  # st_crs(edrr_sf)

# ---- Checking if NLCD landscape works with landscapemetrics ----
  
  # check_landscape(nlcd_raster)
  # check_landscape(nlcd)
  # 
  
  ## Checking the entire nlcd raster layer did not work possibly due to memory issues... the layer being too large
  ## Subsetting to check
  
  options(tigris_use_cache = TRUE)
  Ohio<-states(cb=TRUE) |>
    filter(STUSPS == "OH") |>
    st_transform(crs(nlcd_raster))
  
  
  subset_nlcd<-crop(nlcd_raster, Ohio)
  check_landscape(subset_nlcd)
  
  # layer       crs units   class n_classes OK
  # 1     1 projected     m integer        15  ✔
  # Everything looks OK
  

  
  
# ---- pland (land cover class area %) at class level buffer 1000m ----
  edrr_vect<-vect(edrr_sf)
  metrics<- pblapply(1:nrow(edrr_vect), function(i) {
    tryCatch({
      msd <- sample_lsm(
        landscape = nlcd_raster,
        y = edrr_vect[i],
        size = 1000,
        shape = "circle",
        what = c("lsm_c_pland")
      )
      msd$id <- i
      return( msd)
    }, error = function(e) {
      message(sprintf("Error at point %d: %s", i, e$message))
      NULL
    })
  })
  
 
  
  # Combine results 
  metrics_df <- do.call(rbind, metrics)
  
  # hist(metrics_df$percentage_inside)
  
  # name for wide columns
  metrics_df <- metrics_df %>% mutate(metric_name = paste0("c", class, "_", metric))
  
  # wide format
  metrics_df_wide <- metrics_df %>%
    dplyr::select(id,metric_name, value) %>%
    pivot_wider(
      names_from = metric_name,
      values_from = value,
      values_fill = 0
      
    )
  
  
  # write.csv(metrics_df_wide, "Data/Saved/landscapemetrics_dataframe.csv")
  
  
  metrics_df_wide[, 2:16] <- lapply(metrics_df_wide[, 2:16], as.numeric)
  str(metrics_df_wide)
  metrics_df_wide$total_pland <- rowSums(metrics_df_wide[, c(2:16)], na.rm = TRUE)
  
  
  forest_classes <- metrics_df_wide[, c(1,4,5)]
  
  # Joining with edrr
  
  edrr_with_landcover <- bind_cols(edrr, forest_classes[match(1:nrow(edrr), forest_classes$id), ])

  edrr_with_landcover <- bind_cols(edrr, metrics_df_wide[match(1:nrow(edrr), metrics_df_wide$id), ])
 
 

# ---- Shannon Diversity at landscape level ----
  metrics1<- pblapply(1:nrow(edrr_vect), function(i) {
    tryCatch({
      msd1 <- sample_lsm(
        landscape = nlcd_raster,
        y = edrr_vect[i],
        size = 1000,
        shape = "circle",
        what = "lsm_l_shdi"
      )
      msd1$id <- i
      return( msd1)
    }, error = function(e) {
      message(sprintf("Error at point %d: %s", i, e$message))
      NULL
    })
  })  
  
  
  # Combine results 
  metrics1_df <- do.call(rbind, metrics1)
  # write.csv(metrics1_df, "Data/Saved/edrr_shdi.csv")
  
  shdi <- metrics1_df[, c(4,6)]
  
  edrr_with_landcover <- bind_cols(edrr_with_landcover, shdi[match(1:nrow(edrr_with_landcover), shdi$id), ])
  edrr_with_landcover <- edrr_with_landcover[, -c(27,30)]
  edrr_with_landcover <- edrr_with_landcover %>% rename(shdi = value)
  
  write.csv(edrr_with_landcover, "Data/Saved/edrr_landscape_metrices.csv", row.names = FALSE)
  
  
  
  edrr_with_landcover <- read.csv("Data/Saved/edrr_landscape_metrices.csv")
  
  
# ---- Distribution of landcover classes among trap sites ----  
  
 colnames(edrr_with_landcover)
 
  
  edrr_landcover <- edrr_with_landcover %>% 
                 dplyr:: rename(
                  "Open_Water" = c11_pland,
                  "Developed_Open_Space" = c21_pland , 
                  "Developed_Low_Intensity" = c22_pland,
                  "Developed_Medium_Intensity" = c23_pland,
                  "Developed_High_Intensity" = c24_pland,
                  "Barren_Land" = c31_pland,
                  "Deciduous_Forest" = c41_pland,
                  "Evergreen_Forest" = c42_pland,
                  "Mixed_Forest" = c43_pland,
                  "Shrub_Scrub" = c52_pland,
                  "Grassland_Herbaceous" = c71_pland,
                  "Pasture_Hay" = c81_pland,
                  "Cultivated_Crops" = c82_pland,
                  "Woody_Wetlands" = c90_pland ,
                  "Emergent_Herbaceous_Wetlands" = c95_pland
                 )
                   
                 
  
  colnames(edrr_landcover)
  
   site_landcover <- edrr_landcover %>% 
                   group_by(TrapSite) %>%
                  dplyr::summarise(across(c(Developed_Open_Space, Barren_Land,Deciduous_Forest, Evergreen_Forest, Mixed_Forest,Shrub_Scrub, 
                                            Grassland_Herbaceous,Pasture_Hay, Woody_Wetlands,Developed_Low_Intensity,Developed_Medium_Intensity, 
                                          Developed_High_Intensity, Cultivated_Crops, Emergent_Herbaceous_Wetlands, Open_Water, LC_class), first), .groups="drop")

   
   
   
   head(site_landcover)  
   str(site_landcover)
   
   # site_landcover$total_area <- rowSums(site_landcover[,-c(1, 17)], na.rm = TRUE)
   
   
  
   
  
   site_long <- site_landcover %>%
     pivot_longer(
       cols = -c(TrapSite,LC_class),
       names_to = "LandCover",
       values_to = "Percent"
     ) %>%
     mutate(
       LandCover = factor(
         LandCover,
         levels = c(
           "Open_Water",
           "Developed_Open_Space", 
           "Developed_Low_Intensity",
           "Developed_Medium_Intensity",
           "Developed_High_Intensity",
           "Barren_Land",
           "Deciduous_Forest",
           "Evergreen_Forest",
           "Mixed_Forest",
           "Shrub_Scrub",
           "Grassland_Herbaceous",
           "Pasture_Hay",
           "Cultivated_Crops",
           "Woody_Wetlands",
           "Emergent_Herbaceous_Wetlands"
         )
       )
     )
   
   plot_a <- ggplot(site_long, aes(x = LandCover, y = Percent)) +
     geom_boxplot(fill = "lightgreen") +
     theme_bw() +
     labs(
       x = "Landcover Class",
       y = "Percentage"
     ) +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
   ggsave("Plots/landcover_distribution.jpg", plot = plot_a, width = 8, height = 6, dpi = 600)
 
   
   
   unique(site_landcover$LC_class)
   
   site_landcover$LC_class <- factor(
     site_landcover$LC_class,
     levels = c(
       "Open Water",
       "Developed, Open Space", 
       "Developed, Low Intensity",
       "Developed, Medium Intensity",
       "Developed, High Intensity",
       "Barren Land",
       "Deciduous Forest",
       "Evergreen Forest",
       "Mixed Forest",
       "Shrub/Scrub",
       "Grassland/Herbaceous",
       "Pasture/Hay",
       "Cultivated Crops",
       "Woody Wetlands",
       "Emergent Herbaceous Wetlands"
     )
   )
   


   
   
   
   plot_b <-  ggplot(site_landcover, aes(x = LC_class)) +
     geom_bar(fill = "steelblue") +
     labs(
       x = "Landcover Class",
       y = "Number of Trap Sites"
     ) +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
   ggsave("Plots/landcover_distribution_1.jpg", plot = plot_b, width = 8, height = 6, dpi = 600)
  
# ---- Socio-economic variables ----
  
edrr<-read.csv(file="Data/Saved/edrr_landscape_metrices.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings="", strip.white = T) 
  # Census API key
  #census_api_key("YOUR_API_KEY_HERE", install = TRUE)
  
  readRenviron("~/.Renviron")
  
  
  # Checking if it works
  # get_acs(geography = "state", variables = "B01003_001", year = 2020)
  
  
  edrr_sf <- st_as_sf(edrr, coords = c("Longitude", "Latitude"), crs = 4326) |>
    st_transform(5070) 
  
  
  states <- unique(fips_codes$state)[1:51]  
# ---- Population ----  
  # block groups with total population for all states
  us_block_groups <- map_df(states, function(state_abbr) {
    get_acs(geography = "block group",
            variables = "B01003_001",  # Total population
            year = 2016,
            survey = "acs5",
            state = state_abbr,
            geometry = TRUE,
            output = "wide")
  }) |>
    st_transform(5070)
  
  # The year = 2016 argument specifically requests the data collected between January 2012 and December 2016. 
  #  The 5-year ACS estimates because these are available for small geographic areas as well 
  
  buffers <- st_buffer(edrr_sf, dist = 1000)  
  
  us_block_groups$area_km2 <- st_area(us_block_groups) |> set_units("km^2") |> drop_units()
  us_block_groups$pop_density <- us_block_groups$B01003_001E / us_block_groups$area_km2
  
  
  # plot(us_block_groups$geometry)
  # Convert to a SpatVector
  pop_vect <- vect(us_block_groups)
  
  # Creating an empty raster with desired resolution 
  template <- rast(ext(pop_vect), resolution = 100, crs = crs(pop_vect))
  
  # Rasterize population density
  pop_density_raster <- rasterize(pop_vect, template, field = "pop_density")
  
  # match CRS
  buffers <- st_transform(buffers, crs(pop_density_raster))
  
  # mean population density from raster within each buffer
  buffers$pop_density <- exact_extract(pop_density_raster, buffers, 'mean')
  
  # estimated population in buffer
  buffers$pop_total <- buffers$pop_density * pi * (1^2)  # ~3.14 km²
 
  
  edrr$pop_density<-buffers$pop_density
  edrr$pop_total<- buffers$pop_total

 
  

  
# ---- Income per capita ----
  
  us_block_groups1 <- map_df(states, function(state_abbr) {
    get_acs(geography = "block group",
            variables = "B19301_001",  # per capita income
            year = 2016,
            survey = "acs5",
            state = state_abbr,
            geometry = TRUE,
            output = "wide")
  }) |>
    st_transform(5070)
  
  buffers <- st_buffer(edrr_sf, dist = 1000)
  
  income_vect <- vect(us_block_groups1)
  
  template1 <- rast(ext(income_vect), resolution = 100, crs = crs(income_vect))
  income_raster <- rasterize(income_vect, template1, field = "B19301_001E")
  
  buffers$per_capita_income <- exact_extract(income_raster, buffers, 'mean')

  
  # hist(buffers$per_capita_income)
  buffers$per_capita_income <- as.numeric(buffers$per_capita_income)
  
  edrr$per_capita_income <- buffers$per_capita_income
  
  hist(edrr$per_capita_income)
  
  # write.csv(edrr, "Data/Saved/edrr_landscapemetrics_pop_income.csv", row.names = FALSE)
  # 
  # writeRaster(income_raster, "Data/Saved/income_raster.tif")
  # 

  
 
# ---- PRISM 800m ----
  edrr<-read.csv(file="Data/Saved/edrr_landscapemetrics_pop_income.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings="", strip.white = T) 
 
  prism_set_dl_dir("Data/Processed/prism_monthlys")
  
  # get_prism_monthlys(
  #   type = "ppt",
  #   years = 2009:2019,
  #   mon = c(1:12),
  #   resolution = "800m",
  #   keepZip = FALSE
  # )
  # 
  # get_prism_monthlys(
  #   type = "tmin",
  #   years = 2009:2019,
  #   mon = c(1:12),
  #   resolution = "800m",
  #   keepZip = FALSE
  # )
 
  
# ---- Minimum winter temperature ----
  library(terra)
  library(stringr)
  library(lubridate)
  library(dplyr)
  
 
 
  prism_dir <- "Data/Processed/prism_monthlys"
  
  tmin_files <- list.files(
    prism_dir,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "prism_tmin_us_30s_[0-9]{6}\\.bil$"
  )
  # list.files(prism_dir, recursive = TRUE)
  # Extract YYYYMM and convert to dates
  file_index <- tibble(
    file = tmin_files,
    ym   = str_extract(basename(tmin_files), "[0-9]{6}"),
    date = ymd(paste0(ym, "01"))
  )
  
 
  edrr$TrapYear <- as.numeric(edrr$TrapYear)
  edrr$min_wintemp <- NA_real_
  
  pts <- vect(edrr[, c("Longitude", "Latitude")],
              geom = c("Longitude", "Latitude"),
              crs = "EPSG:4326")
  
  pb <- txtProgressBar(min = 1, max = nrow(edrr), style = 3)
  
 
  # Loop by observation 
 
  for (i in seq_len(nrow(edrr))) {
    
    yr <- edrr$TrapYear[i]
    
    winter_dates <- ymd(c(
      paste0(yr - 1, "1101"),
      paste0(yr - 1, "1201"),
      paste0(yr,     "0101"),
      paste0(yr,     "0201")
    ))
    
    files_4 <- file_index$file[file_index$date %in% winter_dates]
    
    # Safety check
    if (length(files_4) < 4) {
      next
    }
    
    r4 <- terra::rast(files_4)
    pt <- terra::vect(cbind(edrr$Longitude[i], 
         edrr$Latitude[i]),  crs = "EPSG:4326")
    if (!terra::same.crs(pt, r4)) 
        {pt <- terra::project(pt, terra::crs(r4))}
    vals <- terra::extract(r4, pt)
    edrr$min_wintemp[i] <- mean(as.numeric(vals[ , -1]), na.rm = TRUE)
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  
  
# ---- Annual Precipitation ---- 
  
  
  ppt_files <- list.files(
    prism_dir,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "prism_ppt_us_30s_[0-9]{6}\\.bil$"   
  )
  
  # Extract YYYYMM and convert to dates
  file_index <- tibble(
    file = ppt_files,
    ym   = str_extract(basename(ppt_files), "[0-9]{6}"),
    date = ymd(paste0(ym, "01"))
  )
  
  edrr$TrapYear <- as.numeric(edrr$TrapYear)
  edrr$ann_cumppt <- NA_real_   
  
  pts <- vect(edrr[, c("Longitude", "Latitude")],
              geom = c("Longitude", "Latitude"),
              crs = "EPSG:4326")
  
  pb <- txtProgressBar(min = 1, max = nrow(edrr), style = 3)
  
  # Loop by observation 
  for (i in seq_len(nrow(edrr))) {
    
    yr <- edrr$TrapYear[i]
    
    # Previous-year March → current-year February
    ppt_dates <- ymd(c(
      paste0(yr - 1, sprintf("%02d", 3:12), "01"),
      paste0(yr,     sprintf("%02d", 1:2),  "01")
    ))
    
    files_12 <- file_index$file[file_index$date %in% ppt_dates]
    
    # Safety check
    if (length(files_12) < 12) next
    
    r12 <- terra::rast(files_12)
    pt <- terra::vect(
      cbind(edrr$Longitude[i], edrr$Latitude[i]),
      crs = "EPSG:4326"
    )
    
    if (!terra::same.crs(pt, r12)) {
      pt <- terra::project(pt, terra::crs(r12))
    }
    
    vals <- terra::extract(r12, pt)
    
    # Annual cumulative precipitation
    edrr$ann_cumppt[i] <- sum(as.numeric(vals[, -1]), na.rm = TRUE)
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  write.csv(edrr, "Data/saved/edrr_predictors.csv", row.names = FALSE)
  
# ---- Extracting Centroids as Cells ----
  edrr<-read.csv(file="Data/Saved/edrr_predictors.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings="", strip.white = T) 
  
   prism_set_dl_dir("Data/Raw/prism_4km")
   
   get_prism_dailys(
     type = "tmin", 
     minDate = "2010-01-01", 
     maxDate = "2010-01-02", 
     resolution = "4km" 
   )
  my_dates <- seq(as.Date("2010-01-01"), as.Date("2010-01-02"), 1 )
  i <- 1
  prism_raster <- pd_stack(prism_archive_subset("tmin", 
                                                minDate = my_dates[i], maxDate = my_dates[i], 
                                                temp_period = "daily", resolution = "4km"))[[1]]
  
  
  
  edrr_sf <- st_as_sf(edrr, coords = c("Longitude", "Latitude"), crs = crs(prism_raster))
  crs(edrr_sf)
  crs(prism_raster)
  
  
  ## Cell numbers for each point
  cell_num <- cellFromXY(prism_raster, st_coordinates(edrr_sf))
  
  ## Centroids of those cells
  centroids <- xyFromCell(prism_raster, cell_num)
  
  # Add centroids as new columns to your data
  edrr$centroid_lon <- centroids[,1]
  edrr$centroid_lat <- centroids[,2]
  
  
  # For all rows:
  edrr$Cell <- apply(centroids, 1, function(x) paste0("POINT(", x[1], " ", x[2], ")"))
  
  edrr1 <- edrr %>%
    group_by(Cell) %>%
    arrange(Cell, .by_group = TRUE)
  
  edrr$Cell<- as.factor(edrr$Cell)
  length(unique(edrr$Cell)) #755
  length(unique(edrr$TrapSite)) #2257
  
  str(edrr$Cell)
  unique_Cell_ID<-edrr %>%
    distinct(Cell) %>%
    mutate(Cell_ID=paste0("Cell", row_number()))
  
  
  edrr<- edrr %>%
    left_join(unique_Cell_ID, by="Cell")
  
  edrr_Cell1<-edrr[which(edrr$Cell_ID=="Cell1"),]
  table(edrr_Cell1$TrapID)
  edrr<-edrr[, -c(36:38)]
  write.csv(edrr, "Data/saved/edrr_predictors_cell.csv", row.names = FALSE)
  
  
  
  
  
  
  
  
  
# ---- EDRR DATA WITH PREDICTOR VARIABLES ----
  edrr<-read.csv(file="Data/Saved/edrr_predictors_cell.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
  edrr <- edrr %>%
    mutate(Species= dplyr::recode(Species, "Pityogenes kneckteli" = "Pityogenes knechteli"))
 
  

  
  
  
# ---- Extracting trait information from Scolytine Database and Beetles of North America dataset into the EDRR data ----
 
  
  
  
  sco_dist<-read.csv(file="Data/Processed/Scolytinae_Distribution.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings="", strip.white = T) # Scolytinae distribution from the database
  sco_dist<-sco_dist[, 1:6]
  
  
  # Split the string into words
  words <- str_split(sco_dist$Species.name, "\\s+", simplify = TRUE)
  
  # create the ScientificName column
  sco_dist$Species <- ifelse(
    str_detect(words[,3], "^[a-z]"),
    paste(words[,1], words[,2], words[,3]),   # third word is subspecies (lowercase)
    paste(words[,1], words[,2])               # third word is author (uppercase) or missing
  )
  
  sco_dist <-  sco_dist %>%
    mutate(Species= dplyr::recode(Species,
                           "Ips hunteri," = "Ips hunteri"
                          
                            ))
  
  sco_dist_us<-sco_dist[which(sco_dist$Area=="United States"), ] 
  unique(sco_dist_us$Species)
  
  # write.csv(sco_dist_us, "Data/Saved/sco_dist_us.csv")
  
  
  
  sco_dist_us_aggregate<- sco_dist_us %>%
    group_by(Species) %>%
    summarise(
      Species.status=paste(unique(`Species.status`), collapse=",")
      
    )
  
  #write.csv(sco_dist_us_aggregate, "Data/Saved/sco_dist_us_aggregate.csv")  # Species list from Scolytine database with Status (Native/Introduced) information
  
  # EDRR species
  
  edrr_species <- as.data.frame(unique(edrr$Species))
  edrr_species <- edrr_species %>% rename(Species = `unique(edrr$Species)`)
  edrr_species <- edrr_species[order(edrr_species$Species),]
  edrr_species <- as.data.frame(edrr_species)
  
  # Checking if all the species from EDRR data occur in the Scolytine database
  
  common_species <- intersect(edrr_species$edrr_species, sco_dist_us_aggregate$Species)
  common_species <- as.data.frame(common_species) # 260
  
  
  # Species in EDRR data that are not in the Scolytine database
  
  # Species_not_in_database <- setdiff(edrr_species$edrr_species, sco_dist_us_aggregate$Species)
  # Species_in_EDRR_not_in_database <- as.data.frame(Species_not_in_database) # 60

  
  # print(Species_in_EDRR_not_in_database)
  # 
  # write.csv(common_species, "Data/Saved/common_species_edrr_scolytinedatabase.csv")
  # 
  # write.csv(Species_in_EDRR_not_in_database, "Data/Saved/Species_in_edrr_not_in_database.csv")
  
  # Species_in_EDRR_not_in_Scolytine_database
  
  # 1            Ambrosiodmus lecontei Native
  # 2                 Ambrosiodmus sp.
  # 3         Ambrosiophilus nodulosis
  # 4                   Carphoborus sp
  # 5                        Chramesus
  # 6                      Coccotrypes
  # 7                 Conophthorus sp.
  # 8               Corthylus spinifer
  # 9                        Cryphalus
  # 10                      Crypturgus
  # 11         Cyclorhipidion fukiense
  # 12        Cyclorhipidion inarmatum 
  # 13         Dendroctonus brevicomis 
  # 14       Dendroctonus pseudotsugae 
  # 15                      Dryocoetes
  # 16         Euwallacea denticulatus
  # 17                        Hylastes
  # 18                       Hylesinus
  # 19                       Hylocurus
  # 20               Hylocurus torosus
  # 21    Hylurgops rugipennis pinifex
  # 22 Hylurgops rugipennis rugipennis
  # 23         Hylurgops subcostulatus
  # 24                    Hypothenemus
  # 25                    Ips borealis 
  # 26                Ips calligraphus 
  # 27               Ips plastographus 
  # 28                        Micracis
  # 29                    Micracisella
  # 30                     Phloeosinus
  # 31   Phloeosinus scopulorum s.str.
  # 32                    Phloeotribus
  # 33                      Pityoborus
  # 34                      Pityogenes
  # 35            Pityogenes kneckteli
  # 36            Pityogenes plagiatus 
  # 37                    Pityokteines
  # 38            Pityokteines ornatus
  # 39                   Pityophthorus
  # 40               Pityophthorus spp
  # 41                  Premnobius spp
  # 42              Pseudips mexicanus
  # 43          Pseudohylesinus dispar Native
  # 44       Pseudohylesinus nebulosus Native
  # 45             Pseudopityophthorus
  # 46                 Pseudothysanoes
  # 47             Scolytodes schwarzi
  # 48            Scolytogenes jalapae
  # 49                Scolytogenes sp.
  # 50                        Scolytus
  # 51                 Scolytus opacus
  # 52                Theoborus ricini
  # 53                       Thysanoes
  # 54                     Trischidias
  # 55                    Trypodendron
  # 56           Trypodendron borealis
  # 57          Wallacellus denticulus
  # 58              Xyleborus intrusus
  # 59              Xyleborus seriatus Introduced
  # 60                   Xyleborus sp.
  # > 
  
  
  
  
  ## Some species from EDRR dataset are missing on the Scolytine database. Therefore, adding those species (with confirmed status) and status info into the sco_dist_us_aggregate manually before extracting status info
  
  
  additional_species <- data.frame(
    Species = c( "Ambrosiodmus lecontei",
                 "Cyclorhipidion inarmatum",
                 "Hylurgops subcostulatus",
                 "Xyleborus seriatus",
                 "Dendroctonus brevicomis",
                 "Ips borealis",
                 "Ips calligraphus",
                 "Ips plastographus",
                 "Pityogenes plagiatus",
                 "Pseudohylesinus dispar",
                 "Pseudohylesinus nebulosus",
                 "Dendroctonus pseudotsugae",
                 "Xyleborus intrusus",
                 "Pityogenes knechteli",
                 "Pityokteines ornatus",
                 "Hylurgops rugipennis rugipennis",
                 "Hylurgops rugipennis pinifex",
                 "Scolytus opacus"), 
 Species.status = c("Native",
                    "Introduced",
                    "Native",
                    "Introduced",
                    "Native",
                    "Native", 
                    "Native",
                    "Native",
                    "Native",
                    "Native",
                    "Native",
                    "Native",
                    "Native",
                    "Native",
                    "Native",
                    "Native",
                    "Native",
                    "Native"),
 stringsAsFactors = FALSE
 )
  
  sco_species <- rbind(sco_dist_us_aggregate, additional_species)
  
  sco_species$Species.status[sco_species$Species %in% c("Cryptocarenus heveae", "Monarthrum fasciatum", "Pityophthorus juglandis", "Xyloterinus politus")] <- "Native"
 
  write.csv(sco_species, "Data/Saved/species_nativity.csv", row.names = FALSE)
  sco_species<-read.csv(file="Data/Saved/species_nativity.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings="", strip.white = T) 
  table(sco_species$Species.status)
  
  
  edrr_sp <- left_join(edrr, sco_species, by="Species")
  
  edrr_sp$Species.status <- as.factor(edrr_sp$Species.status)
  levels(edrr_sp$Species.status)
  

# ---- Traits ---- 
  
  sco_trt<-read.csv(file="Data/Processed/Scolytinae_General.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings="", strip.white = T) # Scolytinae General dataset from the database with traits
  sco_trt<-sco_trt[, 1:5]
  
  
  
  words1 <- str_split(sco_trt$Species.name, "\\s+", simplify = TRUE)
  
  # create the ScientificName column
  sco_trt$Species <- ifelse(
    str_detect(words1[,3], "^[a-z]"),
    paste(words1[,1], words1[,2], words1[,3]),   
    paste(words1[,1], words1[,2])               
  )
  
  
  
  
  
  
  sco_trt_aggregate<- sco_trt %>%
    group_by(Species) %>%
    summarise(
      Reproduction=paste(unique(Reproduction), collapse=","),
      Feeding_Habits=paste(unique(`Feeding.Habits`), collapse=","),
      
    )
  
  edrr_traits<- left_join(edrr_sp, sco_trt_aggregate, by="Species")
  
  #write.csv(edrr_sco_status_feedhabit, "Data/Saved/edrr_sco_status_feedhabit.csv")
  
  edrr_traits_aggregate<- edrr_traits %>%
    group_by(Species) %>%
    summarise(
      
      Species.status=paste(unique(Species.status), collapse=","),
      Reproduction=paste(unique(Reproduction), collapse = ","),
      Feeding_Habits=paste(unique(Feeding_Habits), collapse = ","),
      
    )
  

  table(edrr_traits_aggregate$Species.status)
  
  statusNA <- edrr_traits_aggregate[which(edrr_traits_aggregate$Species.status=="NA"), ]

  
  # repeat_1<- as.data.frame(table(sco_species$Species))
  # count<-as.data.frame(table(repeat_1$Freq))
# ---- Dataset "Beetles of North America ----
  
  bna<-read.csv(file="Data/Processed/bna.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings="", strip.white = T) # Dataset "Beetles of North America" or bna
  hostfamilies<-as.data.frame(unique(bna$Host.family))
  a<-as.data.frame(sort(hostfamilies$`unique(bna$Host.family)`))
  bna$Species=paste(bna$Genus,bna$Species) # Creating a column with combined Genus and Species names

  # #
  bna_aggregate<- bna %>%
    group_by(Species) %>%
    summarise(
      Statusinfo=paste(unique(Nativeness.in.North.America), collapse=","),
      Feeding_Guild=paste(unique(Feeding.guild), collapse=","),
      
    )
  bna_aggregate <- bna_aggregate[,-2 ]
  
  
    table(bna_aggregate$Feeding_Guild)
   
   additional_guildinfo <- data.frame(
     Species = c( "Anisandrus maiche",
                  "Coptoborus pseudotenuis",
                  "Cyclorhipidion distinguendum",
                  "Cyclorhipidion inarmatum",
                  "Cyclorhipidion pelliculosum",
                  "Euwallacea similis",
                  "Xyleborinus artestriatus",
                  "Xyleborinus octiesdentatus",
                  "Xyleborus seriatus",
                  "Xylosandrus amputatus", 
                  "Trypodendron scabricollis", 
                  "Xyleborus viduus",
                  "Chramesus chapuisi",
                  "Ips borealis borealis", 
                  "Ips borealis lanieri", 
                  "Ips pilifrons pilifrons",
                  "Ips pilifrons utahensis", 
                  "Ips tridens engelmanni",
                  "Phloeosinus scopulorum neomexicanus", 
                  "Phloeosinus taxodii taxodii", 
                  "Phloeosinus vandykei",
                  "Ambrosiodmus lecontei",
                  "Cryphalus ruficollis",
                  "Hylurgops rugipennis rugipennis",
                  "Hylurgops rugipennis pinifex",
                  "Hylurgops subcostulatus",
                  "Pseudohylesinus nebulosus nebulosus",
                  "Trypodendron rufitarsis"), 
     Feeding_Guild = c("ambrosia beetle",
                        "ambrosia beetle",
                        "ambrosia beetle",
                        "ambrosia beetle",
                        "ambrosia beetle",
                        "ambrosia beetle", 
                        "ambrosia beetle",
                        "ambrosia beetle",
                        "ambrosia beetle",
                        "ambrosia beetle", 
                        "ambrosia beetle", 
                        "ambrosia beetle",
                        "bark beetle", 
                        "bark beetle", 
                        "bark beetle", 
                        "bark beetle", 
                        "bark beetle", 
                        "bark beetle", 
                        "bark beetle", 
                        "bark beetle", 
                        "bark beetle",
                       "ambrosia beetle",
                       "bark beetle",
                       "bark beetle",
                       "bark beetle",
                       "bark beetle",
                       "bark beetle",
                       "ambrosia beetle"),
     stringsAsFactors = FALSE
   )
  
 
   
   bna_aggregate <- rbind(bna_aggregate, additional_guildinfo)
   
  
  
  edrr_traits<- left_join(edrr_traits, bna_aggregate,  by="Species")
  
  
  
  
  edrr_traits_agg<- edrr_traits %>%
    group_by(Species) %>%
    summarise(
      
      Species.status=paste(unique(Species.status), collapse=","),
      Feeding_Guild=paste(unique(Feeding_Guild), collapse=","),
      Reproduction=paste(unique(Reproduction), collapse = ","),
      Feeding_Habits=paste(unique(Feeding_Habits), collapse = ",")
      
    )
  
  table(edrr_traits_agg$Species.status)
  
  
  
  
 
   table(edrr_traits_agg$Feeding_Guild)
   
   statusNA_guild <- edrr_traits_agg[which(edrr_traits_agg$Feeding_Guild=="NA"), ]
   
   diff <- setdiff(statusNA_guild$Species,  statusNA$Species)
  diffdf <- as.data.frame(diff) # 3
  print(diff)
  # [1] "Hylocurus biconcavus"   "Hylocurus flaglerensis" "Hylocurus schwarzi"   
  
  
  diff1 <- setdiff( statusNA$Species,statusNA_guild$Species)
   print(diff1)
   # [1] "Pseudips mexicanus" "Theoborus ricini"    Origin is uncertain
   
   
  host_breadth<- bna %>%
    group_by(Species) %>%
    summarise(unique_host_count = n_distinct(Host.family))
  
  #Checking the results above
  # bna_xc<-bna[which(bna$Species=="Xylosandrus compactus"),]
  # length(unique(bna_xc$Host.family)) # 24
  # bna_xg<-bna[bna$Species=="Xylosandrus germanus",]
  # length(unique(bna_xg$Host.family)) # 18

  
  edrr_traits<-left_join(edrr_traits, host_breadth, by="Species")
  
  
  edrr_traits$Host_range <- ifelse(edrr_traits$unique_host_count == 1, "onegenus", "multigenus")
  
  
  
  
 guild<- edrr_traits %>%
    group_by(Species) %>%
    summarise(
      Feeding_Guild=paste(unique(Feeding_Guild), collapse=","),
      
    )
  
  
  
  write.csv(edrr_traits, "Data/Saved/edrr_predictors_and_traits.csv", row.names = FALSE)
  
  
  
  
  
  
  
  
  
  
  
# ---- EDRR DATA WITH PREDICTOR VARIABLES AND BEETLE TRAITS ----
  edrr<-read.csv(file="Data/Saved/edrr_predictors_and_traits.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
 

# ---- Nearest Ports ----
  library(sf)
  
  edrr_sf <- st_as_sf(
    edrr,
    coords = c("Longitude", "Latitude"),
    crs = 4326
  )
  
  
  ports_sf <- st_read("Data/Raw/Principal_Ports.shp")
  # plot(ports_sf$geometry)
  
  edrr_sf <- st_transform(edrr_sf, 5070)   # Albers Equal Area CONUS
  ports_sf <- st_transform(ports_sf, 5070)
  
  ports_cg <- ports_sf %>% filter(TYPE %in% c("Coastal", "Great Lakes"))
  
  
  nearest_port <- st_nearest_feature(edrr_sf, ports_cg)
  
  edrr_sf$nearest_port <- ports_cg$PORTNAME[nearest_port]
  
  edrr_sf$nearest_port_type <- ports_cg$TYPE[nearest_port]
  
  edrr_sf$dist_port_km <- as.numeric(st_distance(edrr_sf,ports_cg[nearest_port, ],by_element = TRUE)) / 1000
  
  
###### All ports #### 

  
  nearest_port_all <- st_nearest_feature(edrr_sf, ports_sf)
  
  edrr_sf$nearest_port_all <- ports_sf$PORTNAME[nearest_port_all]
  
  edrr_sf$nearest_port_type_all <- ports_sf$TYPE[nearest_port_all]
  
  edrr_sf$dist_port_km_all <- as.numeric(st_distance(edrr_sf,ports_sf[nearest_port_all, ],by_element = TRUE)) / 1000
  
  
  # albany <- ports_sf[ports_sf$PORTNAME == "Albany Port District, NY", ]
  # plot(albany$geometry)
  # 
  # st_within(edrr_sf[88, ],  albany,  sparse = FALSE)
######
  edrr_new <- edrr_sf %>%
    st_transform(4326) %>% 
    mutate(
      Longitude = st_coordinates(.)[, 1],
      Latitude  = st_coordinates(.)[, 2]
    ) %>%
    st_drop_geometry()
  # edrr_new <- edrr_new[, -c(46:48)]
  
  write.csv(edrr_new, "Data/Saved/edrr_predictors_and_traits.csv", row.names = FALSE)
  
  # 
  # hist(edrr$dist_port_km)
  # edrr$PortDist <- as.numeric(scale(log1p(edrr$dist_port_km)))
  
# ---- RENAMING AND SCALING PREDICTOR VARIABLES ----
  edrr<-read.csv(file="Data/Saved/edrr_predictors_and_traits.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)

  colnames(edrr)
  
  
    
    edrr <- edrr %>%
      rename(
        Annual_Precipitation = ann_cumppt,
        Minimum_Temperature = min_wintemp,
       LC_Deciduous = c41_pland, 
       LC_Evergreen = c42_pland, 
        Shannon_Div_Landscape= shdi,
       Income_PerCapita = per_capita_income,
       PortDistance = dist_port_km, 
       PortDistance_all = dist_port_km_all
       
      )

    
  
    predictor_variables <- edrr %>% dplyr::select(Annual_Precipitation, Minimum_Temperature, LC_Evergreen,
                                           LC_Deciduous, Shannon_Div_Landscape, pop_density, Income_PerCapita,
                                           PortDistance, PortDistance_all)
    
    head(predictor_variables)
  str(predictor_variables)
   predictor_summary <- edrr %>%
         group_by(TrapSite) %>%
     summarise(across(c(Annual_Precipitation, Minimum_Temperature, LC_Evergreen,
                        LC_Deciduous, Shannon_Div_Landscape, pop_density, Income_PerCapita,
                        PortDistance), first), .groups = "drop")
   
   library(dplyr)
  pred_Var <- predictor_summary %>% dplyr::select(Annual_Precipitation, Minimum_Temperature, LC_Evergreen,
                                                 LC_Deciduous, Shannon_Div_Landscape, 
                                                 pop_density, Income_PerCapita, PortDistance)
   
                                                 
  library(tibble)                                             
  pred_stat <- tibble(
    Variable = names(pred_Var),
    `Mean ± SD` = paste0(
      round(sapply(pred_Var, mean, na.rm = TRUE), 2),
      " ± ",
      round(sapply(pred_Var, sd, na.rm = TRUE), 2)),
      Range = paste0(
        round(sapply(pred_Var, min, na.rm = TRUE), 2),
        " to ",
        round(sapply(pred_Var, max, na.rm = TRUE), 2)
    )
  )
  


   pred_stat
   getwd()
   write.csv(pred_stat, "Data/Saved/predictor_stats.csv", row.names = FALSE)
   ps<-read.csv("Data/Saved/predictor_stats.csv")
   
    par(mfrow = c(3, 3)) 
    for (i in 1:ncol(predictor_variables)) {
      hist(predictor_variables[, i], main = names(predictor_variables)[i], 
           xlab = names(predictor_variables)[i], col = "skyblue")
    }   
    
  
  ## The landcover classes, population density, per capita income, and distance to port are highly skewed
    ## Log transforming 
    
      edrr$Deciduous <- log1p(edrr$LC_Deciduous)
      edrr$Evergreen <- log1p(edrr$LC_Evergreen)
      edrr$POPDEN <- log1p(edrr$pop_density)
      edrr$Income_PerCapita <- log1p(edrr$Income_PerCapita)
      edrr$PortDist <- log1p(edrr$PortDistance) 
      edrr$PortDist_all <- log1p(edrr$PortDistance_all) 
    
    
      
      ### predictor variables after log transforming 
    pred_var_log <- edrr %>% select(Annual_Precipitation, Minimum_Temperature, Evergreen,
                                   Deciduous, Shannon_Div_Landscape, POPDEN, Income_PerCapita,
                                   PortDist, PortDist_all)
    
    par(mfrow = c(3, 3)) 
    for (i in 1:ncol(pred_var_log)) {
      hist(pred_var_log[, i], main = names(pred_var_log)[i], 
           xlab = names(pred_var_log)[i], col = "skyblue")
    }   

  
    
    summary(predictor_variables)
    summary(pred_var_log)
      
 ########## Rescaling the predictors
    
    pred_var_log <- edrr%>% select(Annual_Precipitation, Minimum_Temperature, Evergreen,
                                   Deciduous, Shannon_Div_Landscape, POPDEN, Income_PerCapita,
                                   PortDist, PortDist_all)
    Othervars<-edrr%>% select(-c(Annual_Precipitation, Minimum_Temperature, Evergreen,
                              Deciduous, Shannon_Div_Landscape, POPDEN, Income_PerCapita,
                              PortDist, PortDist_all))
    scaled_predvars<-as.data.frame(scale(pred_var_log)) 
    
    edrr_final<-cbind(Othervars, scaled_predvars)
    
    par(mfrow = c(3, 3)) 
    for (i in 1:ncol(scaled_predvars)) {
      hist(scaled_predvars[, i], main = names(scaled_predvars)[i], 
           xlab = names(scaled_predvars)[i], col = "skyblue")
    }   
    
    write.csv(edrr_final, "Data/Saved/edrr_final_for_analysis.csv", row.names = FALSE)
    
   
    
    
    
    
    
    
    
    
    
    
    
################################################################
################################################################
################################################################
################################################################
# ---- FINAL EDRR DATA FOR ANALYSES ----   
    edrr<-read.csv(file="Data/Saved/edrr_final_for_analysis.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
    
   
    na_counts <- as.data.frame(colSums(is.na(edrr)))
    edrr <- edrr[!is.na(edrr$Income_PerCapita), ]
    
    edrr$TrapYear<- as.factor(edrr$TrapYear)
    edrr$Cell_ID<- as.factor(edrr$Cell_ID)
    edrr$LureName<- as.factor(edrr$LureName)
    edrr$Ecoregion3<- as.factor(edrr$Ecoregion3)
    edrr$TrapSite<- as.factor(edrr$TrapSite)
    
  
    levels(edrr$LureName)
    table(edrr$LureName)
    
    str(edrr)
    count<-as.data.frame(edrr$PestCount)
    count <- count |> arrange(`edrr$PestCount`)
    
    edrr_state_sites <- edrr %>%
      summarise(n_sites = n_distinct(TrapSite), n_year=n_distinct(TrapYear), .by = StateCode)
    edrr_state_sites <- edrr_state_sites %>% arrange(n_sites)
    
    edrr_year_sites <- edrr %>% summarise(n_sites = n_distinct(TrapSite),
                                          .by = TrapYear)
    edrr_year_sites <- edrr_year_sites %>%  arrange(TrapYear)
    mean(edrr_year_sites$n_sites)
    sum(edrr_year_sites$n_sites)
    
    
    edrr_taxa_traits <- edrr %>%
      summarise(Origin=unique(Species.status), Feeding_Guild =unique(Feeding_Guild),Count = sum(PestCount, na.rm = TRUE), .by = Species)
    edrr_taxa_traits <- edrr_taxa_traits %>% arrange(Species)
    species_counts <- edrr_taxa_traits %>%  arrange(desc(Count))
    
    length(unique(edrr$TrapID)) #3708
    length(unique(edrr$TrapSite)) #2248
    length(unique(edrr$StateCode))
    
    write.csv(edrr_taxa_traits, "Data/Saved/edrr_taxa_traits.csv", row.names = FALSE)
    write.csv(species_counts, "Data/Saved/taxa_trait_count.csv", row.names = FALSE)
    
    length(unique(edrr$SiteType))
    

    
###################################################
###################################################
###################################################
################################################### 
    
## ---- UPDATED SUMMARY ----
  
  library(ggh4x)
  edrr_summary <- edrr %>%
    dplyr::filter(Species.status %in% c("Native", "Introduced")) %>%
    dplyr::mutate(Origin = factor(Species.status, levels = c("Native", "Introduced"))) %>%
    dplyr::mutate(
      Guild = dplyr::case_when(
        Feeding_Guild == "ambrosia beetle" ~ "Ambrosia beetle",
        Feeding_Guild == "bark beetle" ~ "Bark beetle",
        Feeding_Guild == "twig beetle" ~ "Twig beetle",
        TRUE ~ "Others"
      ),
      Guild = factor(Guild, levels = c("Bark beetle" ,"Ambrosia beetle", "Twig beetle","Others" ))
    )
  
  
  sum_by <- edrr_summary %>%
    dplyr::group_by(Origin, Guild) %>%

      dplyr::summarise(
        
      Abundance = sum(PestCount, na.rm = TRUE),
      Richness  = dplyr::n_distinct(Species[!is.na(PestCount) & PestCount > 0]),
      .groups = "drop"
    )
  
  total_abundance <- sum(sum_by$Abundance)
  
  sum_long <- sum_by %>%
    pivot_longer(cols = c(Abundance, Richness),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = c("Abundance", "Richness")))
  
  sum_long <- sum_long %>%
    mutate(Origin=dplyr::recode(Origin, "Introduced" = "Non-native"))
  
  
  
  
  
  p_combined <- ggplot(sum_long, aes(x = Origin, y = Value, fill = Guild)) +
    geom_col(width = 0.7, color = "gray20") +
    # Centered labels inside each stacked segment; hide zeros
    geom_text(
      data = subset(sum_long, Value > 0),
      aes(label = Value),                # or label = comma(Value) if using scales
      position = position_stack(vjust = 0.5),
      color = "black",
      size = 3.5
    ) +
    scale_fill_manual(values = c("Ambrosia beetle" = "#D0F0C0", "Bark beetle" = "lightblue",
                                 "Twig beetle" = "#1B9E77","Others" = "grey90")) +
    labs(
      title = "",
      x = NULL,
      y = NULL,
      fill = "Guild"
    ) +
    
    # Facet: two panels (Abundance / Richness)
    ggh4x::facet_wrap2(~ Metric, scales = "free_y") +
    # Give EACH facet its own y-scale formatting
    ggh4x::facetted_pos_scales(
      y = list(
        Metric == "Abundance" ~ scale_y_continuous(labels = label_number(big.mark = ","),  
                                                   expand = expansion(mult = c(0, 0.1))),
        Metric == "Richness"  ~ scale_y_continuous(labels = label_number(accuracy = 1),  
                                                   expand = expansion(mult = c(0, 0.1)))
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "right",
      # ⬇️ "Think box" style: outline around facet strip labels
      strip.background = element_rect(color = "black", fill = "grey75", linewidth = 0.8),
      strip.text = element_text(face = "bold", size = 12),
      # ⬇️ Bigger x-axis labels
      axis.text.x = element_text(size = 14, face = "bold")
    )
  
  
  p_combined
  
  str(sum_by)
  sum_by$Abundance <- as.numeric(sum_by$Abundance)
  sum(sum_by$Abundance, na.rm = TRUE)
  
  sum(edrr$PestCount, na.rm = TRUE)
  
  ggsave("../Plots/summaryplot.jpg", plot = p_combined, width = 10, height = 8, dpi = 600)
  
  

## ---- Grouped Bar ----
  p_combined <- ggplot(
    sum_long,
    aes(x = Origin, y = Value, fill = Guild)
  ) +
    geom_col(
      position = position_dodge(width = 0.85),
      width = 0.75,
      color = "gray20"
    ) +
    
    # Labels above each bar; hide zeros
    geom_text(
      data = subset(sum_long, Value > 0),
      aes(label = Value),
      position = position_dodge(width = 0.85),
      vjust = -0.3,
      hjust = -0.1,
      size = 4.5,
      angle = 45, 
      color = "black"
    ) +
    
    scale_fill_manual(
      values = c(
        "Ambrosia beetle" = "#D0F0C0",
        "Bark beetle"     = "lightblue",
        "Twig beetle"     = "#1B9E77",
        "Others"          = "grey90"
      )
    ) +
    
    labs(
      title = "",
      x = NULL,
      y = NULL,
      fill = "Guild"
    ) +
    
    
    ggh4x::facet_wrap2(~ Metric, scales = "free_y") +
    
    ggh4x::facetted_pos_scales(
      y = list(
        Metric == "Abundance" ~
          scale_y_continuous(
            labels = label_number(big.mark = ","),
            expand = expansion(mult = c(0, 0.15))
          ),
        Metric == "Richness" ~
          scale_y_continuous(
            labels = label_number(accuracy = 1),
            expand = expansion(mult = c(0, 0.15))
          )
      )
    ) +
    
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 18, face = "bold"), 
      legend.text = element_text(size = 16),
      strip.background = element_rect(
        color = "black",
        fill = "grey75",
        linewidth = 0.8
      ),
      strip.text = element_text(face = "bold", size = 18),
      axis.text.x = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 18)
    )
  p_combined
  # current_width  <- dev.size("in")[1]
  # current_height <- dev.size("in")[2]
  ggsave("Plots/summaryplot_grouped.png", plot = p_combined, width = 12, height = 9, units = "in", dpi = 600, bg="white")
########################################################
# ---- Overall Data -----
# ---- Correlation among predictors ----
  colnames(edrr)
  
 pred_vars_edrr <- edrr[, c("Annual_Precipitation","Minimum_Temperature", "Evergreen", "Deciduous",
                       "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist")]
  str(pred_vars_edrr)
  cor_matrix_edrr <- cor(pred_vars_edrr, use = "complete.obs", method = "pearson")
  print(cor_matrix_edrr)
  
  
  high_cor_indices <- which(abs(cor_matrix_edrr) >= 0.7 & abs(cor_matrix_edrr) < 1, arr.ind = TRUE)
  high_cor_pairs <- data.frame(
    Variable1 = rownames(cor_matrix_edrr)[high_cor_indices[, 1]],
    Variable2 = colnames(cor_matrix_edrr)[high_cor_indices[, 2]],
    Correlation = cor_matrix_edrr[high_cor_indices]
  )
  print(high_cor_pairs)

 
 
# ---- Aggregation by Trap ID ----

library(tidyverse)

# Summarize counts
  edrr_long <- edrr %>%
    group_by(TrapID, Species) %>%
    summarise(count = sum(PestCount, na.rm = TRUE), .groups = "drop")
  
  # wide format
  edrr_wide <- edrr_long %>%
    pivot_wider(
      names_from = Species,
      values_from = count,
      values_fill = 0  # Replace NA with 0 
    )
 
  colnames(edrr)
  
  # Extract metadata , keeps first value for each variable
  meta_info <- edrr %>%
    group_by(TrapID) %>%
    summarise(across(c(Cell_ID, SiteLure, StateCode, CountyName,SiteTown,SiteLoc,SiteType, TrapYear, SiteNumber, LureName, FIPS, Longitude, Latitude,
                       TrapSite, landcover, LC_class,  Ecoregion3,  
                       Annual_Precipitation , Minimum_Temperature,  Deciduous, Evergreen,  
                       Shannon_Div_Landscape, POPDEN,  Income_PerCapita, PortDist, PortDist_all), first), .groups = "drop")
                       
                       
  edrr_abundance<- edrr_long %>% left_join(meta_info, by = "TrapID")                    
                        
  edrr_richness<- edrr_wide %>%  left_join(meta_info, by = "TrapID")
  
  
  #Calculating Species Count for each Trap ID
  edrr_richness$SpeciesCount<-rowSums(edrr_richness[, 2:321 ] !=0)
  
  
  
  write.csv(edrr_abundance, "Data/Saved/edrr_abundance.csv", row.names = FALSE)
  write.csv(edrr_richness, "Data/Saved/edrr_richness.csv", row.names = FALSE)

  
  
  mean(edrr_abundance$count)
  var(edrr_abundance$count)
  hist(edrr_abundance$count)
  summary(edrr_abundance$count)
  hist(edrr_richness$SpeciesCount)
 

  
 
 table(edrr_richness$LC_class)
  
  
 
# ---- Checking Correlation on aggregated dataset ----  
  pred_vars_agg <- edrr_abundance[, c("Annual_Precipitation", "Minimum_Temperature", "Evergreen", "Deciduous",
                                    "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist", "PortDist_all")]
  str(pred_vars_agg)
  na_counts <- colSums(is.na(pred_vars_agg))
  print(na_counts)
  cor_matrix_agg <- cor(pred_vars_agg, use = "complete.obs", method = "pearson")
 

  
  high_cor_indices <- which(abs(cor_matrix_agg) >= 0.7 & abs(cor_matrix_agg) < 1, arr.ind = TRUE)
  high_cor_pairs <- data.frame(
    Variable1 = rownames(cor_matrix_agg)[high_cor_indices[, 1]],
    Variable2 = colnames(cor_matrix_agg)[high_cor_indices[, 2]],
    Correlation = cor_matrix_agg[high_cor_indices]
  )
  print(high_cor_pairs)
  
  
  
# ---- Full Model -Abundance ----
  library(glmmTMB)
  library(DHARMa)
  library(car)
   m1<- glmmTMB(count ~ Annual_Precipitation +  Minimum_Temperature + Deciduous + Evergreen + 
                                POPDEN +  Income_PerCapita + Shannon_Div_Landscape + PortDist + LureName + 
                               (1|TrapSite/TrapID) + (1|TrapYear) + (1|Species) + (1|Ecoregion3), 
                               family = nbinom2(),  data=edrr_abundance)
   Anova(m1)
   summary(m1)
   summarym1<- summary(m1)
   coefs_m1 <- as.data.frame(summarym1$coefficients$cond)
   coefs_m1 <- round(coefs_m1, 3)
   coefs_m1$Variable <- rownames(coefs_m1)

###### All ports 
   m1_allPorts<- glmmTMB(count ~ Annual_Precipitation +  Minimum_Temperature + Deciduous + Evergreen + 
                  POPDEN +  Income_PerCapita + Shannon_Div_Landscape + PortDist_all + LureName + 
                  (1|TrapSite/TrapID) + (1|TrapYear) + (1|Species) + (1|Ecoregion3), 
                family = nbinom2(),  data=edrr_abundance)
   Anova(m1_allPorts)
   summary(m1_allPorts)
   
AIC(m1, m1_allPorts)

# ---- Model Diagnostics ----   
#    
   resm1<- simulateResiduals(m1, plot = TRUE)
   
   testOutliers(resm1, type = "bootstrap")
  
   
# ---- Multimodel Inference ----
   
   system.time(update(m1, . ~ .))
   
   library(MuMIn)
   
   options(na.action = "na.fail")

  
   model_set1<- dredge(m1, rank = AIC, trace = TRUE)
   model_averaged<- model.avg(model_set1, subset=delta < 2) 
   modavg<- summary(model_averaged)
   df<-as.data.frame(modavg$coefmat.full)
   CI_df<- as.data.frame(confint(model_averaged, full=T))
   df$CI.min <-CI_df$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
   df$CI.max <-CI_df$`97.5 %`
   df <- round(df, 3)
   df$Variable <- rownames(df)
   write.csv(df, "Data/Saved/Abundance_modavg.csv")
   
   saveRDS(model_set1, file = "Data/Saved/dredge_edrr.rds")

   
   # model_set1 <- readRDS("Data/Saved/dredge_edrr.rds")


# ---- Plot ----
   
   
   
   library(forcats)
   Abundance_overallEDRR<-read.csv(file="Data/Saved/Abundance_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
   Abundance_overallEDRR$Variable <- gsub("cond\\((.*)\\)", "\\1",Abundance_overallEDRR$Variable)
   
  
   Abundance_overallEDRR <- Abundance_overallEDRR %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation", "Minimum_Temperature") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     TRUE ~ "Other"  ))
   Abundance_overallEDRR <- Abundance_overallEDRR[Abundance_overallEDRR$Type != 'Other',]
   
   
   Abundance_overallEDRR <- Abundance_overallEDRR %>% filter(Variable != "(Intercept)")
   
   Abundance_overallEDRR <- Abundance_overallEDRR %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
   Abundance_overallEDRR$Variable <- factor(Abundance_overallEDRR$Variable, 
                                            levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                                           "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
   
  
   Abundance_overallEDRR <- Abundance_overallEDRR %>%
     mutate(Variable = dplyr::recode(Variable,
                              "Annual_Precipitation" = "Annual Precipitation",
                              "Minimum_Temperature" = "Minimum Temperature",
                              "Deciduous" = "Deciduous Forest",
                              "Evergreen" = "Evergreen Forest",
                              "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                              "POPDEN" = "Human Population Density",
                              "Income_PerCapita" = "Per Capita Income",
                              "PortDist" = "Distance from Port"
     ))
    plot1 <- ggplot(Abundance_overallEDRR, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "All Scolytines",

     ) +
      scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 0)) + 
      scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2")) +
      
      annotate("text",x = Inf,y = Inf,label = "A",hjust = 1.4,vjust = 1.3,fontface = "bold",size = 7) +
        
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
        plot.title = element_text(hjust = 0.5, face="bold", size=20),
        strip.text = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
        axis.text.x = element_text(face = "bold", size = 18),
        axis.text.y = element_text(face = "bold", size = 18),
        axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
        axis.title.y = element_text(size = 22, face = "bold", vjust = 2), 
        
        legend.position = c(0.99, 0.02),   # bottom-right INSIDE panel
        legend.justification = c(1, 0),
        legend.background = element_rect(
          fill = "white",
          color = "black",
          linewidth = 0.6
        ),
        legend.key.size = unit(1.2, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold")
       
        
      )+ guides(shape ="none")
    
   
   plot1
  
   ggsave("Plots/Abundance_EDRR.jpg", plot = plot1, width = 10, height = 12, dpi = 600)
# ---- Full MOdel - Richness ----
  
   pred_vars_agg <- edrr_richness[, c("Annual_Precipitation", "Minimum_Temperature", "Evergreen", "Deciduous",
                                       "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist")]
   str(pred_vars_agg)
   na_counts <- colSums(is.na(pred_vars_agg))
   print(na_counts)
   cor_matrix <- cor(pred_vars_agg, use = "complete.obs", method = "pearson")
   print(cor_matrix)
   
   correlation_matrix <- as.data.frame(cor_matrix)
   
   
   
   high_cor_indices <- which(abs(correlation_matrix) >= 0.7 & abs(correlation_matrix) < 1, arr.ind = TRUE)
   high_cor_pairs <- data.frame(
     Variable1 = rownames(correlation_matrix)[high_cor_indices[, 1]],
     Variable2 = colnames(correlation_matrix)[high_cor_indices[, 2]],
     Correlation = correlation_matrix[high_cor_indices]
   )
   print(high_cor_pairs)
   
  mean(edrr_richness$SpeciesCount)
   var(edrr_richness$SpeciesCount)
   
   m2<- glmmTMB(SpeciesCount ~ Annual_Precipitation +  Minimum_Temperature  + Deciduous + Evergreen + 
                              Shannon_Div_Landscape + POPDEN +  Income_PerCapita + PortDist+  LureName + 
                              (1|TrapSite) + (1|TrapYear) + (1|Ecoregion3), family = poisson(),  data=edrr_richness)
   Anova(m2)
   summary(m2)
   summarym2<- summary(m2)
   coefs_m2 <- as.data.frame(summarym2$coefficients$cond)
   coefs_m2 <- round(coefs_m2, 3)
   coefs_m2$Variable <- rownames(coefs_m2)
   
   
   
  
   
# ---- Model Diagnostics ----   
   
   resm2<- simulateResiduals(m2, plot = TRUE)
   testOutliers(resm2, type = "bootstrap")
  
   
   
   
   
# ---- Multimodel Inference ----
   
   
   
   options(na.action = "na.fail")
   model_set2<- dredge(m2, rank = AIC, trace = TRUE)
   model_averaged2<- model.avg(model_set2, subset=delta < 2) 
   modavg2<- summary(model_averaged2)
   df2<-as.data.frame(modavg2$coefmat.full)
   CI_df2<- as.data.frame(confint(model_averaged2, full=T))
   df2$CI.min <-CI_df2$`2.5 %` 
   df2$CI.max <-CI_df2$`97.5 %`
   df2 <- round(df2, 3)
   df2$Variable <- rownames(df2)
   write.csv(df2, "Data/Saved/Richness_modavg.csv")
   

   saveRDS(model_set2, file = "Data/Saved/Richness_dredge.rds")
    
   # model_set2 <- readRDS("Data/Saved/Richness_dredge.rds")
   


# ---- Plot ----

   
  Richness_overallEDRR<-read.csv(file="Data/Saved/Richness_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
  Richness_overallEDRR$Variable <- gsub("cond\\((.*)\\)", "\\1",Richness_overallEDRR$Variable)
   
  Richness_overallEDRR <- Richness_overallEDRR %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation", "Minimum_Temperature") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", 
                     "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     TRUE ~ "Other"  ))
  Richness_overallEDRR <- Richness_overallEDRR[Richness_overallEDRR$Type != 'Other',]
  Richness_overallEDRR <- Richness_overallEDRR %>% filter(Variable != "(Intercept)")
   
  Richness_overallEDRR <- Richness_overallEDRR %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
  Richness_overallEDRR$Variable <- factor(Richness_overallEDRR$Variable, 
                                          levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                                  "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
  Richness_overallEDRR <- Richness_overallEDRR %>%
    mutate(Variable = dplyr::recode(Variable,
                             "Annual_Precipitation" = "Annual Precipitation",
                             "Minimum_Temperature" = "Minimum Temperature",
                             "Deciduous" = "Deciduous Forest",
                             "Evergreen" = "Evergreen Forest",
                             "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                             "POPDEN" = "Human Population Density",
                             "Income_PerCapita" = "Per Capita Income",
                              "PortDist" = "Distance from Port"
    )) 
  
  
   plot2 <- ggplot(Richness_overallEDRR, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "All Scolytines"
     ) +
     scale_shape_manual(values = c("TRUE" = 15, "FALSE" = 0)) + 
     scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2")) +
     annotate("text",x = Inf,y = Inf,label = "A",hjust = 1.3,vjust = 1.3,fontface = "bold",size = 7) +
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major.y = element_blank(), 
       panel.grid.minor.y = element_blank(),
       scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
       plot.title = element_text(hjust = 0.5, face="bold", size=20),
       strip.text = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
       axis.text.x = element_text(face = "bold", size = 18),
       axis.text.y = element_text(face = "bold", size = 18),
       axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
       axis.title.y = element_text(size = 22, face = "bold", vjust = 2),
       legend.position = c(0.99, 0.02),   # bottom-right INSIDE panel
       legend.justification = c(1, 0),
       legend.background = element_rect(
         fill = "white",
         color = "black",
         linewidth = 0.6
       ),
       legend.key.size = unit(1.2, "cm"),
       legend.text = element_text(size = 14),
       legend.title = element_text(size = 16, face = "bold")
     )+ guides(shape="none")
   
   
   
   plot2
 
   ggsave("Plots/Richness_EDRR.jpg", plot = plot2, width = 10, height = 12, dpi = 600)

###################################################
###################################################
###################################################
###################################################      
# ---- Splitting the data based on trait ----
   
   table(edrr$Reproduction)
   table(edrr$Feeding_Habits)
   table(edrr$Host_range)
   table(edrr$Angio_Gymno)
   
   
   edrr$Species.status <- as.factor(edrr$Species.status)
   #edrr <- edrr %>% filter(!is.na(edrr$Species.status))
   edrr_ambrosia<-edrr[which(edrr$Feeding_Guild=="ambrosia beetle"),]
   edrr_bark<-edrr[which(edrr$Feeding_Guild=="bark beetle"),]
   edrr_nonnative<-edrr[which(edrr$Species.status=="Introduced"),]
   edrr_native<-edrr[which(edrr$Species.status=="Native"),]
   edrr_NA<-edrr[is.na(edrr$Species.status),]
   edrr_nn_bb <- edrr_nonnative[which(edrr_nonnative$Feeding_Guild == "bark beetle"),]
   edrr_nn_ab <- edrr_nonnative[which(edrr_nonnative$Feeding_Guild == "ambrosia beetle"),]
   
   edrr_na_bb <- edrr_native[which(edrr_native$Feeding_Guild == "bark beetle"),]
   edrr_na_ab <- edrr_native[which(edrr_native$Feeding_Guild == "ambrosia beetle"),]
   #edrr_nonnative_barkbeetle <- edrr %>% filter(!is.na(Species.status), !is.na(Feeding_Guild),  Species.status == "Introduced" , Feeding_Guild == "bark beetle")
   
   table(edrr_nonnative$Feeding_Guild)
   edrr_nonnative_NA<-edrr_nonnative[is.na(edrr_nonnative$Feeding_Guild),]
   table(edrr$Feeding_Guild)
   length(unique(edrr_native$Species)) #227
   length(unique(edrr_nonnative$Species)) #45
   length(unique(edrr_NA$Species)) #42
  
  nonnative_species_list <- as.data.frame(unique(edrr_nonnative$Species))
  native_species_list <- as.data.frame(unique(edrr_native$Species))
  print(nonnative_species_list)
  
  # unique(edrr_nonnative$Species)
  # 1            Xylosandrus germanus
  # 2       Xylosandrus crassiusculus
  # 3           Xylosandrus compactus
  # 4           Xylosandrus amputatus
  # 5              Xyleborus seriatus
  # 6                Xyleborus pfeili
  # 7             Xyleborus glabratus
  # 8            Xyleborus bispinatus
  # 9           Xyleborinus saxesenii
  # 10     Xyleborinus octiesdentatus
  # 11         Xyleborinus attenuatus
  # 12       Xyleborinus artestriatus
  # 13          Xyleborinus andrewesi
  # 14        Trypodendron domesticum
  # 15              Tomicus piniperda
  # 16            Scolytus schevyrewi
  # 17             Scolytus rugulosus
  # 18         Scolytus multistriatus
  # 19          Pityogenes bidentatus
  # 20            Orthotomicus erosus
  # 21           Hypothenemus erectus
  # 22            Hylurgus ligniperda
  # 23            Hylurgops palliatus
  # 24            Hylastinus obscurus
  # 25                Hylastes opacus
  # 26             Euwallacea validus
  # 27             Euwallacea similis
  # 28         Euwallacea interjectus
  # 29          Euwallacea fornicatus
  # 30          Dryoxylon onoharaense
  # 31         Dryocoetes autographus
  # 32      Dactylotrypes longicollis
  # 33    Cyclorhipidion pelliculosum
  # 34       Cyclorhipidion inarmatum
  # 35   Cyclorhipidion distinguendum
  # 36        Cyclorhipidion bodoanum
  # 37            Crypturgus pusillus
  # 38        Coptoborus pseudotenuis
  # 39       Coccotrypes dactyliperda
  # 40              Cnestus mutilatus
  # 41              Anisandrus maiche
  # 42              Anisandrus dispar
  # 43         Ambrosiophilus atratus
  # 44       Ambrosiodmus rubricollis
  # 45             Ambrosiodmus minor
  
  NA_species_list <- as.data.frame(unique(edrr_NA$Species))
  
  print(NA_species_list)
  

  # unique(edrr_NA$Species)
  # 1                  Xyleborus sp.
  # 2         Wallacellus denticulus
  # 3          Trypodendron borealis
  # 4                   Trypodendron
  # 5                    Trischidias
  # 6                      Thysanoes
  # 7               Theoborus ricini
  # 8                       Scolytus
  # 9               Scolytogenes sp.
  # 10          Scolytogenes jalapae
  # 11           Scolytodes schwarzi
  # 12               Pseudothysanoes
  # 13           Pseudopityophthorus
  # 14            Pseudips mexicanus
  # 15                Premnobius spp
  # 16             Pityophthorus spp
  # 17                 Pityophthorus
  # 18                  Pityokteines
  # 19                    Pityogenes
  # 20                    Pityoborus
  # 21                  Phloeotribus
  # 22 Phloeosinus scopulorum s.str.
  # 23                   Phloeosinus
  # 24                  Micracisella
  # 25                      Micracis
  # 26                  Hypothenemus
  # 27             Hylocurus torosus
  # 28                     Hylocurus
  # 29                     Hylesinus
  # 30                      Hylastes
  # 31       Euwallacea denticulatus
  # 32                    Dryocoetes
  # 33       Cyclorhipidion fukiense
  # 34                    Crypturgus
  # 35                     Cryphalus
  # 36            Corthylus spinifer
  # 37              Conophthorus sp.
  # 38                   Coccotrypes
  # 39                     Chramesus
  # 40                Carphoborus sp
  # 41      Ambrosiophilus nodulosis
  # 42              Ambrosiodmus sp.
  

  
  
   

#################################################################
#################################################################
# ---- Nonnative Aggregated ----
   # Nonnative all
   # Summarize counts
   edrr_long_nn <- edrr_nonnative %>% group_by(TrapID, Species) %>% summarise(count = sum(PestCount, na.rm = TRUE), .groups = "drop")
   
   # wide format
   edrr_wide_nn <- edrr_long_nn %>% pivot_wider(names_from = Species, values_from = count, values_fill = 0)
   
   # Extract metadata , keeps first value for each variable
   meta_info_nn <- edrr_nonnative %>%  group_by(TrapID) %>%
     summarise(across(c(Cell_ID, SiteLure, StateCode, CountyName,SiteTown,SiteLoc,SiteType, TrapYear, SiteNumber, LureName, FIPS, Longitude, 
                        Latitude, TrapSite, landcover, LC_class,  Ecoregion3, Annual_Precipitation , Minimum_Temperature,  Deciduous, Evergreen,  
                        Shannon_Div_Landscape, POPDEN,  Income_PerCapita, PortDist), first), .groups = "drop")
   
   
   
   edrr_abundance_nn <- edrr_long_nn %>% left_join(meta_info_nn, by ="TrapID")
  
   edrr_richness_nn <- edrr_wide_nn %>% left_join(meta_info_nn, by = "TrapID")
   
   
   #Calculating Species Count for each Trap ID
   edrr_richness_nn$SpeciesCount<-rowSums(edrr_richness_nn[, 2:46 ] !=0)
   

   write.csv(edrr_abundance_nn, "Data/Saved/edrr_abundance_nn.csv", row.names = FALSE)
   write.csv(edrr_richness_nn, "Data/Saved/edrr_richness_nn.csv", row.names = FALSE)
   
# ---- Checking correlation among predictors ----   
   pred_vars_nn <- edrr_abundance_nn[, c("Annual_Precipitation","Minimum_Temperature", "Evergreen", "Deciduous",
                                      "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita",  "PortDist")]
   cor_matrix_nn <- cor(pred_vars_nn, use = "complete.obs", method = "pearson")
   correlation_matrix_nn <- as.data.frame(cor_matrix_nn)
  
  
   
   high_cor_indices <- which(abs(correlation_matrix_nn) >= 0.7 & abs(correlation_matrix_nn) < 1, arr.ind = TRUE)
   high_cor_pairs <- data.frame(
     Variable1 = rownames(correlation_matrix_nn)[high_cor_indices[, 1]],
     Variable2 = colnames(correlation_matrix_nn)[high_cor_indices[, 2]],
     Correlation = correlation_matrix_nn[high_cor_indices]
   )
   print(high_cor_pairs)
  
# ---- Min temp and decid forest models ----
  
   mean(edrr_abundance_nn$count)
   var(edrr_abundance_nn$count)
   

   
    mintemp<- glmmTMB(count ~  Minimum_Temperature + Annual_Precipitation + Evergreen + 
                       Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist +  LureName + 
                     (1|TrapSite/TrapID) + (1|TrapYear) + (1|Species) + (1|Ecoregion3), 
                     family = nbinom2(), data=edrr_abundance_nn)
   
   Anova(mintemp)
   summary(mintemp)
   
  
# ---- Model diagnostics ---- 
   
   sim_res_mintemp <- simulateResiduals(mintemp, plot = TRUE)
   testOutliers(sim_res_mintemp, type = "bootstrap")
  
   
   
   decid<- glmmTMB(count ~  Deciduous + Annual_Precipitation + Evergreen + 
                     Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName + 
                       (1|TrapSite/TrapID) + (1|TrapYear) + (1|Species) + (1|Ecoregion3), 
                   family = nbinom2(), data=edrr_abundance_nn)
   
   Anova(decid)
   summary(decid)
   
   
   AIC(mintemp, decid)
   
   
   fml_mintemp <- Reduce(paste, deparse(formula(mintemp)))
   fml_decid   <- Reduce(paste, deparse(formula(decid)))
   aic_vals <- AIC(mintemp, decid)$AIC
   
   
   modelcomp1 <- data.frame(
     Group = c("Non-native Scolytine Abundance", "Non-native Scolytine Abundance"),
     Model = c(
       gsub("count ~ ","", fml_mintemp),
       gsub("count ~ ","", fml_decid)
     ),
     AIC = round(aic_vals, 2),
     Delta_AIC = round(aic_vals - min(aic_vals), 2)
   )
   modelcomp1 <- modelcomp1[order(modelcomp1$AIC), ]
   modelcomp1$Group[duplicated(modelcomp1$Group)] <- ""
   
   names(modelcomp1)[names(modelcomp1) == "Delta_AIC"] <- "ΔAIC"
   
  
# ---- Full Model Abundance ----
  

   
   m1_nn<- glmmTMB(count ~ Annual_Precipitation + Deciduous + Evergreen + 
                                    Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName + 
                                    (1|TrapSite/TrapID) + (1|TrapYear) + (1|Species) + (1|Ecoregion3),
                                    family = nbinom2(), data=edrr_abundance_nn)
   
   Anova(m1_nn)
   summary(m1_nn)

   summarym1_nn<- summary(m1_nn)
   coefs_m1_nn <- as.data.frame(summarym1_nn$coefficients$cond)
   coefs_m1_nn <- round(coefs_m1_nn, 3)
   coefs_m1_nn$Variable <- rownames(coefs_m1_nn)
  
   
  
# ---- Model Diagnostics ----

   
   sim_res_m1_nn <- simulateResiduals(m1_nn, plot = TRUE)
   testOutliers(sim_res_m1_nn, type = "bootstrap")
 
   
   
  

# ---- Multi-model inference- Abundance ---- 
  
   options(na.action = "na.fail")
   dredge_m1_nn <- dredge(m1_nn, rank = AIC, trace = TRUE)
  
   model_averaged_nn<- model.avg(dredge_m1_nn, subset=delta < 2) 
   modavg_nn<- summary(model_averaged_nn)
   df_nn <-as.data.frame(modavg_nn$coefmat.full)
   CI_df_nn<- as.data.frame(confint(model_averaged_nn, full=T))
   df_nn$CI.min <-CI_df_nn$`2.5 %`
   df_nn$CI.max <-CI_df_nn$`97.5 %`
   df_nn <- round(df_nn, 3)
   df_nn$Variable <- rownames(df_nn)
   write.csv(df_nn, "Data/Saved/Abundance_nonnative_modavg.csv")
   
   saveRDS(dredge_m1_nn, file = "Data/Saved/dredge_m1_nn.rds")
# ---- Plot ----  


   Abundance_nonnative<-read.csv(file="Data/Saved/Abundance_nonnative_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
   Abundance_nonnative$Variable <- gsub("cond\\((.*)\\)", "\\1",Abundance_nonnative$Variable)
   
   Abundance_nonnative <- Abundance_nonnative %>%
     bind_rows(
       tibble(
         Variable = "Minimum_Temperature",
         Estimate = 0,
         CI.min = 0,
         CI.max = 0,
         Pr...z.. = NA,
         Type = "Removed",
         p_value = NA,
         Significant = NA
       )
     )
   
   
   Abundance_nonnative <- Abundance_nonnative %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     Variable %in% c("Minimum_Temperature") ~ "Removed", 
     TRUE ~ "Other"  ))
   Abundance_nonnative <- Abundance_nonnative[Abundance_nonnative$Type != 'Other',]
   
   
   Abundance_nonnative <- Abundance_nonnative %>% filter(Variable != "(Intercept)")
   
   Abundance_nonnative <- Abundance_nonnative %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
   Abundance_nonnative$Variable <- factor(Abundance_nonnative$Variable,
                                          levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                                     "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
   
   
   Abundance_nonnative <- Abundance_nonnative %>%
     mutate(Variable = dplyr::recode(Variable,
                              "Annual_Precipitation" = "Annual Precipitation",
                              "Minimum_Temperature" = "Minimum Temperature",
                              "Deciduous" = "Deciduous Forest",
                              "Evergreen" = "Evergreen Forest",
                              "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                              "POPDEN" = "Human Population Density",
                              "Income_PerCapita" = "Per Capita Income",
                              "PortDist" = "Distance from Port"
     ))
   plot_nn_abd <- ggplot(Abundance_nonnative, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "Non-native Scolytines"
     ) +
     scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 2), na.value = 4) +
     scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2", "Removed" = "gray50"),
                        breaks = c("Climatic", "Landscape", "Socioeconomic") ) +
     annotate("text",x = Inf,y = Inf,label = "C",hjust = 1.3,vjust = 1.3,fontface = "bold",size = 7) +
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major.y = element_blank(), 
       panel.grid.minor.y = element_blank(),
       scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
       plot.title = element_text(hjust = 0.5, face="bold", size=20),
       strip.text = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
       axis.text.x = element_text(face = "bold", size = 18),
       axis.text.y = element_blank(),
       axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
       axis.title.y = element_text(size = 22, face = "bold", vjust = 2),
       legend.position = "none",
     )+ guides(shape="none")
   
   
   plot_nn_abd
   
   ggsave("Plots/Abundance_nonnative.jpg", plot = plot_nn_abd, width = 10, height = 12, dpi = 600)
# ---- Full Model Richness ---- 
   pred_vars_nn <- edrr_richness_nn[, c("Annual_Precipitation","Minimum_Temperature", "Evergreen", "Deciduous",
                                         "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist")]
   cor_matrix <- cor(pred_vars_nn, use = "complete.obs", method = "pearson")
   correlation_matrix_nn <- as.data.frame(cor_matrix)
   
   
   
   high_cor_indices <- which(abs(correlation_matrix_nn) >= 0.7 & abs(correlation_matrix_nn) < 1, arr.ind = TRUE)
   high_cor_pairs <- data.frame(
     Variable1 = rownames(correlation_matrix_nn)[high_cor_indices[, 1]],
     Variable2 = colnames(correlation_matrix_nn)[high_cor_indices[, 2]],
     Correlation = correlation_matrix_nn[high_cor_indices]
   )
   print(high_cor_pairs)
   hist(edrr_richness_nn$SpeciesCount)
   mean(edrr_richness_nn$SpeciesCount)
   var(edrr_richness_nn$SpeciesCount)
  
   m2_nn<- glmmTMB(SpeciesCount ~ Annual_Precipitation +  Minimum_Temperature + Deciduous + Evergreen + 
                                 Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName + 
                                (1|TrapSite) + (1|TrapYear) + (1|Ecoregion3), family = poisson(), data=edrr_richness_nn)
   
   
   
   
   Anova(m2_nn)
   summary(m2_nn)
   summarym2_nn<- summary(m2_nn)
   coefs_m2_nn <- as.data.frame(summarym2_nn$coefficients$cond)
   coefs_m2_nn <- round(coefs_m2_nn, 3)
   coefs_m2_nn$Variable <- rownames(coefs_m2_nn)
   write.csv(coefs_m2_nn, "data/saved/Richness~X_nonnativeEDRR.csv")
   
   
   
# ---- Model Diagnostics ----   
   sim_res_m2_nn <- simulateResiduals(m2_nn, plot = TRUE)
   testOutliers(sim_res_m2_nn, type = "bootstrap")
  
# ---- Multi-model inference- Richness ----    
   options(na.action = "na.fail")
   dredge_m2_nn <- dredge(m2_nn, rank = AIC, trace = TRUE)
   model_averaged_nn2<- model.avg(dredge_m2_nn, subset=delta < 2) 
   modavg_nn2<- summary(model_averaged_nn2)
   df_nn2 <-as.data.frame(modavg_nn2$coefmat.full)
   CI_df_nn2<- as.data.frame(confint(model_averaged_nn2, full=T))
   df_nn2$CI.min <-CI_df_nn2$`2.5 %` 
   df_nn2$CI.max <-CI_df_nn2$`97.5 %`
   df_nn2 <- round(df_nn2, 3)
   df_nn2$Variable <- rownames(df_nn2)
   write.csv(df_nn2, "Data/Saved/Richness_nonnative_modavg.csv")
   
    saveRDS(dredge_m2_nn, file = "Data/Saved/dredge_m2_nn.rds")
    
    #dredge_m2_nn <- readRDS("Data/Saved/dredge_m2_nn.rds")


   
   
   
   
   
# ---- Plot ----  
   
   
   Richness_nonnative<-read.csv(file="Data/Saved/Richness_nonnative_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
   Richness_nonnative$Variable <- gsub("cond\\((.*)\\)", "\\1",Richness_nonnative$Variable)
   
   
   Richness_nonnative <- Richness_nonnative %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation", "Minimum_Temperature") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     TRUE ~ "Other"  ))
   Richness_nonnative <- Richness_nonnative[Richness_nonnative$Type != 'Other',]
   
   
   Richness_nonnative <- Richness_nonnative %>% filter(Variable != "(Intercept)")
   
   Richness_nonnative <- Richness_nonnative %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
   Richness_nonnative$Variable <- factor(Richness_nonnative$Variable, 
                                         levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                                      "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
   
   
   Richness_nonnative <- Richness_nonnative %>%
     mutate(Variable = dplyr::recode(Variable,
                              "Annual_Precipitation" = "Annual Precipitation",
                              "Minimum_Temperature" = "Minimum Temperature",
                              "Deciduous" = "Deciduous Forest",
                              "Evergreen" = "Evergreen Forest",
                              "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                              "POPDEN" = "Human Population Density",
                              "Income_PerCapita" = "Per Capita Income",
                              "PortDist" = "Distance from Port"
     ))
   plot_nn_ric <- ggplot(Richness_nonnative, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "Non-native Scolytines"
     ) +
     scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 2)) + 
     scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2")) +
     annotate("text",x = Inf,y = Inf,label = "C",hjust = 1.3,vjust = 1.3,fontface = "bold",size = 7) +
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major.y = element_blank(), 
       panel.grid.minor.y = element_blank(),
       scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
       plot.title = element_text(hjust = 0.5, face="bold", size=20),
       strip.text = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
       axis.text.x = element_text(face = "bold", size = 18),
       axis.text.y = element_blank(),
       axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
       axis.title.y = element_text(size = 22, face = "bold", vjust = 2),
       legend.position = "none",
     )+ guides(shape="none")
   
   
   
   
   plot_nn_ric
   
   ggsave("Plots/Richness_nonnative.jpg", plot = plot_nn_ric, width = 10, height = 12, dpi = 600)
###############################
###############################
# ---- Nonnative Ambrosia Beetles ----
   
   # Summarize counts
   edrr_long_nn_ab <- edrr_nn_ab %>%
     group_by(TrapID, Species) %>%
     summarise(count = sum(PestCount, na.rm = TRUE), .groups = "drop")
   
   # wide format
   edrr_wide_nn_ab <- edrr_long_nn_ab %>%
     pivot_wider(
       names_from = Species,
       values_from = count,
       values_fill = 0  # Replace NA with 0 
     )
   
   # Extract metadata , keeps first value for each variable
   meta_info_nn_ab <- edrr_nn_ab %>%
     group_by(TrapID) %>%
     summarise(across(c(Cell_ID, SiteLure, StateCode, CountyName,SiteTown,SiteLoc,SiteType, TrapYear, SiteNumber, LureName, FIPS, Longitude, Latitude,
                        TrapSite, landcover, LC_class,  Ecoregion3,  
                        Annual_Precipitation , Minimum_Temperature,  Deciduous, Evergreen,  
                        Shannon_Div_Landscape, POPDEN,  Income_PerCapita, PortDist), first), .groups = "drop")
   
   
   edrr_abundance_nn_ab <- edrr_long_nn_ab %>% left_join(meta_info_nn_ab, by = "TrapID")
   
   edrr_richness_nn_ab <- edrr_wide_nn_ab %>%  left_join(meta_info_nn_ab, by = "TrapID")
   
   

   
   
   #Calculating Species Count for each Trap ID
   edrr_richness_nn_ab$SpeciesCount<-rowSums(edrr_richness_nn_ab[, 2:31] !=0)
   
 
   write.csv(edrr_abundance_nn_ab, "Data/Saved/edrr_abundance_nn_ab.csv", row.names = FALSE)
   write.csv(edrr_richness_nn_ab, "Data/Saved/edrr_richness_nn_ab.csv", row.names = FALSE)
  
# ---- Checking correlation among predictors ----   
   pred_vars_nn_ab <- edrr_abundance_nn_ab[, c("Annual_Precipitation","Minimum_Temperature", "Evergreen", "Deciduous",
                                            "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist")]
   cor_matrix_nn_ab <- cor(pred_vars_nn_ab, use = "complete.obs", method = "pearson")
   correlation_matrix_nn_ab <- as.data.frame(cor_matrix_nn_ab)
   
   
   high_cor_indices <- which(abs(correlation_matrix_nn_ab) >= 0.7 & abs(correlation_matrix_nn_ab) < 1, arr.ind = TRUE)
   high_cor_pairs <- data.frame(
     Variable1 = rownames(correlation_matrix_nn_ab)[high_cor_indices[, 1]],
     Variable2 = colnames(correlation_matrix_nn_ab)[high_cor_indices[, 2]],
     Correlation = correlation_matrix_nn_ab[high_cor_indices]
   )
   print(high_cor_pairs)
   
# ---- Checking correlated predictors individually ----
  
   
    mintemp<- glmmTMB(count ~ Minimum_Temperature + Annual_Precipitation + Evergreen + 
                       Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName + (1|TrapSite/TrapID) + (1|Species)
                    + (1|TrapYear) +(1|Ecoregion3), family = nbinom2(),
                    data=edrr_abundance_nn_ab)
   Anova(mintemp)
   summary(mintemp)
   
#----  Model Diagnostics ----
  
   
   sim_res_mintemp <- simulateResiduals(mintemp, plot = TRUE)
   testOutliers(sim_res_mintemp, type = "bootstrap")
   testDispersion(sim_res_mintemp)   
   
   
   
   decid<- glmmTMB(count ~ Deciduous + Annual_Precipitation + Evergreen + 
                      Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName + (1|TrapSite/TrapID) + (1|Species)
                    + (1|TrapYear) +(1|Ecoregion3), family = nbinom2(),
                    data=edrr_abundance_nn_ab)
   Anova(decid)
   summary(decid)
   
   
   AIC(mintemp, decid)
   
   fml_mintemp <- Reduce(paste, deparse(formula(mintemp)))
   fml_decid   <- Reduce(paste, deparse(formula(decid)))
   aic_vals <- AIC(mintemp, decid)$AIC
   
   
   modelcomp2 <- data.frame(
     Group = c("Non-native Ambrosia Beetle Abundance", "Non-native Ambrosia Beetle Abundance"),
     Model = c(
       gsub("count ~ ","", fml_mintemp),
       gsub("count ~ ","", fml_decid)
     ),
     AIC = round(aic_vals, 2),
     Delta_AIC = round(aic_vals - min(aic_vals), 2)
   )
   modelcomp2 <- modelcomp2[order(modelcomp2$AIC), ]
   modelcomp2$Group[duplicated(modelcomp2$Group)] <- ""
   
   names(modelcomp2)[names(modelcomp2) == "Delta_AIC"] <- "ΔAIC"
   
  
   
# ---- Full Model- Abundance----
   
   abundance_nn_ambrosia<- glmmTMB(count ~ Annual_Precipitation +  Deciduous  + Evergreen + 
                                                       Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName +
                                     (1|TrapSite/TrapID) + (1|Species) + (1|TrapYear) + (1|Ecoregion3), family = nbinom2(),
                                   data=edrr_abundance_nn_ab)
   Anova(abundance_nn_ambrosia)
   summary(abundance_nn_ambrosia)
 

   summary_abd_nn_ab<- summary(abundance_nn_ambrosia)
   coefs_m1_nn_ab <- as.data.frame(summary_abd_nn_ab$coefficients$cond)
   coefs_m1_nn_ab <- round(coefs_m1_nn_ab, 3)
   coefs_m1_nn_ab$Variable <- rownames(coefs_m1_nn_ab)
   

# ---- Multi-model inference- Abundance ---- 
   
   options(na.action = "na.fail")
   dredge_abd_nn_ab <- dredge(abundance_nn_ambrosia, rank = AIC, trace = TRUE)
   model_averaged_nn_ab<- model.avg(dredge_abd_nn_ab, subset=delta < 2) 
   modavg_nn_ab<- summary(model_averaged_nn_ab)
   df_nn_ab <-as.data.frame(modavg_nn_ab$coefmat.full)
   CI_df_nn_ab<- as.data.frame(confint(model_averaged_nn_ab, full=T))
   df_nn_ab$CI.min <-CI_df_nn_ab$`2.5 %` 
   df_nn_ab$CI.max <-CI_df_nn_ab$`97.5 %`
   df_nn_ab <- round(df_nn_ab, 3)
   df_nn_ab$Variable <- rownames(df_nn_ab)
   write.csv(df_nn_ab, "Data/Saved/Abundance_nonnativeAmbrosia_modavg.csv")
   

  saveRDS(dredge_abd_nn_ab, "Data/Saved/Abundance_nonnativeAmbrosia_dredge.rds")
  #  dredge_nonnative<-readRDS(file="Data/Saved/Abundance_nonnative_dredge.rds)
  #  
   
   

# ---- Plot ----  
   
   
   Abundance_nonnative_ambrosia<-read.csv(file="Data/Saved/Abundance_nonnativeAmbrosia_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
   Abundance_nonnative_ambrosia$Variable <- gsub("cond\\((.*)\\)", "\\1",Abundance_nonnative_ambrosia$Variable)
   Abundance_nonnative_ambrosia <- Abundance_nonnative_ambrosia %>%
     bind_rows(
       tibble(
         Variable = "Minimum_Temperature",
         Estimate = 0,
         CI.min = 0,
         CI.max = 0,
         Pr...z.. = NA,
         Type = "Removed",
         p_value = NA,
         Significant = NA
       )
     )
   
   Abundance_nonnative_ambrosia <- Abundance_nonnative_ambrosia %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     Variable %in% c("Minimum_Temperature") ~ "Removed", 
     TRUE ~ "Other"  ))
   Abundance_nonnative_ambrosia <- Abundance_nonnative_ambrosia[Abundance_nonnative_ambrosia$Type != 'Other',]
   
   
   Abundance_nonnative_ambrosia <- Abundance_nonnative_ambrosia %>% filter(Variable != "(Intercept)")
   
   Abundance_nonnative_ambrosia <- Abundance_nonnative_ambrosia %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
  
   Abundance_nonnative_ambrosia$Variable <- factor(Abundance_nonnative_ambrosia$Variable, 
                                                   levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                                                                   "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
   
   
   Abundance_nonnative_ambrosia <- Abundance_nonnative_ambrosia %>%
     mutate(Variable = dplyr::recode(Variable,
                              "Annual_Precipitation" = "Annual Precipitation",
                              "Minimum_Temperature" = "Minimum Temperature",
                              "Deciduous" = "Deciduous Forest",
                              "Evergreen" = "Evergreen Forest",
                              "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                              "POPDEN" = "Human Population Density",
                              "Income_PerCapita" = "Per Capita Income",
                              "PortDist" = "Distance from Port"
     ))
   

   
   plot_nn_ambrosia_abd <- ggplot(Abundance_nonnative_ambrosia, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "Non-native Ambrosia Beetles"
     ) +
     scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 2), na.value = 4) +
     scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2", "Removed" = "gray50"),
                        breaks = c("Climatic", "Landscape", "Socioeconomic") ) +
     annotate("text",x = Inf,y = Inf,label = "E",hjust = 1.3,vjust = 1.3,fontface = "bold",size = 7) +
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major.y = element_blank(), 
       panel.grid.minor.y = element_blank(),
       scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
       plot.title = element_text(hjust = 0.5, face="bold", size=20),
       strip.text = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
       axis.text.x = element_text(face = "bold", size = 18),
       axis.text.y = element_blank(),
       axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
       axis.title.y = element_text(size = 22, face = "bold", vjust = 2),
       legend.position = "none",
     )+ guides(shape="none")
   
   
   plot_nn_ambrosia_abd
   
   ggsave("Plots/Abundance_nonnative_ambrosia.jpg", plot = plot_nn_ambrosia_abd, width = 10, height = 12, dpi = 600)
# ---- Full Model- Richness----
   pred_vars_nn_ab <- edrr_richness_nn_ab[, c("Annual_Precipitation","Minimum_Temperature", "Evergreen", "Deciduous",
                                               "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist")]
   cor_matrix_nn_ab <- cor(pred_vars_nn_ab, use = "complete.obs", method = "pearson")
   correlation_matrix_nn_ab <- as.data.frame(cor_matrix_nn_ab)
   
   
   high_cor_indices <- which(abs(correlation_matrix_nn_ab) >= 0.7 & abs(correlation_matrix_nn_ab) < 1, arr.ind = TRUE)
   high_cor_pairs <- data.frame(
     Variable1 = rownames(correlation_matrix_nn_ab)[high_cor_indices[, 1]],
     Variable2 = colnames(correlation_matrix_nn_ab)[high_cor_indices[, 2]],
     Correlation = correlation_matrix_nn_ab[high_cor_indices]
   )
   print(high_cor_pairs)
   
   
   hist(edrr_richness_nn_ab$SpeciesCount)
   mean(edrr_richness_nn_ab$SpeciesCount)
   var(edrr_richness_nn_ab$SpeciesCount)
   
   decid<- glmmTMB(SpeciesCount ~ Deciduous  +  Annual_Precipitation + 
                                    Evergreen + Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName +
                                    (1|TrapSite) + (1|TrapYear) + (1|Ecoregion3), family = poisson(), 
                                  data=edrr_richness_nn_ab)
   
   
   
   Anova(decid)
   summary(decid)
   
   
   ######  Model Diagnostics 
   
   
   sim_res_decid <- simulateResiduals(decid, plot = TRUE)
   testOutliers(sim_res_decid, type = "bootstrap")
   
   mintemp<- glmmTMB(SpeciesCount ~  Minimum_Temperature  + Annual_Precipitation + 
                     Evergreen + Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName +
                     (1|TrapSite) + (1|TrapYear) + (1|Ecoregion3), family = poisson(), 
                   data=edrr_richness_nn_ab)
   
   
   
   Anova(mintemp)
   summary(mintemp)
   
   
   AIC(mintemp, decid)
   
   fml_mintemp <- Reduce(paste, deparse(formula(mintemp)))
   fml_decid   <- Reduce(paste, deparse(formula(decid)))
   aic_vals <- AIC(mintemp, decid)$AIC
   
   
   modelcomp3 <- data.frame(
     Group = c("Non-native Ambrosia Beetle Richness", "Non-native Ambrosia Beetle Richness"),
     Model = c(
       gsub("SpeciesCount ~ ","", fml_mintemp),
       gsub("SpeciesCount ~ ","", fml_decid)
     ),
     AIC = round(aic_vals, 2),
     Delta_AIC = round(aic_vals - min(aic_vals), 2)
   )
   modelcomp3 <- modelcomp3[order(modelcomp3$AIC), ]
   modelcomp3$Group[duplicated(modelcomp3$Group)] <- ""
   
   names(modelcomp3)[names(modelcomp3) == "Delta_AIC"] <- "ΔAIC"
   
   
   
  
   
   richness_nn_ambrosia<- glmmTMB(SpeciesCount ~ Annual_Precipitation +  Deciduous  + 
                                     Evergreen + Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName +
                                     (1|TrapSite) + (1|TrapYear) + (1|Ecoregion3), family = poisson(), 
                                        data=edrr_richness_nn_ab)
   
   Anova(richness_nn_ambrosia)
   summary(richness_nn_ambrosia)
   
   saveRDS(richness_nn_ambrosia, "Data/Saved/edrr_richness_nn_ab.rds")
   
   summary_ric_nn_ab<- summary(richness_nn_ambrosia)
   coefs_m2_nn_ab <- as.data.frame(summary_ric_nn_ab$coefficients$cond)
   coefs_m2_nn_ab <- round(coefs_m2_nn_ab, 3)
   coefs_m2_nn_ab$Variable <- rownames(coefs_m2_nn_ab)
 

# ---- Multi-model inference- Richness ---- 
   
   options(na.action = "na.fail")
   dredge_ric_nn_ab <- dredge(richness_nn_ambrosia, rank = AIC, trace = TRUE)
   model_averaged_nn2_ab<- model.avg(dredge_ric_nn_ab, subset=delta < 2) 
   modavg_nn2_ab<- summary(model_averaged_nn2_ab)
   df_nn2_ab <-as.data.frame(modavg_nn2_ab$coefmat.full)
   CI_df_nn2_ab<- as.data.frame(confint(model_averaged_nn2_ab, full=T))
   df_nn2_ab$CI.min <-CI_df_nn2_ab$`2.5 %` 
   df_nn2_ab$CI.max <-CI_df_nn2_ab$`97.5 %`
   df_nn2_ab <- round(df_nn2_ab, 3)
   df_nn2_ab$Variable <- rownames(df_nn2_ab)
   write.csv(df_nn2_ab, "Data/Saved/Richness_nonnativeAmbrosia_modavg.csv")
   
  
    saveRDS(dredge_ric_nn_ab, "Data/Saved/Richness_nonnativeAmbrosia_dredge.rds")
   # dredge_nonnative<-readRDS(file="Data/Saved/Abundance_nonnative_dredge.csv")
   # 
   
   
   
# ---- Plot ----  
   
   
   Richness_nonnative_ambrosia<-read.csv(file="Data/Saved/Richness_nonnativeAmbrosia_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
   Richness_nonnative_ambrosia$Variable <- gsub("cond\\((.*)\\)", "\\1",Richness_nonnative_ambrosia$Variable)
   
   Richness_nonnative_ambrosia <- Richness_nonnative_ambrosia %>%
     bind_rows(
       tibble(
         Variable = "Minimum_Temperature",
         Estimate = 0,
         CI.min = 0,
         CI.max = 0,
         Pr...z.. = NA,
         Type = "Removed",
         p_value = NA,
         Significant = NA
       )
     )
   Richness_nonnative_ambrosia <- Richness_nonnative_ambrosia %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     Variable %in% c("Minimum_Temperature") ~ "Removed", 
     TRUE ~ "Other"  ))
   Richness_nonnative_ambrosia <- Richness_nonnative_ambrosia[Richness_nonnative_ambrosia$Type != 'Other',]
   
   
   Richness_nonnative_ambrosia <- Richness_nonnative_ambrosia %>% filter(Variable != "(Intercept)")
   
   Richness_nonnative_ambrosia <- Richness_nonnative_ambrosia %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
   Richness_nonnative_ambrosia$Variable <- factor(Richness_nonnative_ambrosia$Variable, 
                                                  levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                                                  "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
   
   
   Richness_nonnative_ambrosia <- Richness_nonnative_ambrosia %>%
     mutate(Variable = dplyr::recode(Variable,
                              "Annual_Precipitation" = "Annual Precipitation",
                              "Minimum_Temperature" = "Minimum Temperature",
                              "Deciduous" = "Deciduous Forest",
                              "Evergreen" = "Evergreen Forest",
                              "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                              "POPDEN" = "Human Population Density",
                              "Income_PerCapita" = "Per Capita Income", 
                              "PortDist" = "Distance from Port"
     ))
   plot_nn_ambrosia_ric <- ggplot(Richness_nonnative_ambrosia, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "Non-native Ambrosia Beetles"
     ) +
     scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 2), na.value = 4) + 
     scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2", "Removed" = "gray50")) +
     annotate("text",x = Inf,y = Inf,label = "E",hjust = 1.3,vjust = 1.3,fontface = "bold",size = 7) +
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major.y = element_blank(), 
       panel.grid.minor.y = element_blank(),
       scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
       plot.title = element_text(hjust = 0.5, face="bold", size=20),
       strip.text = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
       axis.text.x = element_text(face = "bold", size = 18),
       axis.text.y = element_blank(),
       axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
       axis.title.y = element_text(size = 22, face = "bold", vjust = 2),
       legend.position = "none",
     )+ guides(shape="none")
   
   
   plot_nn_ambrosia_ric
   
   ggsave("Plots/Richness_nonnative_ambrosia.jpg", plot = plot_nn_ambrosia_ric, width = 10, height = 12, dpi = 600)
   
############################################
#############################################
# ---- Nonnative Bark Beetles ----  
   
   # Summarize counts
   edrr_long_nn_bb <- edrr_nn_bb %>%
     group_by(TrapID, Species) %>%
     summarise(count = sum(PestCount, na.rm = TRUE), .groups = "drop")
   
   # wide format
   edrr_wide_nn_bb <- edrr_long_nn_bb %>%
     pivot_wider(
       names_from = Species,
       values_from = count,
       values_fill = 0  # Replace NA with 0 
     )
   
   # Extract metadata , keeps first value for each variable
   meta_info_nn_bb <- edrr_nn_bb %>%
     group_by(TrapID) %>%
     summarise(across(c(Cell_ID, SiteLure, StateCode, CountyName,SiteTown,SiteLoc,SiteType, TrapYear, SiteNumber, LureName, FIPS, Longitude, Latitude,
                        TrapSite, landcover, LC_class,  Ecoregion3,  
                        Annual_Precipitation , Minimum_Temperature,  Deciduous, Evergreen,  
                        Shannon_Div_Landscape, POPDEN,  Income_PerCapita, PortDist), first), .groups = "drop")
   
   
   
   
   edrr_abundance_nn_bb <- edrr_long_nn_bb %>% left_join(meta_info_nn_bb, by = "TrapID")
   
   edrr_richness_nn_bb <- edrr_wide_nn_bb %>%  left_join(meta_info_nn_bb, by = "TrapID")
   
  length(unique(edrr_nn_bb$Species)) 
 
  
   
   
   #Calculating Species Count for each Trap ID
  edrr_richness_nn_bb$SpeciesCount<-rowSums(edrr_richness_nn_bb[, 2:12] !=0)
   
  
  write.csv(edrr_abundance_nn_bb, "Data/Saved/edrr_abundance_nn_bb.csv", row.names = FALSE)
  write.csv(edrr_richness_nn_bb, "Data/Saved/edrr_richness_nn_bb.csv", row.names = FALSE)
 length(unique(edrr_richness_nn_bb$TrapSite))
        
 str(edrr_richness_nn_bb)
   
# ---- Checking correlation among predictors ----   
   pred_vars_nn_bb <- edrr_abundance_nn_bb[, c("Annual_Precipitation","Minimum_Temperature", "Evergreen", "Deciduous",
                                            "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist")]
   cor_matrix_nn_bb <- cor(pred_vars_nn_bb, use = "complete.obs", method = "pearson")
   correlation_matrix_nn_bb <- as.data.frame(cor_matrix_nn_bb)
   
   
   high_cor_indices <- which(abs(correlation_matrix_nn_bb) >= 0.7 & abs(correlation_matrix_nn_bb) < 1, arr.ind = TRUE)
   high_cor_pairs <- data.frame(
     Variable1 = rownames(correlation_matrix_nn_bb)[high_cor_indices[, 1]],
     Variable2 = colnames(correlation_matrix_nn_bb)[high_cor_indices[, 2]],
     Correlation = correlation_matrix_nn_bb[high_cor_indices]
   )
   print(high_cor_pairs)
   

# ---- Full Model- Abundance----
  
   mean(edrr_abundance_nn_bb$count)
   var(edrr_abundance_nn_bb$count)
   hist(edrr_abundance_nn_bb$count)
   abundance_nn_bark<- glmmTMB(count ~ Annual_Precipitation +  Minimum_Temperature  + Deciduous  + 
                                              Evergreen +  Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName +
                                             (1|TrapSite/TrapID) + (1|TrapYear) + (1|Species) + (1|Ecoregion3),
                               family = nbinom2(), data=edrr_abundance_nn_bb)
   
   Anova(abundance_nn_bark)
   summary(abundance_nn_bark)
   
   
   
   summary_abd_nn_bb<- summary(abundance_nn_bark)
   coefs_m1_nn_bb <- as.data.frame(summary_abd_nn_bb$coefficients$cond)
   coefs_m1_nn_bb <- round(coefs_m1_nn_bb, 3)
   coefs_m1_nn_bb$Variable <- rownames(coefs_m1_nn_bb)
   write.csv(coefs_m1_nn_bb, "data/saved/Abundance~X_nonnativeAmbrosia.csv")
   
# ---- Model Diagnostics ----

   
   sim_res_m1_nn_bb <- simulateResiduals(abundance_nn_bark, plot = TRUE)
   testOutliers(sim_res_m1_nn_bb, type = "bootstrap")
  
   
# ---- Multi-model inference- Abundance ---- 
   
   options(na.action = "na.fail")
   dredge_abd_nn_bb <- dredge(abundance_nn_bark, rank = AIC, trace = TRUE)
   model_averaged_nn_bb<- model.avg(dredge_abd_nn_bb, subset=delta < 2) 
   modavg_nn_bb<- summary(model_averaged_nn_bb)
   df_nn_bb <-as.data.frame(modavg_nn_bb$coefmat.full)
   CI_df_nn_bb<- as.data.frame(confint(model_averaged_nn_bb, full=T))
   df_nn_bb$CI.min <-CI_df_nn_bb$`2.5 %` 
   df_nn_bb$CI.max <-CI_df_nn_bb$`97.5 %`
   df_nn_bb <- round(df_nn_bb, 3)
   df_nn_bb$Variable <- rownames(df_nn_bb)
   write.csv(df_nn_bb, "Data/Saved/Abundance_nonnativeBark_modavg.csv")
   

  saveRDS(dredge_abd_nn_bb, "Data/Saved/Abundance_nonnativeBark_dredge.rds")
   # Abundance_nonnativeBark_dredge<- readRDS("Data/Saved/Abundance_nonnativeBark_dredge.rds")
 #   
 #   
   
   
   

# ---- Plot ----  
   
   
   Abundance_nonnative_bark<-read.csv(file="Data/Saved/Abundance_nonnativeBark_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
   Abundance_nonnative_bark$Variable <- gsub("cond\\((.*)\\)", "\\1",Abundance_nonnative_bark$Variable)
   
   
   Abundance_nonnative_bark <- Abundance_nonnative_bark %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation", "Minimum_Temperature") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     TRUE ~ "Other"  ))
   Abundance_nonnative_bark <- Abundance_nonnative_bark[Abundance_nonnative_bark$Type != 'Other',]
   
   
   Abundance_nonnative_bark <- Abundance_nonnative_bark %>% filter(Variable != "(Intercept)")
   
   Abundance_nonnative_bark <- Abundance_nonnative_bark %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
   Abundance_nonnative_bark$Variable <- factor(Abundance_nonnative_bark$Variable,
                                               levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                                      "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
   
   
   Abundance_nonnative_bark <- Abundance_nonnative_bark %>%
     mutate(Variable = dplyr::recode(Variable,
                              "Annual_Precipitation" = "Annual Precipitation",
                              "Minimum_Temperature" = "Minimum Temperature",
                              "Deciduous" = "Deciduous Forest",
                              "Evergreen" = "Evergreen Forest",
                              "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                              "POPDEN" = "Human Population Density",
                              "Income_PerCapita" = "Per Capita Income",
                              "PortDist" = "Distance from Port"
     ))
   plot_nn_bark_abd <- ggplot(Abundance_nonnative_bark, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "Non-native Bark Beetles"
     ) +
     scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 2)) + 
     scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2")) +
     annotate("text",x = Inf,y = Inf,label = "G",hjust = 1.3,vjust = 1.3,fontface = "bold",size = 7) +
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major.y = element_blank(), 
       panel.grid.minor.y = element_blank(),
       scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
       plot.title = element_text(hjust = 0.5, face="bold", size=20),
       strip.text = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
       axis.text.x = element_text(face = "bold", size = 18),
       axis.text.y = element_blank(),
       axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
       axis.title.y = element_text(size = 22, face = "bold", vjust = 2),
       legend.position = "none",
     )+ guides(shape="none")
   
   
   plot_nn_bark_abd
   
   ggsave("Plots/Abundance_nonnative_bark.jpg", plot = plot_nn_bark_abd, width = 10, height = 12, dpi = 600)
# ---- Full Model- Richness----
 
   
   pred_vars_nn_bb <- edrr_richness_nn_bb[, c("Annual_Precipitation","Minimum_Temperature", "Evergreen", "Deciduous",
                                               "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist")]
   cor_matrix <- cor(pred_vars_nn_bb, use = "complete.obs", method = "pearson")
   correlation_matrix_nn_bb <- as.data.frame(cor_matrix)
   
   
   mean(edrr_richness_nn_bb$SpeciesCount)
   var(edrr_richness_nn_bb$SpeciesCount)
   hist(edrr_richness_nn_bb$SpeciesCount)
 richness_nn_bark<- glmmTMB(SpeciesCount ~ Annual_Precipitation +  Minimum_Temperature  + Deciduous + Evergreen + 
                              Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName + 
                              (1|TrapSite) + (1|TrapYear) + (1|Ecoregion3), family = compois(), data=edrr_richness_nn_bb)
 Anova(richness_nn_bark)  
 summary(richness_nn_bark)
   
   
   
 # richness_nn_bark_p<- glmmTMB(SpeciesCount ~ Annual_Precipitation +  Minimum_Temperature  + Deciduous + Evergreen + 
 #                              Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName + 
 #                              (1|TrapSite) + (1|TrapYear) + (1|Ecoregion3), family = poisson(), data=edrr_richness_nn_bb)
 # Anova(richness_nn_bark_p)  
 # summary(richness_nn_bark_p)
 # 
 #   AIC(richness_nn_bark, richness_nn_bark_p)
   
   
   summary_ric_nn_bb<- summary(richness_nn_bark)
   coefs_m2_nn_bb <- as.data.frame(summary_ric_nn_bb$coefficients$cond)
   coefs_m2_nn_bb <- round(coefs_m2_nn_bb, 3)
   coefs_m2_nn_bb$Variable <- rownames(coefs_m2_nn_bb)
   write.csv(coefs_m2_nn_bb, "data/saved/Richness~X_nonnativeBark.csv")
   
   
  
   
# ---- Model Diagnostics ----
   
   
   sim_res_m2_nn_bb <- simulateResiduals(richness_nn_bark, plot = TRUE)
   sim_res_m2_nn_bb_p <- simulateResiduals(richness_nn_bark_p, plot = TRUE)
   
   testOutliers(sim_res_m2_nn_bb, type = "bootstrap")
   
  
# ---- Multi-model inference- Richness ---- 
   
   options(na.action = "na.fail")
   dredge_ric_nn_bb <- dredge(richness_nn_bark, rank = AIC, trace = TRUE)
   model_averaged_nn2_bb<- model.avg(dredge_ric_nn_bb, subset=delta < 2) 
   modavg_nn2_bb<- summary(model_averaged_nn2_bb)
   df_nn2_bb <-as.data.frame(modavg_nn2_bb$coefmat.full)
   CI_df_nn2_bb<- as.data.frame(confint(model_averaged_nn2_bb, full=T))
   df_nn2_bb$CI.min <-CI_df_nn2_bb$`2.5 %` 
   df_nn2_bb$CI.max <-CI_df_nn2_bb$`97.5 %`
   df_nn2_bb <- round(df_nn2_bb, 3)
   df_nn2_bb$Variable <- rownames(df_nn2_bb)
   write.csv(df_nn2_bb, "Data/Saved/Richness_nonnativeBark_modavg.csv")
   
   attr(dredge_ric_nn_bb, "warnings")
   top_models <- subset(dredge_ric_nn_bb, delta < 2)
   summary(top_models)

   saveRDS(dredge_ric_nn_bb, "Data/Saved/Richness_nonnativeBark_dredge.rds")
   #dredge_ric_nn_bb<-readRDS(file="Data/Saved/Richness_nonnativeBark_dredge.rds")
   
   
   
   
   
# ---- Plot ----  
   
   
   Richness_nonnative_bark<-read.csv(file="Data/Saved/Richness_nonnativeBark_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
   Richness_nonnative_bark$Variable <- gsub("cond\\((.*)\\)", "\\1",Richness_nonnative_bark$Variable)
   
   
   Richness_nonnative_bark <- Richness_nonnative_bark %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation", "Minimum_Temperature") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     TRUE ~ "Other"  ))
   Richness_nonnative_bark <- Richness_nonnative_bark[Richness_nonnative_bark$Type != 'Other',]
   
   
   Richness_nonnative_bark <- Richness_nonnative_bark %>% filter(Variable != "(Intercept)")
   
   Richness_nonnative_bark <- Richness_nonnative_bark %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
   Richness_nonnative_bark$Variable <- factor(Richness_nonnative_bark$Variable, 
                                              levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                                          "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
   
   
   Richness_nonnative_bark <- Richness_nonnative_bark %>%
     mutate(Variable = dplyr::recode(Variable,
                              "Annual_Precipitation" = "Annual Precipitation",
                              "Minimum_Temperature" = "Minimum Temperature",
                              "Deciduous" = "Deciduous Forest",
                              "Evergreen" = "Evergreen Forest",
                              "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                              "POPDEN" = "Human Population Density",
                              "Income_PerCapita" = "Per Capita Income",
                              "PortDist" = "Distance from Port"
     ))
   plot_nn_bark_ric <- ggplot(Richness_nonnative_bark, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "Non-native Bark Beetles"
     ) +
     scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 2)) + 
     scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2")) +
     annotate("text",x = Inf,y = Inf,label = "G",hjust = 1.3,vjust = 1.3,fontface = "bold",size = 7) +
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major.y = element_blank(), 
       panel.grid.minor.y = element_blank(),
       scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
       plot.title = element_text(hjust = 0.5, face="bold", size=20),
       strip.text = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
       axis.text.x = element_text(face = "bold", size = 18),
       axis.text.y = element_blank(),
       axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
       axis.title.y = element_text(size = 22, face = "bold", vjust = 2),
       legend.position = "none",
     )+ guides(shape="none")
   
   plot_nn_bark_ric
   
   ggsave("Plots/Richness_nonnative_bark.jpg", plot = plot_nn_bark_ric, width = 10, height = 12, dpi = 600)
   
##################################
###############################
# ---- Native Aggregated ----    
   
   edrr_long_na <- edrr_native %>%
     group_by(TrapID, Species) %>%
     summarise(count = sum(PestCount, na.rm = TRUE), .groups = "drop")
   
   # wide format
   edrr_wide_na <- edrr_long_na %>%
     pivot_wider(
       names_from = Species,
       values_from = count,
       values_fill = 0  # Replace NA with 0 
     )
   
   # Extract metadata , keeps first value for each variable
   meta_info_na <- edrr_native %>%
     group_by(TrapID) %>%
     summarise(across(c(Cell_ID, SiteLure, StateCode, CountyName,SiteTown,SiteLoc,SiteType, TrapYear, SiteNumber, LureName, FIPS, Longitude, Latitude,
                        TrapSite, landcover, LC_class,  Ecoregion3,  
                        Annual_Precipitation , Minimum_Temperature,  Deciduous, Evergreen,  
                        Shannon_Div_Landscape, POPDEN,  Income_PerCapita, PortDist), first), .groups = "drop")
   
   
   edrr_abundance_na <- edrr_long_na %>% left_join(meta_info_na, by = "TrapID")
  
   edrr_richness_na <- edrr_wide_na %>% left_join(meta_info_na, by = "TrapID")
   
   
   
 
   
   
   #Calculating Species Count for each Trap ID
   edrr_richness_na$SpeciesCount<-rowSums(edrr_richness_na[, 2:234] !=0)
   write.csv(edrr_abundance_na, "Data/Saved/edrr_abundance_na.csv", row.names = FALSE)
   write.csv(edrr_richness_na, "Data/Saved/edrr_richness_na.csv", row.names = FALSE)
  
   
   
# ---- Checking correlation among predictors ----   
   pred_vars_na <- edrr_abundance_na[, c("Annual_Precipitation","Minimum_Temperature", "Evergreen", "Deciduous",
                                      "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist")]
   cor_matrix_na <- cor(pred_vars_na, use = "complete.obs", method = "pearson")
   correlation_matrix_na <- as.data.frame(cor_matrix_na)
   
   
   high_cor_indices <- which(abs(correlation_matrix_na) >= 0.7 & abs(correlation_matrix_na) < 1, arr.ind = TRUE)
   high_cor_pairs <- data.frame(
     Variable1 = rownames(correlation_matrix_na)[high_cor_indices[, 1]],
     Variable2 = colnames(correlation_matrix_na)[high_cor_indices[, 2]],
     Correlation = correlation_matrix_na[high_cor_indices]
   )
   print(high_cor_pairs)
   

   
# ---- Full Model Abundance ----
   
   
  
   m1_na<- glmmTMB(count ~ Annual_Precipitation +  Minimum_Temperature  + Deciduous + Evergreen + 
                     Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName + 
                     (1|TrapSite/TrapID) + (1|TrapYear) + (1|Species) + (1|Ecoregion3) , 
                   family = nbinom2(), data=edrr_abundance_na)
   
   Anova(m1_na)
   summary(m1_na)
   summarym1_na<- summary(m1_na)
   coefs_m1_na <- as.data.frame(summarym1_na$coefficients$cond)
   coefs_m1_na<- round(coefs_m1_na, 3)
   coefs_m1_na$Variable <- rownames(coefs_m1_na)
   write.csv(coefs_m1_na, "data/saved/Abundance~X_nonnativeEDRR.csv")
   
# ---- Model Diagnostics ----
   
   
   sim_res_m1_na <- simulateResiduals(m1_na, plot = TRUE)
   testOutliers(sim_res_m1_na, type = "bootstrap")
   
   
# ---- Multi-model inference- Abundance ---- 
   
   options(na.action = "na.fail")
   dredge_m1_na <- dredge(m1_na, rank = AIC, trace = TRUE)
   
   model_averaged_na<- model.avg(dredge_m1_na, subset=delta < 2) 
   modavg_na<- summary(model_averaged_na)
   df_na <-as.data.frame(modavg_na$coefmat.full)
   CI_df_na<- as.data.frame(confint(model_averaged_na, full=T))
   df_na$CI.min <-CI_df_na$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
   df_na$CI.max <-CI_df_na$`97.5 %`
   df_na <- round(df_na, 3)
   df_na$Variable <- rownames(df_na)
   write.csv(df_na, "Data/Saved/Abundance_native_modavg.csv")
   
   saveRDS(dredge_m1_na, file = "Data/Saved/dredge_m1_na.rds")
   # dredge_m1_na <- readRDS("Data/Saved/dredge_m1_na.rds")
   
# ---- Plot ----  
   
   
   Abundance_native<-read.csv(file="Data/Saved/Abundance_native_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
   Abundance_native$Variable <- gsub("cond\\((.*)\\)", "\\1",Abundance_native$Variable)
   
   
   Abundance_native <- Abundance_native %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation", "Minimum_Temperature") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     TRUE ~ "Other"  ))
   Abundance_native <- Abundance_native[Abundance_native$Type != 'Other',]
   
   
   Abundance_native <- Abundance_native %>% filter(Variable != "(Intercept)")
   
   Abundance_native <- Abundance_native %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
   Abundance_native$Variable <- factor(Abundance_native$Variable, 
                                       levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                                                                   "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
   
   
   Abundance_native <- Abundance_native %>%
     mutate(Variable = dplyr::recode(Variable,
                              "Annual_Precipitation" = "Annual Precipitation",
                              "Minimum_Temperature" = "Minimum Temperature",
                              "Deciduous" = "Deciduous Forest",
                              "Evergreen" = "Evergreen Forest",
                              "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                              "POPDEN" = "Human Population Density",
                              "Income_PerCapita" = "Per Capita Income",
                              "PortDist" = "Distance from Port"
     ))
   
   plot_na_abd <- ggplot(Abundance_native, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "Native Scolytines"
     ) +
     scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) + 
     scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2")) +
     annotate("text",x = Inf,y = Inf,label = "B",hjust = 1.3,vjust = 1.3,fontface = "bold",size = 7) +
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major.y = element_blank(), 
       panel.grid.minor.y = element_blank(),
       scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
       plot.title = element_text(hjust = 0.5, face="bold", size=20),
       strip.text = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
       axis.text.x = element_text(face = "bold", size = 18),
       axis.text.y = element_text(face = "bold", size = 18),
       axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
       axis.title.y = element_text(size = 22, face = "bold", vjust = 2),
       legend.position = "none",
     )+ guides(shape="none")
   
   
   plot_na_abd
   
   ggsave("Plots/Abundance_native.jpg", plot = plot_na_abd, width = 10, height = 12, dpi = 600)
   
# ---- Full Model Richness ---- 
   pred_vars_na <- edrr_richness_na[, c("Annual_Precipitation","Minimum_Temperature", "Evergreen", "Deciduous",
                                      "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist")]
   cor_matrix <- cor(pred_vars_na, use = "complete.obs", method = "pearson")
   correlation_matrix_na <- as.data.frame(cor_matrix)
   
   
  
   m2_na<- glmmTMB(SpeciesCount ~ Annual_Precipitation +  Minimum_Temperature  + Deciduous  + 
                     Evergreen + Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName +
                     (1|TrapSite) + (1|TrapYear) + (1|Ecoregion3), family = poisson(), data=edrr_richness_na)
   Anova(m2_na)
   summary(m2_na)
   
  
   
   
   summarym2_na<- summary(m2_na)
   coefs_m2_na <- as.data.frame(summarym2_na$coefficients$cond)
   coefs_m2_na <- round(coefs_m2_na, 3)
   coefs_m2_na$Variable <- rownames(coefs_m2_na)
   write.csv(coefs_m2_na, "data/saved/Richness~X_nativeEDRR.csv")
   
# ---- Model Diagnostics ----   
   sim_res_m2_na <- simulateResiduals(m2_na, plot = TRUE)
   testOutliers(sim_res_m2_na, type = "bootstrap")
   
  
# ---- Multi-model inference- Richness ----    
   options(na.action = "na.fail")
   dredge_m2_na <- dredge(m2_na, rank = AIC, trace = TRUE)
   model_averaged_na2<- model.avg(dredge_m2_na, subset=delta < 2) 
   modavg_na2<- summary(model_averaged_na2)
   df_na2 <-as.data.frame(modavg_na2$coefmat.full)
   CI_df_na2<- as.data.frame(confint(model_averaged_na2, full=T))
   df_na2$CI.min <-CI_df_na2$`2.5 %` 
   df_na2$CI.max <-CI_df_na2$`97.5 %`
   df_na2 <- round(df_na2, 3)
   df_na2$Variable <- rownames(df_na2)
   write.csv(df_na2, "Data/Saved/Richness_native_modavg.csv")
   
  
   saveRDS(dredge_m2_na, "Data/Saved/Richness_native_dredge.rds")
  
   
   
   
   
   
   
# ---- Plot ----  
   
   
   Richness_native<-read.csv(file="Data/Saved/Richness_native_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
   Richness_native$Variable <- gsub("cond\\((.*)\\)", "\\1",Richness_native$Variable)
   
   
   Richness_native <- Richness_native %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation", "Minimum_Temperature") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     TRUE ~ "Other"  ))
   Richness_native <- Richness_native[Richness_native$Type != 'Other',]
   
   
   Richness_native <- Richness_native %>% filter(Variable != "(Intercept)")
   
   Richness_native <- Richness_native %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
   Richness_native$Variable <- factor(Richness_native$Variable, 
                                      levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                              "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
   
   
   Richness_native <- Richness_native %>%
     mutate(Variable = dplyr::recode(Variable,
                              "Annual_Precipitation" = "Annual Precipitation",
                              "Minimum_Temperature" = "Minimum Temperature",
                              "Deciduous" = "Deciduous Forest",
                              "Evergreen" = "Evergreen Forest",
                              "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                              "POPDEN" = "Human Population Density",
                              "Income_PerCapita" = "Per Capita Income",
                              "PortDist" = "Distance from Port"
     ))
   plot_na_ric <- ggplot(Richness_native, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "Native Scolytines"
     ) +
     scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) + 
     scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2")) +
     annotate("text",x = Inf,y = Inf,label = "B",hjust = 1.3,vjust = 1.3,fontface = "bold",size = 7) +
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major.y = element_blank(), 
       panel.grid.minor.y = element_blank(),
       scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
       plot.title = element_text(hjust = 0.5, face="bold", size=20),
       strip.text = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
       axis.text.x = element_text(face = "bold", size = 18),
       axis.text.y = element_text(face = "bold", size = 18),
       axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
       axis.title.y = element_text(size = 22, face = "bold", vjust = 2),
       legend.position = "none",
     )+ guides(shape="none")
   
   
   plot_na_ric
   
   ggsave("Plots/Richness_native.jpg", plot = plot_na_ric, width = 10, height = 12, dpi = 600)
   
############################################
#############################################
# ---- Native Abrosia Beetles ----
   
  
   
   # Summarize counts
   edrr_long_na_ab <- edrr_na_ab %>%
     group_by(TrapID, Species) %>%
     summarise(count = sum(PestCount, na.rm = TRUE), .groups = "drop")
   
   # wide format
   edrr_wide_na_ab <- edrr_long_na_ab %>%
     pivot_wider(
       names_from = Species,
       values_from = count,
       values_fill = 0  # Replace NA with 0 
     )
   
   # Extract metadata , keeps first value for each variable
   meta_info_na_ab <- edrr_na_ab %>%
     group_by(TrapID) %>%
     summarise(across(c(Cell_ID, SiteLure, StateCode, CountyName,SiteTown,SiteLoc,SiteType, TrapYear, SiteNumber, LureName, FIPS, Longitude, Latitude,
                        TrapSite, landcover, LC_class,  Ecoregion3,  
                        Annual_Precipitation , Minimum_Temperature,  Deciduous, Evergreen,  
                        Shannon_Div_Landscape, POPDEN,  Income_PerCapita, PortDist), first), .groups = "drop")
   
   
   
   edrr_abundance_na_ab <- edrr_long_na_ab %>% left_join(meta_info_na_ab, by = "TrapID")

   edrr_richness_na_ab <- edrr_wide_na_ab %>% left_join(meta_info_na_ab, by = "TrapID")
   
  
   
   
   #Calculating Species Count for each Trap ID
   edrr_richness_na_ab$SpeciesCount<-rowSums(edrr_richness_na_ab[, 2:40] !=0)
   
  
   write.csv(edrr_abundance_na_ab, "Data/Saved/edrr_abundance_na_ab.csv", row.names = FALSE)
   write.csv(edrr_richness_na_ab, "Data/Saved/edrr_richness_na_ab.csv", row.names = FALSE)
   
   
  
   
   
# ---- Checking correlation among predictors ----   
   pred_vars_na_ab <- edrr_abundance_na_ab[, c("Annual_Precipitation","Minimum_Temperature", "Evergreen", "Deciduous",
                                            "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist")]
   cor_matrix_na_ab <- cor(pred_vars_na_ab, use = "complete.obs", method = "pearson")
   correlation_matrix_na_ab <- as.data.frame(cor_matrix_na_ab)
   # gg_correlation <- ggcorr(pred_vars_na_ab, label = TRUE)
   # gg_correlation
   
   
   high_cor_indices <- which(abs(correlation_matrix_na_ab) >= 0.7 & abs(correlation_matrix_na_ab) < 1, arr.ind = TRUE)
   high_cor_pairs <- data.frame(
     Variable1 = rownames(correlation_matrix_na_ab)[high_cor_indices[, 1]],
     Variable2 = colnames(correlation_matrix_na_ab)[high_cor_indices[, 2]],
     Correlation = correlation_matrix_na_ab[high_cor_indices]
   )
   print(high_cor_pairs)
   
# ---- Individual models ----
   mintemp<- glmmTMB(count ~  Minimum_Temperature + Annual_Precipitation + Evergreen + 
                       Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName + (1|TrapSite/TrapID) + (1|Species)
                     + (1|TrapYear) +(1|Ecoregion3), family = nbinom2(), data=edrr_abundance_na_ab)
   
   Anova(mintemp)
   summary(mintemp)
   
   
   sim_res_mintemp <- simulateResiduals(mintemp, plot = TRUE)
   testDispersion(sim_res_mintemp)
   testOutliers(sim_res_mintemp, type = "bootstrap")
   
   
   decid<- glmmTMB(count ~  Deciduous + Annual_Precipitation + Evergreen + 
                     Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName + (1|TrapSite/TrapID) + (1|Species)
                   + (1|TrapYear) +(1|Ecoregion3), family = nbinom2(), data=edrr_abundance_na_ab)
   
   Anova(decid)
   summary(decid)
   
   
   AIC(mintemp, decid)
   
   
   fml_mintemp <- Reduce(paste, deparse(formula(mintemp)))
   fml_decid   <- Reduce(paste, deparse(formula(decid)))
   aic_vals <- AIC(mintemp, decid)$AIC
   
   
   modelcomp4 <- data.frame(
     Group = c("Native Ambrosia Beetle Abundance", "Native Ambrosia Beetle Abundance"),
     Model = c(
       gsub("count ~ ","", fml_mintemp),
       gsub("count ~ ","", fml_decid)
     ),
     AIC = round(aic_vals, 2),
     Delta_AIC = round(aic_vals - min(aic_vals), 2)
   )
   modelcomp4 <- modelcomp4[order(modelcomp4$AIC), ]
   modelcomp4$Group[duplicated(modelcomp4$Group)] <- ""
   
   names(modelcomp4)[names(modelcomp4) == "Delta_AIC"] <- "ΔAIC"
   
   
# ---- Full Model- Abundance ----
   
   abundance_na_ambrosia<- glmmTMB(count ~ Annual_Precipitation + Deciduous  + 
                                     Evergreen + Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName +
                                     (1|TrapSite/TrapID) + (1|TrapYear) +(1|Species) + (1|Ecoregion3),
                                   family = nbinom2(), data=edrr_abundance_na_ab)
   
   Anova(abundance_na_ambrosia)
   summary(abundance_na_ambrosia)
   
   
   
   summary_abd_na_ab<- summary(abundance_na_ambrosia)
   coefs_m1_na_ab <- as.data.frame(summary_abd_na_ab$coefficients$cond)
   coefs_m1_na_ab <- round(coefs_m1_na_ab, 3)
   coefs_m1_na_ab$Variable <- rownames(coefs_m1_na_ab)
   write.csv(coefs_m1_na_ab, "data/saved/Abundance~X_nativeAmbrosia.csv")
   
# ---- Model Diagnostics ----
  
   
   sim_res_m1_na_ab <- simulateResiduals(abundance_na_ambrosia, plot = TRUE)
   testOutliers(sim_res_m1_na_ab, type = "bootstrap")

   
# ---- Multi-model inference- Abundance ---- 
   
   options(na.action = "na.fail")
   dredge_abd_na_ab <- dredge(abundance_na_ambrosia, rank = AIC, trace = TRUE)
   model_averaged_na_ab<- model.avg(dredge_abd_na_ab, subset=delta < 2) 
   modavg_na_ab<- summary(model_averaged_na_ab)
   df_na_ab <-as.data.frame(modavg_na_ab$coefmat.full)
   CI_df_na_ab<- as.data.frame(confint(model_averaged_na_ab, full=T))
   df_na_ab$CI.min <-CI_df_na_ab$`2.5 %` 
   df_na_ab$CI.max <-CI_df_na_ab$`97.5 %`
   df_na_ab <- round(df_na_ab, 3)
   df_na_ab$Variable <- rownames(df_na_ab)
   write.csv(df_na_ab, "Data/Saved/Abundance_nativeAmbrosia_modavg.csv")
   
 
   saveRDS(dredge_abd_na_ab, "Data/Saved/Abundance_nativeAmbrosia_dredge.rds")
  
   
   
   
   
# ---- Plot ----  
   
   
   Abundance_native_ambrosia<-read.csv(file="Data/Saved/Abundance_nativeAmbrosia_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
   Abundance_native_ambrosia$Variable <- gsub("cond\\((.*)\\)", "\\1",Abundance_native_ambrosia$Variable)
   Abundance_native_ambrosia <- Abundance_native_ambrosia %>%
     bind_rows(
       tibble(
         Variable = "Minimum_Temperature",
         Estimate = 0,
         CI.min = 0,
         CI.max = 0,
         Pr...z.. = NA,
         Type = "Removed",
         p_value = NA,
         Significant = NA
       )
     )
   
   Abundance_native_ambrosia <- Abundance_native_ambrosia %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     Variable %in% c("Minimum_Temperature") ~ "Removed", 
     TRUE ~ "Other"  ))
   Abundance_native_ambrosia <- Abundance_native_ambrosia[Abundance_native_ambrosia$Type != 'Other',]
   
   
   Abundance_native_ambrosia <- Abundance_native_ambrosia %>% filter(Variable != "(Intercept)")
   
   Abundance_native_ambrosia <- Abundance_native_ambrosia %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
   Abundance_native_ambrosia$Variable <- factor(Abundance_native_ambrosia$Variable, 
                                                levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                                              "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
   
   
   Abundance_native_ambrosia <- Abundance_native_ambrosia %>%
     mutate(Variable = dplyr::recode(Variable,
                              "Annual_Precipitation" = "Annual Precipitation",
                              "Minimum_Temperature" = "Minimum Temperature",
                              "Deciduous" = "Deciduous Forest",
                              "Evergreen" = "Evergreen Forest",
                              "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                              "POPDEN" = "Human Population Density",
                              "Income_PerCapita" = "Per Capita Income",
                              "PortDist" = "Distance from Port"
     ))
   plot_na_ambrosia_abd <- ggplot(Abundance_native_ambrosia, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "Native Ambrosia Beetles"
     ) +
     scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), na.value = 4) +
     scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2", "Removed" = "gray50"),
                        breaks = c("Climatic", "Landscape", "Socioeconomic") ) +
     annotate("text",x = Inf,y = Inf,label = "D",hjust = 1.3,vjust = 1.3,fontface = "bold",size = 7) +
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major.y = element_blank(), 
       panel.grid.minor.y = element_blank(),
       scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
       plot.title = element_text(hjust = 0.5, face="bold", size=20),
       strip.text = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
       axis.text.x = element_text(face = "bold", size = 18),
       axis.text.y = element_text(face = "bold", size = 18),
       axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
       axis.title.y = element_text(size = 22, face = "bold", vjust = 2),
       legend.position = "none",
     )+ guides(shape="none")
   
   plot_na_ambrosia_abd
   
   ggsave("Plots/Abundance_native_ambrosia.jpg", plot = plot_na_ambrosia_abd, width = 10, height = 12, dpi = 600)
   
# ---- Full Model- Richness ----
   pred_vars_na_ab <- edrr_richness_na_ab[, c("Annual_Precipitation","Minimum_Temperature", "Evergreen", "Deciduous",
                                               "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist")]
   cor_matrix <- cor(pred_vars_na_ab, use = "complete.obs", method = "pearson")
   correlation_matrix_na_ab <- as.data.frame(cor_matrix)
   
   hist(edrr_richness_na_ab$SpeciesCount)
   mean(edrr_richness_na_ab$SpeciesCount)
   var(edrr_richness_na_ab$SpeciesCount)
   
   
  richness_na_ambrosia<- glmmTMB(SpeciesCount ~ Annual_Precipitation +  Minimum_Temperature  + Deciduous  + 
                                   Evergreen + Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName +
                                   (1|TrapSite) + (1|TrapYear) + (1|Ecoregion3), family = poisson(), data=edrr_richness_na_ab)
   
   Anova(richness_na_ambrosia)
   summary(richness_na_ambrosia)
   
   
   
   summary_ric_na_ab<- summary(richness_na_ambrosia)
   coefs_m2_na_ab <- as.data.frame(summary_ric_na_ab$coefficients$cond)
   coefs_m2_na_ab <- round(coefs_m2_na_ab, 3)
   coefs_m2_na_ab$Variable <- rownames(coefs_m2_na_ab)
  
# ---- Model Diagnostics ----
   
   
   sim_res_m2_na_ab <- simulateResiduals(richness_na_ambrosia, plot = TRUE)
   testOutliers(sim_res_m2_na_ab, type = "bootstrap")
   
# ---- Multi-model inference- Richness ---- 
   
   options(na.action = "na.fail")
   dredge_ric_na_ab <- dredge(richness_na_ambrosia, rank = AIC, trace = TRUE)
   model_averaged_na2_ab<- model.avg(dredge_ric_na_ab, subset=delta < 2) 
   modavg_na2_ab<- summary(model_averaged_na2_ab)
   df_na2_ab <-as.data.frame(modavg_na2_ab$coefmat.full)
   CI_df_na2_ab<- as.data.frame(confint(model_averaged_na2_ab, full=T))
   df_na2_ab$CI.min <-CI_df_na2_ab$`2.5 %` 
   df_na2_ab$CI.max <-CI_df_na2_ab$`97.5 %`
   df_na2_ab <- round(df_na2_ab, 3)
   df_na2_ab$Variable <- rownames(df_na2_ab)
   write.csv(df_na2_ab, "Data/Saved/Richness_nativeAmbrosia_modavg.csv")
   
  
   saveRDS(dredge_ric_na_ab, "Data/Saved/Richness_nativeAmbrosia_dredge.rds")
  
  # dredge_ric_na_ab <-  readRDS("Data/Saved/Richness_nativeAmbrosia_dredge.csv")
   
   
# ---- Plot ----  
   
   
   Richness_native_ambrosia<-read.csv(file="Data/Saved/Richness_nativeAmbrosia_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
   Richness_native_ambrosia$Variable <- gsub("cond\\((.*)\\)", "\\1",Richness_native_ambrosia$Variable)
   
   
   Richness_native_ambrosia <- Richness_native_ambrosia %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation", "Minimum_Temperature") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     TRUE ~ "Other"  ))
   Richness_native_ambrosia <- Richness_native_ambrosia[Richness_native_ambrosia$Type != 'Other',]
   
   
   Richness_native_ambrosia <- Richness_native_ambrosia %>% filter(Variable != "(Intercept)")
   
   Richness_native_ambrosia <- Richness_native_ambrosia %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
   Richness_native_ambrosia$Variable <- factor(Richness_native_ambrosia$Variable, 
                                               levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                                            "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
   
   
   Richness_native_ambrosia <- Richness_native_ambrosia %>%
     mutate(Variable = dplyr::recode(Variable,
                              "Annual_Precipitation" = "Annual Precipitation",
                              "Minimum_Temperature" = "Minimum Temperature",
                              "Deciduous" = "Deciduous Forest",
                              "Evergreen" = "Evergreen Forest",
                              "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                              "POPDEN" = "Human Population Density",
                              "Income_PerCapita" = "Per Capita Income",
                              "PortDist" = "Distance from Port"
     ))
   plot_na_ambrosia_ric <- ggplot(Richness_native_ambrosia, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "Native Ambrosia Beetles"
     ) +
     scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) + 
     scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2")) +
     annotate("text",x = Inf,y = Inf,label = "D",hjust = 1.3,vjust = 1.3,fontface = "bold",size = 7) +
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major.y = element_blank(), 
       panel.grid.minor.y = element_blank(),
       scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
       plot.title = element_text(hjust = 0.5, face="bold", size=20),
       strip.text = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
       axis.text.x = element_text(face = "bold", size = 18),
       axis.text.y = element_text(face = "bold", size = 18),
       axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
       axis.title.y = element_text(size = 22, face = "bold", vjust = 2),
       legend.position = "none",
     )+ guides(shape="none")
   
   
   plot_na_ambrosia_ric
   
   ggsave("Plots/Richness_native_ambrosia.jpg", plot = plot_na_ambrosia_ric, width = 10, height = 12, dpi = 600)
   
##################################################################
#######################################################################
# ---- Native Bark Beetles ----  
   
   # Summarize counts
   edrr_long_na_bb <- edrr_na_bb %>%
     group_by(TrapID, Species) %>%
     summarise(count = sum(PestCount, na.rm = TRUE), .groups = "drop")
   
   # wide format
   edrr_wide_na_bb <- edrr_long_na_bb %>%
     pivot_wider(
       names_from = Species,
       values_from = count,
       values_fill = 0  # Replace NA with 0 
     )
   
   # Extract metadata , keeps first value for each variable
   meta_info_na_bb <- edrr_na_bb %>%
     group_by(TrapID) %>%
     summarise(across(c(Cell_ID, SiteLure, StateCode, CountyName,SiteTown,SiteLoc,SiteType, TrapYear, SiteNumber, LureName, FIPS, Longitude, Latitude,
                        TrapSite, landcover, LC_class,  Ecoregion3,  
                        Annual_Precipitation , Minimum_Temperature,  Deciduous, Evergreen,  
                        Shannon_Div_Landscape, POPDEN,  Income_PerCapita, PortDist), first), .groups = "drop")
   
   
   
   edrr_abundance_na_bb <- edrr_long_na_bb %>% left_join(meta_info_na_bb, by = "TrapID")
   edrr_richness_na_bb <- edrr_wide_na_bb %>% left_join(meta_info_na_bb, by = "TrapID")
   
   
   #Calculating Species Count for each Trap ID
   edrr_richness_na_bb$SpeciesCount<-rowSums(edrr_richness_na_bb[, 2:132] !=0)
   
   write.csv(edrr_abundance_na_bb, "Data/Saved/edrr_abundance_na_bb.csv", row.names = FALSE)
   write.csv(edrr_richness_na_bb, "Data/Saved/edrr_richness_na_bb.csv", row.names = FALSE)

   
# ---- Checking correlation among predictors ----   
   pred_vars_na_bb <- edrr_abundance_na_bb[, c("Annual_Precipitation","Minimum_Temperature", "Evergreen", "Deciduous",
                                            "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist")]
   cor_matrix_na_bb <- cor(pred_vars_na_bb, use = "complete.obs", method = "pearson")
   correlation_matrix_na_bb <- as.data.frame(cor_matrix_na_bb)
   
   
   
   high_cor_indices <- which(abs(correlation_matrix_na_bb) >= 0.7 & abs(correlation_matrix_na_bb) < 1, arr.ind = TRUE)
   high_cor_pairs <- data.frame(
     Variable1 = rownames(correlation_matrix_na_bb)[high_cor_indices[, 1]],
     Variable2 = colnames(correlation_matrix_na_bb)[high_cor_indices[, 2]],
     Correlation = correlation_matrix_na_bb[high_cor_indices]
   )
   print(high_cor_pairs)
   

   
# ---- Full Model- Abundance ----
   
   abundance_na_bark<- glmmTMB(count ~ Annual_Precipitation +  Minimum_Temperature  + Deciduous  + 
                                     Evergreen  + Shannon_Div_Landscape + POPDEN + Income_PerCapita + PortDist + LureName+
                                     (1|TrapSite/TrapID) + (1|TrapYear) +(1|Species) + (1|Ecoregion3), 
                               family = nbinom2(), data=edrr_abundance_na_bb)
   
  Anova(abundance_na_bark)
    summary(abundance_na_bark)
   
   
   
   summary_abd_na_bb<- summary(abundance_na_bark)
   coefs_m1_na_bb <- as.data.frame(summary_abd_na_bb$coefficients$cond)
   coefs_m1_na_bb <- round(coefs_m1_na_bb, 3)
   coefs_m1_na_bb$Variable <- rownames(coefs_m1_na_bb)
  
# ---- Model Diagnostics ----
   
   
   sim_res_m1_na_bb <- simulateResiduals(abundance_na_bark, plot = TRUE)
   testOutliers(sim_res_m1_na_bb, type = "bootstrap")
   
   
# ---- Multi-model inference- Abundance ---- 
   
   options(na.action = "na.fail")
   dredge_abd_na_bb <- dredge(abundance_na_bark, rank = AIC, trace = TRUE)
   model_averaged_na_bb<- model.avg(dredge_abd_na_bb, subset=delta < 2) 
   modavg_na_bb<- summary(model_averaged_na_bb)
   df_na_bb <-as.data.frame(modavg_na_bb$coefmat.full)
   CI_df_na_bb<- as.data.frame(confint(model_averaged_na_bb, full=T))
   df_na_bb$CI.min <-CI_df_na_bb$`2.5 %` 
   df_na_bb$CI.max <-CI_df_na_bb$`97.5 %`
   df_na_bb <- round(df_na_bb, 3)
   df_na_bb$Variable <- rownames(df_na_bb)
   write.csv(df_na_bb, "Data/Saved/Abundance_nativeBark_modavg.csv")
   
  
   saveRDS(dredge_abd_na_bb, "Data/Saved/Abundance_nativeBark_dredge.rds")
  
   
# ---- Plot ----  
   
   
   Abundance_native_bark<-read.csv(file="Data/Saved/Abundance_nativeBark_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
   Abundance_native_bark$Variable <- gsub("cond\\((.*)\\)", "\\1",Abundance_native_bark$Variable)
   
   
   Abundance_native_bark <- Abundance_native_bark %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation", "Minimum_Temperature") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     TRUE ~ "Other"  ))
   Abundance_native_bark <- Abundance_native_bark[Abundance_native_bark$Type != 'Other',]
   
   
   Abundance_native_bark <- Abundance_native_bark %>% filter(Variable != "(Intercept)")
   
   Abundance_native_bark <- Abundance_native_bark %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
   Abundance_native_bark$Variable <- factor(Abundance_native_bark$Variable, 
                                            levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                                                                             "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
   
   
   Abundance_native_bark <- Abundance_native_bark %>%
     mutate(Variable = dplyr::recode(Variable,
                              "Annual_Precipitation" = "Annual Precipitation",
                              "Minimum_Temperature" = "Minimum Temperature",
                              "Deciduous" = "Deciduous Forest",
                              "Evergreen" = "Evergreen Forest",
                              "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                              "POPDEN" = "Human Population Density",
                              "Income_PerCapita" = "Per Capita Income",
                              "PortDist" = "Distance from Port"
     ))
   plot_na_bark_abd <- ggplot(Abundance_native_bark, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "Native Bark Beetles"
     ) +
     scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) + 
     scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2")) +
     annotate("text",x = Inf,y = Inf,label = "F",hjust = 1.3,vjust = 1.3,fontface = "bold",size = 7) +
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major.y = element_blank(), 
       panel.grid.minor.y = element_blank(),
       scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
       plot.title = element_text(hjust = 0.5, face="bold", size=20),
       strip.text = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
       axis.text.x = element_text(face = "bold", size = 18),
       axis.text.y = element_text(face = "bold", size = 18),
       axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
       axis.title.y = element_text(size = 22, face = "bold", vjust = 2),
       legend.position = "none",
     )+ guides(shape="none")
   
   plot_na_bark_abd
   
   ggsave("Plots/Abundance_native_bark.jpg", plot = plot_na_bark_abd, width = 10, height = 12, dpi = 600)
   
# ---- Full Model- Richness ----
   pred_vars_na_bb <- edrr_richness_na_bb[, c("Annual_Precipitation","Minimum_Temperature", "Evergreen", "Deciduous",
                                               "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist")]
   cor_matrix <- cor(pred_vars_na_bb, use = "complete.obs", method = "pearson")
   correlation_matrix_na_bb <- as.data.frame(cor_matrix)
   
   hist(edrr_richness_na_bb$SpeciesCount)
   mean(edrr_richness_na_bb$SpeciesCount)
   var(edrr_richness_na_bb$SpeciesCount)
   
   
   richness_na_bark<- glmmTMB(SpeciesCount ~ Annual_Precipitation +  Minimum_Temperature  + Deciduous  + 
                                    Evergreen + Shannon_Div_Landscape +  POPDEN + Income_PerCapita + PortDist + LureName +
                                    (1|TrapSite) + (1|TrapYear) + (1|Ecoregion3), 
                              family = poisson(), data=edrr_richness_na_bb)
   
   Anova(richness_na_bark)
   summary(richness_na_bark)
   
   
   
    summary_ric_na_bb<- summary(richness_na_bark)
   coefs_m2_na_bb <- as.data.frame(summary_ric_na_bb$coefficients$cond)
   coefs_m2_na_bb <- round( coefs_m2_na_bb, 3)
   coefs_m2_na_bb$Variable <- rownames( coefs_m2_na_bb)
   # write.csv(coefs_m1_nn_ab, "data/saved/Abundance~X_nonnativeAmbrosia.csv")
   
   
# ---- Model Diagnostics ----
   
   
   sim_res_m2_na_bb <- simulateResiduals(richness_na_bark, plot = TRUE)
   testOutliers(sim_res_m2_na_bb, type = "bootstrap")
   
# ---- Multi-model inference- Richness ---- 
   
   options(na.action = "na.fail")
   dredge_ric_na_bb <- dredge(richness_na_bark, rank = AIC, trace = TRUE)
   model_averaged_na2_bb<- model.avg(dredge_ric_na_bb, subset=delta < 2) 
   modavg_na2_bb<- summary(model_averaged_na2_bb)
   df_na2_bb <-as.data.frame(modavg_na2_bb$coefmat.full)
   CI_df_na2_bb<- as.data.frame(confint(model_averaged_na2_bb, full=T))
   df_na2_bb$CI.min <-CI_df_na2_bb$`2.5 %` 
   df_na2_bb$CI.max <-CI_df_na2_bb$`97.5 %`
   df_na2_bb <- round(df_na2_bb, 3)
   df_na2_bb$Variable <- rownames(df_na2_bb)
   write.csv(df_na2_bb, "Data/Saved/Richness_nativeBark_modavg.csv")
   
 
   saveRDS(dredge_ric_na_bb, "Data/Saved/Richness_nativeBark_dredge.rds")
   
   
   
# ---- Plot ----  
   
   
   Richness_native_bark<-read.csv(file="Data/Saved/Richness_nativeBark_modavg.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, strip.white = T)
   Richness_native_bark$Variable <- gsub("cond\\((.*)\\)", "\\1",Richness_native_bark$Variable)
   
   
   Richness_native_bark <- Richness_native_bark %>% mutate(Type = case_when(
     Variable %in% c("Annual_Precipitation", "Minimum_Temperature") ~ "Climatic",
     Variable %in% c("Deciduous", "Evergreen", "Shannon_Div_Landscape") ~ "Landscape", 
     Variable %in% c("POPDEN", "Income_PerCapita", "PortDist") ~ "Socioeconomic", 
     TRUE ~ "Other"  ))
   Richness_native_bark <- Richness_native_bark[Richness_native_bark$Type != 'Other',]
   
   
   Richness_native_bark <- Richness_native_bark %>% filter(Variable != "(Intercept)")
   
   Richness_native_bark <- Richness_native_bark %>% mutate(
     p_value = Pr...z.. ,
     Significant = ifelse(p_value < 0.05, TRUE, FALSE)
   )
   
   Richness_native_bark$Variable <- factor(Richness_native_bark$Variable, 
                                           levels = c("Annual_Precipitation", "Minimum_Temperature","Deciduous", "Evergreen",  
                                                        "Shannon_Div_Landscape", "POPDEN", "Income_PerCapita", "PortDist"))
   
   
   Richness_native_bark <- Richness_native_bark %>%
     mutate(Variable = dplyr::recode(Variable,
                              "Annual_Precipitation" = "Annual Precipitation",
                              "Minimum_Temperature" = "Minimum Temperature",
                              "Deciduous" = "Deciduous Forest",
                              "Evergreen" = "Evergreen Forest",
                              "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                              "POPDEN" = "Human Population Density",
                              "Income_PerCapita" = "Per Capita Income",
                              "PortDist" = "Distance from Port"
     ))
   plot_na_bark_ric <- ggplot(Richness_native_bark, aes(x = Estimate, y = fct_rev(Variable) , color = Type)) +
     geom_point(aes(shape=Significant),size = 5,  stroke=1.5, fill="white") +
     geom_errorbar(aes(xmin = CI.min, xmax = CI.max), width = 0.1, linewidth = 1, orientation = "y") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
     # facet_wrap( ~Type , scales = "free_y", ncol = 1, strip.position="top") +
     theme_minimal(base_size = 14) +
     labs(
       x = "Coefficients",
       y = "Predictors",
       color = "Predictor Type",
       shape = NULL, 
       title = "Native Bark Beetles"
     ) +
     scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) + 
     scale_color_manual(values = c("Climatic" = "purple4", "Landscape" = "palegreen3", "Socioeconomic" = "orange2")) +
     annotate("text",x = Inf,y = Inf,label = "F",hjust = 1.3,vjust = 1.3,fontface = "bold",size = 7) +
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major.y = element_blank(), 
       panel.grid.minor.y = element_blank(),
       scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))),
       plot.title = element_text(hjust = 0.5, face="bold", size=20),
       strip.text = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, linewidth = 2), 
       axis.text.x = element_text(face = "bold", size = 18),
       axis.text.y = element_text(face = "bold", size = 18),
       axis.title.x = element_text(size = 22, face = "bold", vjust = -1),
       axis.title.y = element_text(size = 22, face = "bold", vjust = 2),
       legend.position = "none",
     )+ guides(shape="none")
   
   
   plot_na_bark_ric
   
   ggsave("Plots/Richness_native_bark.jpg", plot = plot_na_bark_ric, width = 10, height = 12, dpi = 600)
   
######################################################
######################################################
######################################################
   
# ---- Combined abundance plot ----
   
   library(patchwork)
   
   
  
   
   layout <- "
  AA
  BC
  DE
  FG
"
   
   
   # Combine using wrap_plots with axis collection
   abundance_combined_all <- wrap_plots(
     plot1,
     plot_na_abd, 
     plot_nn_abd, 
     plot_na_ambrosia_abd, 
     plot_nn_ambrosia_abd, 
     plot_na_bark_abd, 
     plot_nn_bark_abd, 
     design = layout
   ) + 
     plot_layout(axis_titles = "collect" , axes = "collect")+ 
     plot_annotation(
       title = "Abundance",
       
       theme = theme(plot.title = element_text(hjust = 0.6, vjust= 1, size = 24, face = "bold"))
     ) & 
     theme(
           plot.margin = margin(t = 10, r = 30, b = 10, l = 30), 
           axis.title.x = element_text(vjust = -3), 
           axis.title.y = element_text(vjust = 10),
           
     )
   
   
   
   abundance_combined_all
   
   
   ggsave("Plots/abundance_combined_all.jpg", plot = abundance_combined_all, width = 20, height = 24, dpi = 600)
   
# ---- Combined richness plot ----
   
   
   
   
   
   layout <- "
  AA
  BC
  DE
  FG
"
   
   
   
   # Combine using wrap_plots with axis collection
   richness_combined_all <- wrap_plots(
     plot2, 
     plot_na_ric, 
     plot_nn_ric, 
     plot_na_ambrosia_ric, 
     plot_nn_ambrosia_ric, 
     plot_na_bark_ric, 
     plot_nn_bark_ric, 
     design = layout
   ) + 
     plot_layout(axis_titles = "collect")+ 
     plot_annotation(
       title = "Species Richness",
    
       theme = theme(plot.title = element_text(hjust = 0.6, vjust= 1, size = 24, face = "bold"))
     ) & 
     theme(
           plot.margin = margin(t = 10, r = 30, b = 10, l = 30), 
           axis.title.x = element_text(vjust = -3), 
           axis.title.y = element_text(vjust = 10),
           
     )
   
   
   
   richness_combined_all
   
   
   ggsave("Plots/richness_combined_all.jpg", plot = richness_combined_all, width = 18, height = 22, dpi = 600)
   
   
###########################################################

# ---- Plotting the trap locations over US map ----  
   
  
   edrr_sf <- st_as_sf(edrr, coords = c("Longitude", "Latitude"), crs = "EPSG:4326")  
   us_map<- st_read("Data/Raw/cb_2018_us_state_5m.shp")
   edrr_sf<-st_transform(edrr_sf, st_crs(us_map))
   
   
   # No grid
   trapsmap <- ggplot() +
     geom_sf(data = us_map, fill = "white", color = "black") +
     
     geom_sf(data = edrr_sf, shape=21, fill = "skyblue", color= "black", size = 3,  alpha = 0.04) +
     
     coord_sf(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
     
     theme_minimal() +
     # theme(panel.grid = element_blank()) +
     labs(
       title = "",
       x = "", y = ""
     )
   trapsmap
   
   # With grid
   trapsmap1 <- ggplot() +
     geom_sf(data = us_map, fill = "white", color = "black") +
     
     geom_sf(data = edrr_sf, shape=21, fill = "darkorange", color= "black", size = 3,  alpha = 0.04) +
     
     coord_sf(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
     
     theme_minimal() +

     labs(
       title = "",
       x = "", y = ""
     )
   
   trapsmap1
    
   ggsave("Plots/trapsmap.jpg", plot = trapsmap, dpi = 600)
   
   ggsave("Plots/trapsmap1.jpg", plot = trapsmap1, dpi = 600)
   
   

  
  
   
# ---- Abundance table for model averages ----
  
   
   library(dplyr)
   
   process_abundance <- function(file_path, group_name) {
     df <- read.csv(file_path)
     
     df$Variable <- gsub("cond\\((.*)\\)", "\\1", df$Variable)
     df <- df[, -c(1, 4)]
     df <- df %>% relocate(7)
     df <- data.frame(Group = c(group_name, rep("", nrow(df) - 1)),df)
     df <- df %>% arrange(Variable)
     return(df)
   }
   
   
   Abundance_modavg <- process_abundance("Data/Saved/Abundance_modavg.csv","All Scolytines")
   
   Abundance_nonnative_modavg <- process_abundance("Data/Saved/Abundance_nonnative_modavg.csv","Non-native Scolytines")
   
   Abundance_nonnativeAmbrosia_modavg <- process_abundance("Data/Saved/Abundance_nonnativeAmbrosia_modavg.csv","Non-native Ambrosia Beetles")
   
   Abundance_nonnativeBark_modavg <- process_abundance("Data/Saved/Abundance_nonnativeBark_modavg.csv","Non-native Bark Beetles")
   
   Abundance_native_modavg <- process_abundance("Data/Saved/Abundance_native_modavg.csv","Native Scolytines")
   
   Abundance_nativeAmbrosia_modavg <- process_abundance("Data/Saved/Abundance_nativeAmbrosia_modavg.csv","Native Ambrosia Beetles")
   
   Abundance_nativeBark_modavg <- process_abundance("Data/Saved/Abundance_nativeBark_modavg.csv","Native Bark Beetles")
   
  
   
   combined_abundance_table <- rbind(Abundance_modavg,
                     Abundance_native_modavg, Abundance_nonnative_modavg, 
                     Abundance_nativeAmbrosia_modavg,    Abundance_nonnativeAmbrosia_modavg, 
                     Abundance_nativeBark_modavg, Abundance_nonnativeBark_modavg)
   
   
   combined_abundance_table <- combined_abundance_table %>%
     mutate(Variable = dplyr::recode(Variable,
                                     "(Int)" = "Intercept",
                                     "Annual_Precipitation" = "Annual Precipitation",
                                     "Minimum_Temperature" = "Minimum Temperature",
                                     "Deciduous" = "Deciduous Forest",
                                     "Evergreen" = "Evergreen Forest",
                                     "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                                     "POPDEN" = "Human Population Density",
                                     "Income_PerCapita" = "Per Capita Income",
                                     "PortDist" = "Proximity to Port"
     ))
   
   
   combined_abundance_table <- combined_abundance_table %>% 
     rename(
       'Predictor Variable' = Variable,
       'Model-Averaged Estimate' = Estimate,
       SE = Std..Error, 
       'z-value' = z.value, 
       'p-value' = Pr...z.., 
       LCL = CI.min,
       UCL = CI.max
     )   
   
   
   write.csv(combined_abundance_table, "Data/Saved/combined_abundance_table.csv", row.names = FALSE)
# ---- Richness table for model averages ----
   
   
   library(dplyr)
   
   process_richness <- function(file_path, group_name) {
     df <- read.csv(file_path)
     
     df$Variable <- gsub("cond\\((.*)\\)", "\\1", df$Variable)
     df <- df[, -c(1, 4)]
     df <- df %>% relocate(7)
     df <- data.frame(Group = c(group_name, rep("", nrow(df) - 1)),df)
     df <- df %>% arrange(Variable)
     return(df)
   }
   
   
  Richness_modavg <- process_richness("Data/Saved/Richness_modavg.csv","All Scolytines")
   
  Richness_nonnative_modavg <- process_richness("Data/Saved/Richness_nonnative_modavg.csv","Non-native Scolytines")
   
  Richness_nonnativeBark_modavg <- process_richness("Data/Saved/Richness_nonnativeBark_modavg.csv","Non-native Bark Beetles")
   
  Richness_native_modavg <- process_richness("Data/Saved/Richness_native_modavg.csv","Native Scolytines")
   
  Richness_nativeAmbrosia_modavg <- process_richness("Data/Saved/Richness_nativeAmbrosia_modavg.csv","Native Ambrosia Beetles")
   
  Richness_nativeBark_modavg <- process_richness("Data/Saved/Richness_nativeBark_modavg.csv","Native Bark Beetles")
   
  Richness_nonnativeAmbrosia_modavg <- process_richness("Data/Saved/Richness_nonnativeAmbrosia_modavg.csv","Non-native Ambrosia Beetles")
  
 
  
  #  Richness_nonnativeAmbrosia_modavg <- read.csv("Data/Saved/Richness_nonnativeAmbrosia_modavg.csv")
  # 
  # 
  # Richness_nonnativeAmbrosia_modavg$Variable <- gsub("cond\\((.*)\\)", "\\1", Richness_nonnativeAmbrosia_modavg$Variable)
  # Richness_nonnativeAmbrosia_modavg <- Richness_nonnativeAmbrosia_modavg[, -1]
  # Richness_nonnativeAmbrosia_modavg <- Richness_nonnativeAmbrosia_modavg %>% relocate(7)
  # Richness_nonnativeAmbrosia_modavg <- data.frame(Group = c("Non-native Ambrosia Beetles", rep("", nrow(Richness_nonnativeAmbrosia_modavg) - 1)),Richness_nonnativeAmbrosia_modavg)
  # Richness_nonnativeAmbrosia_modavg <- Richness_nonnativeAmbrosia_modavg %>% arrange(Variable)
   
   combined_richness_table <- rbind(Richness_modavg,
                                    Richness_native_modavg, Richness_nonnative_modavg, 
                                    Richness_nativeAmbrosia_modavg,   Richness_nonnativeAmbrosia_modavg, 
                                    Richness_nativeBark_modavg, Richness_nonnativeBark_modavg)
   
   
   combined_richness_table <- combined_richness_table %>%
     mutate(Variable = dplyr::recode(Variable,
                                     "(Int)" = "Intercept",
                                     "(Intercept)" = "Intercept",
                                     "Annual_Precipitation" = "Annual Precipitation",
                                     "Minimum_Temperature" = "Minimum Temperature",
                                     "Deciduous" = "Deciduous Forest",
                                     "Evergreen" = "Evergreen Forest",
                                     "Shannon_Div_Landscape" = "Shannon's Diversity (Landscape)",
                                     "POPDEN" = "Human Population Density",
                                     "Income_PerCapita" = "Per Capita Income",
                                     "PortDist"  = "Proximity to Port"
     ))

   
   combined_richness_table <- combined_richness_table %>% 
     rename(
       'Predictor Variable' = Variable,
       'Model-Averaged Estimate' = Estimate,
       SE = Std..Error, 
       'z-value' = z.value, 
       'p-value' = Pr...z.., 
       LCL = CI.min,
       UCL = CI.max
       )   
   
   
   write.csv(combined_richness_table, "Data/Saved/combined_richness_table.csv", row.names = FALSE)
   
   
# ---- Combined table for models with separately fit correlated predictors ----
   model_comparison <- rbind(modelcomp1, modelcomp2, modelcomp3, modelcomp4)
   
   write.csv(model_comparison, "Data/Saved/model_comparison_table.csv", row.names = FALSE)
   