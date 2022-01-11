# Plotting a map in R
# Hannah Buckland

# modified from exisiting scripts for vVMSG

# Read in packages for dealing with spatial data
library(rnaturalearth)
library(sf)
library(tidyverse)


#### Data processing ####
# Read in data base of Mazama localities

MLDB <- read.csv("data/MLDB_latlong.csv")

# Now convert this to a spatial data set for plotting
MLDB_sf <- st_as_sf(MLDB, 
                   coords = c("Longitude", "Latitude"), crs = 4326, 
                   agr = "constant") # convert to spatial dataset

# Read in data for plotting map outlines of states 
statesdata <- ne_states(country =c("United States of America","Canada"),
                        returnclass = 'sf')

# Set up a dataframe with the location of Crater Lake/Mount Mazama
volcano <- data.frame(name="Mazama",
                      latitude = 42.9446,
                      longitude = -122.1090) # data frame with location of mount Mazama

# convert to this lat long dataframe to a spatial dataset
volc_sf <- st_as_sf(volcano,
                    coords = c("longitude","latitude"),crs = 4326, 
                    agr = "constant") 

# read in a shapefiles of the Mazama isopachs from Buckland et al. 2020
onecm_iso <- st_read("data/1cm_isopach.shp") # 1cm hand drawn
spline_isos <- st_read("data/distal200.shp") # Buckland et al. 2020 spline isopachs


#### Plotting code #####
loc_map <- ggplot(data=statesdata) +
  geom_sf(fill = "#F4F9E9",
          colour = "#2F323A",
          size=0.2) +
  geom_sf(data=onecm_iso,
          fill=NA, 
          colour = "grey58",
          lty = 2) +
  geom_sf(data=spline_isos,
          fill=NA,
          colour = "grey58") +
  geom_sf(data=MLDB_sf,
          colour="black",
          size = 1) +
  geom_sf(data=volc_sf,
          pch=24,
          size=4,
          colour = "black",
          fill = "red") +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Localities in the Mazama Locality Database") + 
  coord_sf(xlim=c(-125,-105),ylim=c(38,55)) +
  theme(panel.background = element_rect(fill="#B4D6D3"),
        legend.position = "right")


ggsave("plots/MLSB_localities.png",loc_map,width=5,height=7,dpi=300)