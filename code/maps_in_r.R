# Plotting a map in R
# Hannah Buckland

# modified from exisiting scripts for vVMSG

# Read in packages for dealing with spatial data
library(rnaturalearth)
library(sf)
library(tidyverse)

# Read in data base of Mazama localities

MLDB <- read.csv("../data/MLDB_latlong.csv")

# Now convert this to a spatial data set for plotting
MLDB_sf <- st_as_sf(MLDB, 
                   coords = c("Longitude", "Latitude"), crs = 4326, 
                   agr = "constant") # convert to spatial dataset

# Read in data for plotting map outlines of countries
statesdata <- ne_states(country =c("United States of America","Canada","Mexico"),
                        returnclass = 'sf')

volcano <- data.frame(name="Mazama",
                      latitude = 42.9446,
                      longitude = -122.1090) # data frame with location of mount Mazama

volc_sf <- st_as_sf(volcano,
                    coords = c("longitude","latitude"),crs = 4326, 
                    agr = "constant") # convert to spatial dataset

ggplot(data=statesdata) +
  geom_sf(fill = "#F4F9E9",colour = "#2F323A",size=0.2) +
  geom_sf(data=MLDB_sf,
          aes(colour=Locality_MLDB,
              size=n)) +
  geom_sf_text(data = loc_sf,
               aes(label=Locality_MLDB),
               nudge_x = 0.25,
               nudge_y = 0.25,
               size=3) +
  geom_sf(data=volc_sf,
          pch=17,
          size=4,
          colour = "black") +
  scale_size_binned(breaks= c(5,10,25,50,100)) +
  labs(color="MLDB Number", size = "Number of Analyses") +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Visualising the spatial distribution of Mazama glass analyses") + 
  coord_sf(xlim=c(-125,-115),ylim=c(40.5,50)) +
  theme(panel.background = element_rect(fill="#B4D6D3"),
        legend.position = "right")
