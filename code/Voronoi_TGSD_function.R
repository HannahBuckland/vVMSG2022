## Function to calculate the TGSD of a deposit
# Function adapted from the TOTGS Matlab code by Biass and Bonadonna, 2014 - https://github.com/e5k/TOTGS/tree/v2.1

# Author - Dr Hannah Buckland
# 06/09/2021 (Adapted for vVMSG2022 GitHub repo on 26/01/2022)

library(rgdal)
library(ggvoronoi)
library(rgeos)
library(tidyverse)
library(viridis)
library(patchwork)
library(rnaturalearth)
library(sf)
library(here)
library(data.table)
library(rworldmap)

#### Enter file names for the csv file with GSD at each locality (datapoints),
#### the csv file containing the lat long of the zero line (zeros) and 
#### the lat long of your volcano (volc_ll) ####
TGSD_input <- list(datapoints=here("data","MSH_TGSD_dataset.csv"),
                   zeros=here("data","MSH_TGSD_zeroline.csv"),
                   volc_ll = data.frame(vol_lat=46.2,vol_long=-122.18),
                   utm_zone = "")

# Here please input Y or N depending on whether you want the outputs to save (be careful of overwriting)
save_condition <- "Y"

################### Shouldn't need to change below this line ##################

#### Function to read in the GSDs and zero points ####
readdat <- function(TGSD_input){
  datapoints <- TGSD_input$datapoints
  zeros <- TGSD_input$zeros
  
  llpoints <- read.csv(datapoints) # read in the GSD data
  llzeros <- read.csv(zeros) # read in the lat long of the zero line
  
  dataread <- list(llpoints,llzeros) # combine into a list
  names(dataread) <- c("llpoints","llzeros")
  return(dataread)
}

voronoi_input <- readdat(TGSD_input) # apply the readdat function to the input

#### Function to set the UTM zone for the spatial calculations ####
utmzone <- function(GSDpoints){
  
  # Each UTM zone is 6 degrees wide, so you can get the zone number by looking at the longitude
  utmzone <- median(floor((GSDpoints$long + 180) / 6) + 1) # here we use the median as data can span multiple zones
  # now test whether the zone is positive or negative (north or south)
  if(min(GSDpoints$lat) > 0) {
    proj_utm <- utmzone + 32600
  } else {
    proj_utm <- utmzone + 32700
    }
  proj_CRS <- st_crs(proj_utm)$proj4string
  return(proj_CRS)
}

utm_convert <- utmzone(voronoi_input$llpoints) # Apply utm function to GSD points file

#### Voronoi TGSD calculation function ####

TGSD_voronoi <- function(input){
  
  llpoints <- input$llpoints # isolate the GSD data
  
  llpoints <- llpoints %>%
    select(lat,long,gm2,starts_with("phi")) # selects only the columns we need for analysis
  
  llzeros <- input$llzeros 
  
  fullpoints <- bind_rows(llpoints,llzeros) # combine the data and zeros for Voronoi tesslation
  fullpoints[is.na(fullpoints)] <- 0 # replaces NA values for zero line with zero
  
  # Need to convert from geographic coordinates to projected for area calculations
  ll <- cbind(fullpoints$long, fullpoints$lat)
  ll <- sp::SpatialPoints(ll, proj4string = CRS("+proj=longlat"))
  
  utm <- sp::spTransform(ll, CRS(utm_convert)) # transform from lat long to UTM (projected)
  
  vdf <- as.data.frame(cbind(coordinates(utm),fullpoints$gm2)) # get utm coordinates and mass into data frame
  colnames(vdf) <- c("x","y","z") 
  
  vordat_spdf <- ggvoronoi::voronoi_polygon(vdf, x="x",y="y") # voronoi tesselation
  vor_area <- rgeos::gArea(vordat_spdf,byid = TRUE) #calculate area of each voronoi cell
  vordat_df <- as.data.frame(cbind(vordat_spdf$x,vordat_spdf$y,vordat_spdf$z,vor_area)) #assign area to lat long pair
  colnames(vordat_df) <- c("x","y","z","area")
  
  TGSDcalc <- subset(vordat_df, z > 0) #filter to exclude zero points
  
  calcpoints <- cbind(llpoints,TGSDcalc) #combine GSD dataframe with area dataframe
  calcpoints$totarea <- sum(calcpoints$area)
  calcpoints$prop <- calcpoints$area / calcpoints$totarea #weighting by proportion of voronoi cell area
  
  gsdpoints <- calcpoints[,grepl("phi",names(calcpoints))] #isolate the GSDs
  gsdmass <- gsdpoints/100 * calcpoints$gm2 #get percentage of mass
  gsdweighted <- gsdmass * calcpoints$prop
  
  TGSD <- colMeans(gsdweighted)
  totmass <- sum(TGSD)
  TGSD_perc <- TGSD/totmass*100
  
  phi_cols <- llpoints %>%
    select(starts_with("phi")) # select only the grain size columns
  
  phi_bins <- colnames(phi_cols) # use colnames to get bin spacing
  phi_bins <- data.table::tstrsplit(phi_bins,"phi")[[2]] # remove prefix
  phi_bins <- as.numeric(gsub("\\_", "-",phi_bins)) # convert to numeric and replace underscore with minus sign
  
  resPDF <- data.frame(cbind(phi_bins,TGSD_perc,cumsum(TGSD_perc)))
  colnames(resPDF) <- c("phi","PDF","CDF")
  rownames(resPDF) <- NULL
  
  # interpolate median and quartiles
  medphi <- approx(resPDF$CDF, resPDF$phi_bins, xout=50, method = "linear", ties = mean)
  phi84 <- approx(resPDF$CDF, resPDF$phi_bins, xout=84, method = "linear", ties = mean)
  phi16 <- approx(resPDF$CDF, resPDF$phi_bins, xout=16, method = "linear", ties = mean)
  phi95 <- approx(resPDF$CDF, resPDF$phi_bins, xout=95, method = "linear", ties = mean)
  phi5 <- approx(resPDF$CDF, resPDF$phi_bins, xout=5, method = "linear", ties = mean)
  
  # calculate the Folk and Ward standard deviation 
  sigmaphi <- ((phi84$y-phi16$y)/4) + ((phi95$y+phi5$y)/6.6)
  
  outVoronoi <- list(vordat_spdf,
                     resPDF,
                     fullpoints,
                     round(medphi$y,4),
                     round(sigmaphi,4))
  names(outVoronoi) <- c("SPDF","PDF","input","Md50","StDev")
  return(outVoronoi)
  
}

out_voronoi<-TGSD_voronoi(voronoi_input) # apply the Voronoi function to the input

#### Function to convert Voronoi output into spatial data, useful for plotting ####

spatial_proc <- function(out_voronoi){
  
  map <- sf::st_as_sf(out_voronoi$SPDF) # convert to sf object
  map_utm <- sf::st_set_crs(map,CRS(utm_convert))
  points <- sf::st_as_sf(out_voronoi$input, 
                     coords = c("long", "lat"), crs = 4326, 
                     agr = "constant")
  dfv <- out_voronoi$PDF
  
  out <- list(map_utm,points,dfv)
  names(out) <- c("map_sf","points_sf","TGSD")
  return(out)
}

sf_out <-spatial_proc(out_voronoi) # apply spatial function to Voronoi output


#### Plotting functions to produce histogram of TGSD and map of Voronoi cells ####

# Function to find out what countries are required for the basemap
coords2country <- function(zeropoints){
  countriesSP <- rworldmap::getMap(resolution='low') # get data from world map
  llonly <- zeropoints %>%
    select(long,lat) # order the longitude and latitude of the zeros
  
  #setting CRS directly to that from rworldmap
  pointsSP <- sp::SpatialPoints(llonly, proj4string=CRS(proj4string(countriesSP)))  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices <- sp::over(pointsSP, countriesSP)
  
  # return the ADMIN names of each country
  countries <- indices %>%
    distinct(ADMIN)
  
  country_list <- as.matrix(countries)
  return(country_list)
}

countries_out <- coords2country(voronoi_input$llzeros)

# Plotting function
TGSD_plots <- function(input){
  statesdata <- rnaturalearth::ne_states(country = countries_out,
                                         returnclass = 'sf')
  
  # get location of volcano/vent for plotting
  volcano <- TGSD_input$volc_ll
  
  # convert to spatial object
  volc_sf <- sf::st_as_sf(volcano,
                      coords = c("vol_long","vol_lat"),crs = 4326,
                      agr = "constant")
  
  # set theme for consistent plot style
  theme <-
    theme_set(theme(legend.position = "none",
                    legend.title = element_text(size =9),
                    legend.text = element_text(size =8),
                    legend.key.size = unit(0.3,"cm"),
                    legend.background = element_rect(fill = "grey99",colour="grey20"),
                    axis.text=element_text(colour="black"),
                    plot.background = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_rect(colour="black",fill=NA),
                    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")))
  
  # set axis breaks for difference grain size scales (linear and log/phi)
  mm_plotting <- data.frame(breaks= -log2(c(1e-3,10e-3,100e-3,1000e-3,10000e-3)),
                            labels = c(0.001,0.01,0.1,1,10))
  phi_plotting <- data.frame(breaks = seq(min(sf_out$TGSD$phi),max(sf_out$TGSD$phi),
                                          by=2))
  
  y_ax_fact <- 100/ceiling(max(sf_out$TGSD$PDF)) # factor for getting CDF and PDF on same plot
  plot_limits <- data.frame(xmin=ceiling(min(voronoi_input$llzeros$long))-2,
                            xmax=ceiling(max(voronoi_input$llzeros$long))+2,
                            ymin=ceiling(min(voronoi_input$llzeros$lat))-2,
                            ymax=ceiling(max(voronoi_input$llzeros$lat))+2)
  
  
  TGSDplt <- ggplot(data=sf_out$TGSD) +
    geom_col(aes(x=phi,y=PDF),fill= "#999999") +
    geom_line(aes(x=phi,y=CDF/y_ax_fact)) +
    scale_y_continuous(expand = c(0,0),
                       name = "Frequency %",
                       sec.axis = sec_axis(~.*y_ax_fact,
                                           name="Cumulative %")) +
    scale_x_reverse(breaks = phi_plotting$breaks, 
                    labels = phi_plotting$breaks,
                    name = "Grain size (phi)",
                    sec.axis = sec_axis(trans = ~.*1, 
                                        name = "Grain Size (mm)", 
                                        breaks = mm_plotting$breaks,
                                        labels = mm_plotting$labels)) +
    theme(aspect.ratio = 1)
  
  vormap <- ggplot(data=statesdata) +
    geom_sf(fill="grey60") +
    geom_sf(data=sf_out$map_sf %>% filter(z==0),
            fill = NA,colour = "grey80") + 
    geom_sf(data=sf_out$map_sf %>% filter(z>0),
            aes(fill=log(z)),alpha=0.7,colour = "black") +
    geom_sf(data=sf_out$points_sf,size=0.5) +
    geom_sf(data=volc_sf, fill="red",colour="black",pch=24,size=2) +
    scale_fill_viridis(option = "viridis",direction=-1,
                       guide=guide_colorbar(direction="horizontal")) +
    scale_x_continuous(breaks = seq(-124,-112,by=4)) +
    scale_y_continuous(breaks = seq(44,50,by=2)) +
    coord_sf(xlim=c(plot_limits$xmin,plot_limits$xmax),
             ylim=c(plot_limits$ymin,plot_limits$ymax)) +
    labs(fill = expression("Log(g/m"^2*")")) +
    theme(legend.position = "bottom")
  
  out_plots <- TGSDplt + vormap
  return(out_plots)
  
}

results_plotted <- TGSD_plots()

out_voronoi$PDF
results_plotted


#### User interaction for saving ####

if(save_condition == "Y") {
  write.csv(out_voronoi$PDF,file(here("data","TGSD_out.csv")),row.names = FALSE)
  ggsave(file=(here("plots","TGSD_plots.pdf")),results_plotted)
  
  print("Outputs have been saved in `plots` and `data`")
} else{
  print("Outputs not saved")
}



