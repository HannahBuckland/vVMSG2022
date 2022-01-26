## Script to execute the Voronoi tesselation method of calculating the TGSD of a deposit
# Method adapted from the TOTGS Matlab code by Biass and Bonadonna, 2014 - https://github.com/e5k/TOTGS/tree/v2.1
# The individual functions used can be found in "Voronoi_TGSD_functions.R"

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
TGSD_input <- list(
  datapoints = here("data", "MSH_TGSD_dataset.csv"),
  zeros = here("data", "MSH_TGSD_zeroline.csv"),
  volc_ll = data.frame(vol_lat = 46.2, vol_long = -122.18),
  utm_zone = ""
)

# Here please input Y or N depending on whether you want the outputs to save (be careful of overwriting)
save_condition <- "Y"


################### ################### ################### ###################
################### Shouldn't need to change below this line# ##################
################### ################### ################### ###################

source(here("code","Voronoi_TGSD_functions.R"))

voronoi_input <-
  readdat(TGSD_input) # apply the readdat function to the input

utm_convert <-
  utmzone(voronoi_input$llpoints) # Apply utm function to GSD points file

out_voronoi <-
  TGSD_voronoi(voronoi_input) # apply the Voronoi function to the input

sf_out <-
  spatial_proc(out_voronoi) # apply spatial function to Voronoi output

countries_out <- 
  coords2country(voronoi_input$llzeros)

results_plotted <- TGSD_plots(TGSD_input)

out_voronoi$PDF
results_plotted

#### Conditional saving ####

if (save_condition == "Y") {
  write.csv(out_voronoi$PDF, file(here("data", "TGSD_out.csv")), row.names = FALSE)
  ggsave(file = (here("plots", "TGSD_plots.pdf")), results_plotted)
  
  print("Outputs have been saved in `plots` and `data`")
} else{
  print("Outputs not saved")
}

