# Function to clean probe data

# Hannah Buckland
# Updated for vVMSG2022

data_clean <- function(x){
  
  mean_TOTAL <- mean(x$TOTAL) # calculates the mean of the analytical totals
  sd_TOTAL <- sd(x$TOTAL) # calculates the standard deviation of the analytical totlals
  
  probe_dat_clean <- x %>%
    filter(TOTAL <= mean_TOTAL + 2*sd_TOTAL & TOTAL >= mean(TOTAL) - 2*sd(TOTAL)) %>% # excludes analyses with anomalous totals
    filter(Al2O3 < 25) # removes analyses that you likely hit a crystal
  
  return(probe_dat_clean) # returns the cleaned data frame
}