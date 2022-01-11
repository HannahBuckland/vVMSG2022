# Function to clean probe data

# Hannah Buckland
# Updated for vVMSG2022

data_clean <- function(x){
  
  mean_TOTAL <- mean(x$TOTAL)
  sd_TOTAL <- sd(x$TOTAL)
  
  probe_dat_clean <- x %>%
    filter(TOTAL <= mean_TOTAL + 2*sd_TOTAL & TOTAL >= mean(TOTAL) - 2*sd(TOTAL)) %>%
    filter(Al2O3 < 25) 
  
  return(probe_dat_clean)
}