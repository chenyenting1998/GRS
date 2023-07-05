#######################
# Setup - Cruise color
#######################

# Description:
# Assigning color-blind-friendly colors to each cruise for data visualization

# Author: Yen-Ting Chen
# Date of creation: Unknown
# Date of last modification: 2023/07/05

cruise_color <- c("OR1-1219" = "#DC3220", # red
                  "OR1-1242" = "#0C7BDC") # blue

save(cruise_color, file = "data/cruise_color.RData")