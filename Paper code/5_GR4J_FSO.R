# Function Space Optimization case study with d-GR4J
# Moritz Feigl, 2019
#


# 0. Load everything needed 
setwd("FSO_paper")
source("Functions/FSO_functions.R")
FSO_setup()

# 1. Optimization, testing & diagnostic plots
grid <- expand.grid("Test_number" = c(1.1, 1.2, 2.4, 2.5, 2.6, 4.4, 4.5, 4.6),
                   "Optimizer" = c("GA", "DDS", "PSO"),
                   "run" = c(1:5),
                   stringsAsFactors = FALSE)

mapply(FSO,
      Optimizer = grid$Optimizer,
      Test_number = grid$Test_number,
      run = grid$run,
      MoreArgs = list("iterations" = 5000,
                      "training_basins" = training_basins, 
                      'test_basins' = test_basins))


# 2. Plot results
source("Paper code/FSO_plots") 
