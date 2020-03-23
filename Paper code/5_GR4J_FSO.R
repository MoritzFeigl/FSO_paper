#
# Function Space Optimization
# Case study GR4J Optimization & testing
#
# Moritz Feigl
#


# 0. Load everything needed 
setwd("D:/Dropbox/Diss/CF_Grammar/Code/GR4J_case_study")
source("Functions/FSO_functions.R")
FSO_setup()

# 1. Optimization, testing & diagnostic plots
grid <- expand.grid("Test_number" = c(1.1, 1.2, 2.4, 2.5, 2.6),
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
FSO_plot(Test_number = 3.1)
