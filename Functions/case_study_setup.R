# Setup for d-GR4J case study
# Moritz Feigl, 2019
#

# 1. define "True" functions -------------------------------------------------------------
transfer_functions <- list(GR4Jx1 = "0.5 + evi*1.5 + exp(bdim)*0.9", # production store
                           GR4Jx2 = "0",          # interaction direct runoff and routing store
                           GR4Jx3 = "-1.3  - log(slope)",  # routing store
                           GR4Jx4 = "-log(hand)*0.2 - slope*1.5 - 1.5") # length of UH

# 2. define gmean ------------------------------------------------------------------------
gmean_parameter <- 1

# 3. define test and training  basins ----------------------------------------------------
training_basins <- c(1, 3, 4, 8, 27, 30, 35, 38, 46, 47, 50, 51, 57, 60, 69, 72, 75, 84, 
                     90, 93, 95, 96, 98, 100, 101, 103, 111) # 27 training basins
test_basins <- c(2, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                 21, 22, 23, 24, 25, 26, 28, 29, 31, 32, 33, 34, 36, 37, 39, 40,
                 41,  42, 43, 44, 45, 48, 49, 52, 53, 54, 55, 56, 58, 59, 61, 62,
                 63, 64, 65, 66, 67, 68, 70, 71, 73, 74, 76, 77, 78, 79, 80, 81,
                 82, 83, 85, 86, 87, 88, 89, 91, 92, 94, 97, 99, 102, 104, 105, 106,
                 107, 108, 109, 110, 112)



