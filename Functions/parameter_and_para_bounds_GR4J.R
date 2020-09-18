# get parameter table and bounds for d-GR4J case study
# Moritz Feigl, 2019
#


# 2 KM GR4J Parameters
para <- read.table("GR4J_distributed/input/para_Mur_GR4J_true.txt",
                   skip = 2, sep = " ", header = FALSE)
colnames(para) <- read.table("GR4J_distributed/input/para_Mur_GR4J_true.txt",
                             skip = 1, sep = " ", header = FALSE, nrows = 1, stringsAsFactors = FALSE)
# 1 KM GR4J Parameters
para_1km <- read.table("GR4J_distributed_1km/input/para_Mur_GR4J_true.txt",
                       skip = 2, sep = " ", header = FALSE)
colnames(para_1km) <- read.table("GR4J_distributed_1km/input/para_Mur_GR4J_true.txt",
                                 skip = 1, sep = " ", header = FALSE, nrows = 1, stringsAsFactors = FALSE)
# Parameter bounds
bounds <- read.table("GR4J_distributed/parameter.txt",
                     skip = 1, sep = "\t", header = FALSE)
parameter_bounds <- list(GR4Jx1 = as.numeric(bounds[7, c(3, 4)]), # production store
                         GR4Jx2 = as.numeric(bounds[8, c(3, 4)]),          # transboundary flow
                         GR4Jx3 = as.numeric(bounds[9, c(3, 4)]),  # routing store
                         GR4Jx4 = as.numeric(bounds[10, c(3, 4)])) #length of the UH