#
# Function Space Optimization
# Case study GR4J Setup
#
# Moritz Feigl
#


setwd("FSO_paper")

if (!file.exists("True parameters")){
  dir.create("True parameters")
}
if (!file.exists("True parameters/Plots")){
  dir.create("True parameters/Plots")
}

# 0. Load Data and define help functions -----------------------------------------------
# get True functions
source("Functions/case_study_setup.R")
# Load helper functions
source("Functions/FSO_functions.R")
# load l0 layer
l0 <- load_sp_mur(scale = TRUE, na.approx = TRUE, 
                  only_training_basins = FALSE, full_dataset = FALSE)
# load 2km parameter file
para <- read.table("GR4J_distributed/input/para_Mur_GR4J_start.txt",
                   skip = 2, sep = " ", header = FALSE)
colnames(para) <- read.table("GR4J_distributed/input/para_Mur_GR4J_start.txt",
                             skip = 1, sep = " ", header = FALSE, nrows = 1, stringsAsFactors = FALSE)
para <- para[, names(para) != "NA"]
# Load 1km Parameter file
para_1km <- read.table("GR4J_distributed_1km/input/Para_Mur_optNSE_1km_24h.txt",
                       skip = 2, sep = " ", header = FALSE)
colnames(para_1km) <- read.table("GR4J_distributed_1km/input/Para_Mur_optNSE_1km_24h.txt",
                                 skip = 1, sep = " ", header = FALSE, nrows = 1, stringsAsFactors = FALSE)
para_1km <- para_1km[, names(para_1km) != "NA"]
# Load parameter bounds
bounds <- read.table("GR4J_distributed/parameter.txt",
                     skip = 1, sep = "\t", header = FALSE)
parameter_bounds <- list(GR4Jx1 = as.numeric(bounds[7, c(3, 4)]), # production store
                         GR4Jx2 = as.numeric(bounds[8, c(3, 4)]),          # transboundary flow
                         GR4Jx3 = as.numeric(bounds[9, c(3, 4)]),  # routing store
                         GR4Jx4 = as.numeric(bounds[10, c(3, 4)]))


# 1. Raster & Plots of true parameters -------------------------------------------------
for(parameter in c("x1", "x3", "x4")){
  true_250m <- raster_from_tf(tf = transfer_functions[[eval(paste0("GR4J", parameter))]], 
                              tf_bounds = parameter_bounds[[eval(paste0("GR4J", parameter))]])
  true_1km <- raster_from_tf(tf = transfer_functions[[eval(paste0("GR4J", parameter))]], 
                             tf_bounds = parameter_bounds[[eval(paste0("GR4J", parameter))]],
                             aggregate = TRUE, km1 = TRUE, gmean_parameter = gmean_parameter)
  true_2km <- raster_from_tf(tf = transfer_functions[[eval(paste0("GR4J", parameter))]], 
                             tf_bounds = parameter_bounds[[eval(paste0("GR4J", parameter))]],
                             aggregate = TRUE, km1 = FALSE, gmean_parameter = gmean_parameter)
  for(rasters in c("250m", "1km", "2km")){
    writeRaster(get(paste0("true_", rasters)), 
                paste0("True parameters/", parameter, "_", rasters, ".asc"), 
                format = "ascii", overwrite = TRUE)
    png(paste0("True parameters/Plots/", parameter, "_", rasters, ".png"),
        width = 1200, height = 800)
    plot(get(paste0("true_", rasters)))
    dev.off()
    
  }
}


# 2. Produce correct 1km and 2km parameter sets ----------------------------------------

# Change columns to be like the 2km version
para_1km <- para_1km[, names(para_1km) %in% names(para)]
missing_columns <- names(para)[!(names(para) %in% names(para_1km))]
add_columns <- as.data.frame(matrix(NA, nrow = nrow(para_1km), ncol = length(missing_columns)))
# set in values that are not in the 1km version
names(add_columns) <- missing_columns
add_columns$Div_TONZ <- para_1km$TONZ
add_columns[, c("Div_LowThr", "Div_UpThr")] <- 999
add_columns[, "CTGLAC_"] <- 0.015
add_columns[, "GLFAKTGLAC_"] <- 0.005
add_columns[, "ICECOVGLAC"] <- 1
add_columns[, "TABGLAC_"] <- 4
add_columns[, "GR4Jx1"] <- 1000
add_columns[, "GR4Jx2"] <- 0
add_columns[, "GR4Jx3"] <- 500
add_columns[, "GR4Jx4"] <- 4
add_columns[, "GR4JiniSt1"] <- 0.6
add_columns[, "GR4JiniSt2"] <- 0.8
para_1km <- cbind(para_1km, add_columns)
para_1km <- para_1km[, eval(names(para))]
cat("\n", file = "GR4J_distributed_1km/input/para_Mur_GR4J_true.txt")
write.table(para_1km[, 1:172], "GR4J_distributed_1km/input/para_Mur_GR4J_true.txt", append = TRUE,
            row.names = FALSE, quote = FALSE)

para_1km <- read.table("GR4J_distributed_1km/input/para_Mur_GR4J_true.txt",
                       skip = 2, sep = " ", header = FALSE)
colnames(para_1km) <- read.table("GR4J_distributed_1km/input/para_Mur_GR4J_true.txt",
                                 skip = 1, sep = " ", header = FALSE, nrows = 1, stringsAsFactors = FALSE)

if(ncol(para) != ncol(para_1km)) cat("Somethings wrong!")

# 
# # Produce correctly aggregated parameters for 2km resolution GR4J
km1_to_km2 <- feather::read_feather("Basin infos/km1_to_km2.feather")
km1_to_km2 <- km1_to_km2[order(km1_to_km2$X1km_NZ), ]
para_1km <- para_1km[order(para_1km$NZ_), ]
para_1to2 <- cbind(para_1km, km1_to_km2[, c("X1km_NZ", "X2km_NZ")])

#only aggregate parameter that make sense
para_2km <- aggregate(. ~ X2km_NZ, para_1to2[, -c(1:13, 173)], FUN = gmean, p = gmean_parameter)
para_2km <- para_2km[order(para_2km$X2km_NZ), ]

para_2km <- merge(para[, 1:13], para_2km, by.x = "NZ_", by.y = "X2km_NZ")
para_2km <- para_2km[order(para_2km$NZ_), ]
#para_2km <- cbind(para[, 1:13], para_2km)
para_2km$X2km_NZ <- NULL
para_2km[, 1:3] <- para_2km[, c(2, 3, 1)]
names(para_2km)[1:3] <- c("NB_", "IZ_", "NZ_")

# round to the same number of decimals as the para_1km columns
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
for(i in 14:172){
  decs <- decimalplaces(para_1km[, i])
  para_2km[, i] <- round(para_2km[, i], decs)
}
cat("\n", file = "GR4J_distributed/input/para_Mur_GR4J_true.txt")
write.table(para_2km, "GR4J_distributed/input/para_Mur_GR4J_true.txt", append = TRUE,
            row.names = FALSE, quote = FALSE)


# Load 2km parameters
para <- read.table("GR4J_distributed/input/para_Mur_GR4J_true.txt",
                   skip = 2, sep = " ", header = FALSE)
colnames(para) <- read.table("GR4J_distributed/input/para_Mur_GR4J_true.txt",
                             skip = 1, sep = " ", header = FALSE, nrows = 1, stringsAsFactors = FALSE)
para <- para[, names(para) != "NA"]

# 3. Create Synthetic parameters for 2km GR4J ------------------------------------------
new_gr4j_para <- create_GR4J_para(transfer_functions, l0, parameter_bounds, 
                                  gmean_parameter = gmean_parameter)

# merge new parameters with parameter file
para_new <- merge(para, new_gr4j_para, by = c("NB_", "NZ_"), suffixes =  c(".old", ""),
                  all.x = TRUE)
para_new[which(is.na(para_new), arr.ind = TRUE)] <- 0 # parameter of unused basins -> 0
para_new <- para_new[, -grep(".old", names(para_new))]
states <- para_new[, grep("iniSt", names(para_new))]
para_new <- para_new[, -grep("iniSt", names(para_new))]
para_new <- cbind(para_new, states, stringsAsFactors = FALSE)
para_new <- para_new[order(para_new$NZ_), ]

para_new[, 1:3] <- para_new[, c(1, 3, 2)]
names(para_new)[2:3] <- c("IZ_", "NZ_")

cat("\n", file = "GR4J_distributed/input/para_Mur_GR4J_true.txt")
write.table(para_new, "GR4J_distributed/input/para_Mur_GR4J_true.txt",
            append = TRUE, row.names = FALSE, quote = FALSE)

# 4. Q & states for 2km GR4J -----------------------------------------------------------

# save output as new Qobs file

# **************************  Change in defaults ******************************
# run GR4J_parallel for all basins with:
# para_Mur_GR4J_true.txt parameter file
# full time series: 2003 - 2012
# true Qobs
# all basins
# ************************************************************************
default <- readLines("GR4J_distributed/input/defaults.txt")
# set para file
default[14] <- "para_Mur_GR4J_true.txt"
default[15] <- "para_Mur_GR4J_fsOptim.txt"
# set dates
default[32] <- "2003 01 02 00 00" # start date
default[33] <- "2009 08 31 00 00" # backup
default[36] <- "2012 12 31 00 00" # end date
default[37] <- "2009 08 31 00 00" # backup
# set spin-up
default[40] <- "241"
default[41] <- "2433"
# set Datafile
default[7] <- "QObs_24h.txt" # used
default[8] <- "QObs_24h_synth_GR4J.txt" # backup
# set output type
default[26] <- "1"
# write new default file
writeLines(default, con = "GR4J_distributed/input/defaults.txt")

# Start Model run
setwd("GR4J_distributed/")
sys::exec_wait("start_GR4J_case_study_full_run.bat",
               std_out = "GR4J_output.txt")
setwd("..")

# Generate synthetic Q timeseries
synth_obs <- read.table("GR4J_distributed/output/GR4J.runoff", skip = 17, header = TRUE)
synth_obs <- synth_obs[, -c(grep("QOBS", names(synth_obs)),
                            grep("Qloc", names(synth_obs)))]
synth_obs[which(is.na(synth_obs), arr.ind = TRUE)] <- -999
cat("synthetic_Qobs_Mur_24h_mean_UTC   2003-2012\n",
    file = "GR4J_distributed/input/QObs_24h_synth_GR4J.txt")
for(i in 1:112){
  cat(paste0(gsub("QSIM", "QOBS", names(synth_obs)[5+i]),  "\n"),
      file = "GR4J_distributed/input/QObs_24h_synth_GR4J.txt", append = TRUE)
}
cat("##############################################\n",
    file = "GR4J_distributed/input/QObs_24h_synth_GR4J.txt", append = TRUE)
write.table(synth_obs, "GR4J_distributed/input/QObs_24h_synth_GR4J.txt",
            append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
# Generate true states time series
# read state results
state_names <- read.table("GR4J_distributed/output/COSERO.plus1", 
                          header=FALSE, fill=TRUE, skip = 1, nrows = 1, 
                          stringsAsFactors = FALSE)
state1_ind <- grep(state_names, pattern= "BW0GEB")
state2_ind <- grep(state_names, pattern= "BW3GEB")
state_classes <- rep("NULL", length(state_names))
state_classes[c(state1_ind, state2_ind)] <- "numeric"
state_classes[1:3] <- "integer"
states <- read.table("GR4J_distributed/output/COSERO.plus1", 
                     header=FALSE, fill=FALSE, skip = 243,
                     colClasses = state_classes)
state_classes[c(1:3, state1_ind, state2_ind)] <- "character"

names(states) <- read.table("GR4J_distributed/output/COSERO.plus1", 
                            header=FALSE, fill=FALSE, skip = 1,nrows = 1,
                            colClasses = state_classes)
# relevanten states
#BW0GEB_B == GR4JState1, BW3GEB_B == GR4JState2
state1_ind2 <- grep(names(states), pattern= "BW0GEB")
state2_ind2 <- grep(names(states), pattern= "BW3GEB")
names(states)[state1_ind2] <- gsub(names(states)[state1_ind2], 
                                   pattern = "BW0GEB", replacement = "GR4JState1")
names(states)[state2_ind2] <- gsub(names(states)[state2_ind2], 
                                   pattern = "BW3GEB", replacement = "GR4JState2")
# write to feather
feather::write_feather(states[, c(1:3, state1_ind2)], "True parameters/GR4J_state1.feather")
feather::write_feather(states[, c(1:3,state2_ind2)], "True parameters/GR4J_state2.feather")
# save for comparison
states1_2km <- states[, c(1:3, state1_ind2)]
states2_2km <- states[, c(1:3, state2_ind2)]


# 5. Create Synthetic parameters for 1km GR4J ------------------------------------------
new_gr4j_para_1km <- create_GR4J_para(transfer_functions, l0, parameter_bounds, 
                                      gmean_parameter = gmean_parameter, km1 = TRUE)

# merge new parameters with parameter file
para_new_1km <- merge(para_1km, new_gr4j_para_1km, by = c("NB_", "NZ_"), suffixes =  c(".old", ""),
                      all.x = TRUE)
para_new_1km[which(is.na(para_new_1km), arr.ind = TRUE)] <- 0 # parameter of unused basins -> 0
para_new_1km <- para_new_1km[, -grep(".old", names(para_new_1km))]
states_1km <- para_new_1km[, grep("iniSt", names(para_new_1km))]
para_new_1km <- para_new_1km[, -grep("iniSt", names(para_new_1km))]
para_new_1km <- cbind(para_new_1km, states_1km, stringsAsFactors = FALSE)
para_new_1km <- para_new_1km[order(para_new_1km$NZ_), ]

para_new_1km[, 1:3] <- para_new_1km[, c(1, 3, 2)]
names(para_new_1km)[2:3] <- c("IZ_", "NZ_")

cat("\n", file = "GR4J_distributed_1km/input/para_Mur_GR4J_true.txt")
write.table(para_new_1km, "GR4J_distributed_1km/input/para_Mur_GR4J_true.txt",
            append = TRUE, row.names = FALSE, quote = FALSE)

# 6. Q & states for 1km GR4J -----------------------------------------------------------

# save output as new Qobs file

# **************************  Change in defaults ******************************
# run GR4J_parallel for all basins with:
# para_Mur_GR4J_true.txt parameter file
# full time series: 2003 - 2012
# true Qobs
# all basins
# ************************************************************************
default <- readLines("GR4J_distributed_1km/input/defaults.txt")
# set para file
default[16] <- "para_Mur_GR4J_true.txt"
default[17] <- "para_Mur_GR4J_fsOptim.txt"
# set dates
default[31] <- "2003 01 02 00 00" # start date
default[32] <- "2009 08 31 00 00" # backup
default[35] <- "2012 12 31 00 00" # end date
default[36] <- "2009 08 31 00 00" # backup
# set spin-up time
default[39] <- "241"
default[40] <- "2433"
# set Datafile
default[7] <- "QObs_24h.txt" # used
default[8] <- "QObs_24h_synth_GR4J.txt" # backup
# set output type
default[25] <- "1"
# write new default file
writeLines(default, con = "GR4J_distributed_1km/input/defaults.txt")

# Start Model run
setwd("GR4J_distributed_1km/")
sys::exec_wait("start_GR4J_case_study_full_run.bat",
               std_out = "GR4J_output.txt")
setwd("..")

# Generate parameters and start model
synth_obs <- read.table("GR4J_distributed_1km/output/GR4J.runoff", skip = 17, header = TRUE)
synth_obs <- synth_obs[, -c(grep("QOBS", names(synth_obs)),
                            grep("Qloc", names(synth_obs)))]
synth_obs[which(is.na(synth_obs), arr.ind = TRUE)] <- -999
cat("synthetic_Qobs_Mur_24h_mean_UTC   2003-2012\n",
    file = "GR4J_distributed_1km/input/QObs_24h_synth_GR4J.txt")
for(i in 1:112){
  cat(paste0(gsub("QSIM", "QOBS", names(synth_obs)[5+i]),  "\n"),
      file = "GR4J_distributed_1km/input/QObs_24h_synth_GR4J.txt", append = TRUE)
}
cat("##############################################\n",
    file = "GR4J_distributed_1km/input/QObs_24h_synth_GR4J.txt", append = TRUE)

write.table(synth_obs, "GR4J_distributed_1km/input/QObs_24h_synth_GR4J.txt",
            append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

# Generate true states time series
# read state results
state_names <- read.table("GR4J_distributed_1km/output/COSERO.plus1", 
                          header=FALSE, fill=TRUE, skip = 1, nrows = 1, 
                          stringsAsFactors = FALSE)
state1_ind <- grep(state_names, pattern= "BW0GEB")
state2_ind <- grep(state_names, pattern= "BW3GEB")
state_classes <- rep("NULL", length(state_names))
state_classes[c(state1_ind, state2_ind)] <- "numeric"
state_classes[1:3] <- "integer"
states <- read.table("GR4J_distributed_1km/output/COSERO.plus1", 
                     header=FALSE, fill=FALSE, skip = 243,
                     colClasses = state_classes)
state_classes[c(1:3, state1_ind, state2_ind)] <- "character"

names(states) <- read.table("GR4J_distributed_1km/output/COSERO.plus1", 
                            header=FALSE, fill=FALSE, skip = 1,nrows = 1,
                            colClasses = state_classes)

# relevanten states
#BW0GEB_B == GR4JState1, BW3GEB_B == GR4JState2
state1_ind2 <- grep(names(states), pattern= "BW0GEB")
state2_ind2 <- grep(names(states), pattern= "BW3GEB")
names(states)[state1_ind2] <- gsub(names(states)[state1_ind2], 
                                   pattern = "BW0GEB", replacement = "GR4JState1")
names(states)[state2_ind2] <- gsub(names(states)[state2_ind2], 
                                   pattern = "BW3GEB", replacement = "IGR4JState2")
# write to feather
feather::write_feather(states[, c(1:3, state1_ind2)], "True parameters/GR4J_state1_1km.feather")
feather::write_feather(states[, c(1:3,state2_ind2)], "True parameters/GR4J_state2_1km.feather")

# save for comparison
states1_1km <- states[, c(1:3, state1_ind2)]
states2_1km <- states[, c(1:3, state2_ind2)]



# 7. Check iything was done correctly --------------------------------------------------

# Check if 1km and 2km are roughly the same
km2Q <- read.table( "GR4J_distributed/input/QObs_24h_synth_GR4J.txt", skip = 114, header = FALSE)
km1Q <- read.table("GR4J_distributed_1km/input/QObs_24h_synth_GR4J.txt", skip = 114, header = FALSE)


for(i in 1:112){
  cat("\nBasin", i, "\n")
  QNSE <- NSE(km2Q[-c(1:241), 5+i], km1Q[-c(1:241), 5+i])
  stNSE1 <-  NSE(states1_1km[, 3+i], states1_2km[, 3+i])
  stNSE2 <- NSE(states2_1km[, 3+i], states2_2km[, 3+i])
  if(sum(c(QNSE) < 0.97) > 0) {
    stop(paste0("The Q difference between 1km and 2km in basin ", i, " is too large!\nQ NSE = ", QNSE))
  }
  cat("NSE Q:", QNSE, "\n")
  cat("NSE State 1:", stNSE1, "\n")
  cat("NSE State 2:", stNSE2, "\n")
}


plot(km2Q[-c(1:241), 6], type = "l", ylab = "Q", xlab = "t") # Basin 1
plot(km2Q[-c(1:241), 105], type = "l", ylab = "Q", xlab = "t") # Basin 100
plot(km2Q[-c(1:241), 60], type = "l", ylab = "Q", xlab = "t") 
