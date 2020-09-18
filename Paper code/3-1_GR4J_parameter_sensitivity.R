# Estimate parameter sensitivity for d-GR4J
# Moritz Feigl, 2019
#

# 0. Setup -------------------------------------------------------------------------------
setwd("FSO_paper")
source("Functions/FSO_functions.R")
library(ggplot2)
FSO_setup()
# create folder
if (!file.exists("Basin infos")){
  dir.create("Basin infos")
}
# examined basins
basin_subset <- c(1, 100)
# load l0 layer
l0 <- load_sp_mur(scale = TRUE, na.approx = TRUE, 
                  only_training_basins = FALSE, full_dataset = FALSE)
para <- read.table("GR4J_distributed/input/para_Mur_GR4J_true.txt",
                   skip = 2, sep = " ", header = FALSE)
colnames(para) <- read.table("GR4J_distributed/input/para_Mur_GR4J_true.txt",
                             skip = 1, sep = " ", header = FALSE, nrows = 1, stringsAsFactors = FALSE)
para <- para[, names(para) != "NA"]
bounds <- read.table("GR4J_distributed/parameter.txt",
                     skip = 1, sep = "\t", header = FALSE)
parameter_bounds <- list(GR4Jx1 = as.numeric(bounds[7, c(3, 4)]), 
                         GR4Jx2 = as.numeric(bounds[8, c(3, 4)]),          
                         GR4Jx3 = as.numeric(bounds[9, c(3, 4)]),  
                         GR4Jx4 = as.numeric(bounds[10, c(3, 4)]))
# get real mean parameter values for basins
real_paras <- data.frame(Basins = basin_subset, x1 = NA, x3 = NA, x4 = NA)
for(basin in 1:length(basin_subset)){
  x1 <- mean(para[para$NB_ == basin_subset[basin], "GR4Jx1"])
  x3 <- mean(para[para$NB_ == basin_subset[basin], "GR4Jx3"])
  x4 <- mean(para[para$NB_ == basin_subset[basin], "GR4Jx4"])
  real_paras[basin, -1] <- c(x1, x3, x4)
}
print(real_paras)

# 1. Model setup -----------------------------------------------------------------------
# Parameter sensitivity for basin_subset: 1, 100
default <- readLines("GR4J_distributed/input/defaults.txt")
# set para file
default[14] <- "para_Mur_GR4J_MC.txt"
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
default[7] <- "QObs_24h_synth_GR4J.txt" # used
default[8] <- "QObs_24h.txt" # backup
# set output type
default[26] <- "1"
# write new default file
writeLines(default, con = "GR4J_distributed/input/defaults.txt")
# steps:
# 1. sample parameters for basin_subset
# 2. Run model
# 3. Get NSE
# 4. write results
paras_x1 <- eval(parse(text = paste0(parameter_bounds[["GR4Jx1"]], collapse = ":")))
paras_x3 <- eval(parse(text = paste0(parameter_bounds[["GR4Jx3"]], collapse = ":")))
paras_x4 <- eval(parse(text = paste0("((", parameter_bounds[["GR4Jx4"]][1], 
                                     "*100):(",
                                     parameter_bounds[["GR4Jx4"]][2], 
                                     "*100))/100")))
# create dir for saving results
dir.create(paste0("Basin infos/parameter sensitivity/sens_analysis_long"))

# 2. MC simulation ---------------------------------------------------------------------
for(i in 1624:5000){
  results <- data.frame(Basins = basin_subset, x1 = NA, x3 = NA, x4 = NA, NSE = NA)
  for(basin in 1:length(basin_subset)){
    x1 <- sample(paras_x1, 1)
    x3 <- sample(paras_x3, 1)
    x4 <- sample(paras_x4, 1)
    para[para$NB_ == basin_subset[basin], "GR4Jx1"] <- x1
    para[para$NB_ == basin_subset[basin], "GR4Jx3"] <- x3
    para[para$NB_ == basin_subset[basin], "GR4Jx4"] <- x4
    results[basin, -c(1, 5)] <- c(x1, x3, x4)
  }
  cat("\n", file = "GR4J_distributed/input/para_Mur_GR4J_MC.txt")
  write.table(para, "GR4J_distributed/input/para_Mur_GR4J_MC.txt", append = TRUE,
              row.names = FALSE, quote = FALSE)
  
  # Start Model run
  setwd("GR4J_distributed/")
  model_run <- try(sys::exec_wait("start_GR4J_case_study_MC.bat",
                                  std_out = "GR4J_output.txt"))
  while(class(model_run) == "try-error"){
    model_run <- try(sys::exec_wait("start_GR4J_case_study_MC.bat",
                                    std_out = "GR4J_output.txt"))
  }
  setwd("..")
  # Get statistics and calculate losses
  statistics <- read.table("GR4J_distributed/output/statistics_gr4j_Mur.txt",
                           skip = 21, header = TRUE)
  results$NSE <- statistics$NSE
  cat("\nrun", i, "\n")
  print(results)
  if(i == 1){
    all_results <- results
  } else {
    all_results <- rbind(all_results, results)
  }
  # Every 100 iterations produce set of sensitivity plots
  if((i %% 100) == 0){
    setwd(paste0("Basin infos/parameter sensitivity/sens_analysis"))
    for(basin in c(1, 100)){
      ggplot(all_results[all_results$Basins == basin,], aes(x = x1, y = NSE)) + 
        geom_point() + ylim(0, 1) + 
        geom_vline(xintercept = real_paras[real_paras$Basins == basin, "x1"]) +
        annotate(geom = "text", x =  real_paras[real_paras$Basins == basin, "x1"] + 290, 
                 y = 0., label = "real mean basin value") +
        labs(x = "X1", title = paste0("Basin ", basin, " - Parameter sensitivity X1")) +
        ggsave(paste0("Basin_", basin, "_X1_sensitivity.png"), height = 7, width = 7, unit = "in")
      
      ggplot(all_results[all_results$Basins == basin,], aes(x = x3, y = NSE)) + 
        geom_point() + ylim(0, 1) + 
        geom_vline(xintercept = real_paras[real_paras$Basins == basin, "x3"]) +
        annotate(geom = "text", x =  real_paras[real_paras$Basins == basin, "x3"] + 70, 
                 y = 0., label = "real mean basin value") + 
        labs(x = "X3", title = paste0("Basin ", basin, " - Parameter sensitivity X3"))
      ggsave(paste0("Basin_", basin, "_X3_sensitivity.png"), height = 7, width = 7, unit = "in")
      
      ggplot(all_results[all_results$Basins == basin,], aes(x = x4, y = NSE)) +
        geom_point() + ylim(0, 1) + 
        geom_vline(xintercept = real_paras[real_paras$Basins == basin, "x4"]) +
        annotate(geom = "text", x =  real_paras[real_paras$Basins == basin, "x4"] + 0.84, 
                 y = 0., label = "real mean basin value") +  
        labs(x = "X4", title = paste0("Basin ", basin, " - Parameter sensitivity X4"))
      ggsave(paste0("Basin_", basin, "_X4_sensitivity.png"), height = 7, width = 7, unit = "in")
      
      ggplot(all_results[all_results$Basins == basin & all_results$NSE > 0,], 
             aes(x = x1, y = x3, col = NSE)) + 
        geom_point() +
        scale_colour_gradientn(colors = c("red", "deeppink2", 
                                          "darkolivegreen1", "darkolivegreen"),
                               limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
        labs(x = "X1", title = paste0("Basin ", 1, " - Parameter sensitivity X1/X3")) +
        geom_point(aes(x = real_paras$x1[real_paras$Basins == basin], 
                       y = real_paras$x3[real_paras$Basins == basin]), 
                   col = "black", size = 5) +
        ggsave(paste0("Basin_", basin, "_NSE_field.png"), height = 4, width = 9, unit = "in")
    }
    feather::write_feather(all_results, "MC_sensitivity_analysis.feather")
    setwd("../../..")
  }
}

# 3. FAST sensitivity ------------------------------------------------------------------
# sensitivity analysis with Fourier Amplitude Sensitivity Test (FAST)
library(fast)
library(dplyr)
library(ggplot2)
library(tidyverse)
par_fast <- fast_parameters(
  minimum = c(1, 1, 0.5),
  maximum = c(4000, 1000, 6),
  names = c("X1", "X3", "X4")
) %>% as_tibble()
for(i in 1:nrow(par_fast)){
  results <- data.frame(Basins = basin_subset, x1 = NA, x3 = NA, x4 = NA, NSE = NA)
  for(basin in 1:length(basin_subset)){
    x1 <- par_fast[i, 1]
    x3 <- par_fast[i, 2]
    x4 <- par_fast[i, 3]
    para[para$NB_ == basin_subset[basin], "GR4Jx1"] <- x1
    para[para$NB_ == basin_subset[basin], "GR4Jx3"] <- x3
    para[para$NB_ == basin_subset[basin], "GR4Jx4"] <- x4
    results[basin, -c(1, 5)] <- c(x1, x3, x4)
  }
  cat("\n", file = "GR4J_distributed/input/para_Mur_GR4J_MC.txt")
  write.table(para, "GR4J_distributed/input/para_Mur_GR4J_MC.txt", append = TRUE,
              row.names = FALSE, quote = FALSE)
  # Start Model run
  setwd("GR4J_distributed/")
  model_run <- try(sys::exec_wait("start_GR4J_case_study_MC.bat",
                                  std_out = "GR4J_output.txt"))
  while(class(model_run) == "try-error"){
    model_run <- try(sys::exec_wait("start_GR4J_case_study_MC.bat",
                                    std_out = "GR4J_output.txt"))
  }
  setwd("..")
  # Get statistics and calculate losses
  statistics <- read.table("GR4J_distributed/output/statistics_gr4j_Mur.txt",
                           skip = 21, header = TRUE)
  results$NSE <- statistics$NSE
  cat("\nrun", i, "\n")
  print(results)
  if(i == 1){
    all_results <- results
  } else {
    all_results <- rbind(all_results, results)
  }
}
basin1_nse <- all_results[all_results$Basins== 100, "NSE"]
sens_fast <- sensitivity(basin1_nse, 3)
result_fast <- tibble(parameter = names(par_fast),
                      fast      = sens_fast) %>%
  mutate(parameter = factor(parameter) %>% fct_reorder(., fast))
ggplot(data = result_fast) +
  geom_col(aes(x = parameter, y = fast)) +
  xlab("Parameter") +
  ylab("Sensitivity") +
  coord_flip() +
  theme_bw() +
  ggsave("Basin infos/parameter sensitivity/fast_parameter_sensitivity.png", 
         width = 7, height = 7, units = "in")

