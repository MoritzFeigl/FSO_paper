# Select synthetic transfer functions
# Moritz Feigl, 2019
#

setwd("FSO_paper")

if (!file.exists("True parameters")){
  dir.create("True parameters")
}
if (!file.exists("True parameters/Plots")){
  dir.create("True parameters/Plots")
}

# 1. Load Data and define help functions -----------------------------------------------
# get True functions
source("Functions/case_study_setup.R")
# Load helper functions
source("Functions/FSO_functions.R")

# load l0 layer
l0 <- load_sp_mur(scale = TRUE, na.approx = TRUE, 
                  only_training_basins = FALSE, full_dataset = FALSE)

bounds <- read.table("GR4J_distributed/parameter.txt",
                     skip = 1, sep = "\t", header = FALSE)
parameter_bounds <- list(GR4Jx1 = as.numeric(bounds[7, c(3, 4)]), 
                         GR4Jx2 = as.numeric(bounds[8, c(3, 4)]),          
                         GR4Jx3 = as.numeric(bounds[9, c(3, 4)]),  
                         GR4Jx4 = as.numeric(bounds[10, c(3, 4)]))
transfer_functions <- list(GR4Jx1 = "-1.5 + evi*1.5 + exp(bdim)*0.9", # production store
                           GR4Jx2 = "0",          # interaction direct runoff and routing store
                           GR4Jx3 = "-1.3  - log(slope)",  # routing store
                           GR4Jx4 = "-log(hand)*0.2 - slope*1.5 - 1.5") # length of UH

# 2. Analyze chosen transfer functions -------------------------------------------------
# check the resulting parameter distributions
# X1 production store: interpretation as soil water storage -> evi and bdim
summary(
  rescale(evaluate_function_from_string(transfer_functions[["GR4Jx1"]], l0 = l0),
          to = parameter_bounds[["GR4Jx1"]])
)
# X3 routing store: slope and hand define how water is routed in catchment
summary(
  rescale(evaluate_function_from_string(transfer_functions[["GR4Jx3"]], l0 = l0), 
          to = parameter_bounds[["GR4Jx3"]])
)
# X4 unit hydrograph length
summary(
  rescale(evaluate_function_from_string(transfer_functions[["GR4Jx4"]], l0 = l0), 
          to = parameter_bounds[["GR4Jx4"]])
)

# 3. Search TFs in VAE data ------------------------------------------------------------
# check if there is something resembling this kind of function in the training/test data
transfer_functions <- feather::read_feather("../3_Function space/Data/functions_simple_onlyfunctions.feather")
tfs <- sapply(transfer_functions[, 2], gsub, pattern = " ", replacement = "")
# p1: 0.5 + exp(bdim)* 0.9 + evi * 1.5
p1_search <- grep("numeric*exp(bdim)", tfs, fixed = TRUE)
tfs_p1 <- tfs[p1_search]
p1_search2 <- grep("-evi*numeric", tfs_p1, fixed = TRUE)
tfs_p1[p1_search2]
# p3: -1.3  - log(slope)
p3_search <- grep("numeric+log(slope)", tfs, fixed = TRUE)
tfs_p3 <- tfs[p3_search]
chars <- sapply(tfs_p3, nchar)
p3_search2 <- grep("numeric+", tfs_p3, fixed = TRUE)
tfs_p3_2 <- tfs_p3[p3_search2]
p3_search3 <- grep("slope^2", tfs_p3_2, fixed = TRUE)
tfs_p3_2[p3_search3]
# p4: -1.5 - log(hand) * 0.2 - slope * 1.5
p4_search <- grep("numeric*log(hand)", tfs, fixed = TRUE)
tfs_p4 <- tfs[p4_search]
p4_search2 <- grep("-numeric*slope", tfs_p4, fixed = TRUE)
tfs_p4_2 <- tfs_p4[p4_search2]
p4_search3 <- grep("+numeric", tfs_p4_2, fixed = TRUE)
tfs_p4_2[p4_search3]

# 4. Check basins ----------------------------------------------------------------------
# Check if training basins cover all characteristics of catchment
l0$'train_test' <- 'test'
l0$'train_test'[l0$nb %in% training_basins] <- 'train'
library(ggplot2)
library(ggpubr)
library(gridExtra)
# Elevation
l0_agg <- aggregate(elevation ~ nb, l0, mean)
l0_agg <- merge(l0_agg, l0[, c("nb", "train_test")], all.x = TRUE)
p1 <- ggplot(l0_agg, aes(factor(train_test), elevation)) + geom_boxplot() + ylab("elevation") + xlab("")
# slope
l0_agg <- aggregate(slope ~ nb, l0, mean)
l0_agg <- merge(l0_agg, l0[, c("nb", "train_test")], all.x = TRUE)
p2 <- ggplot(l0_agg, aes(factor(train_test), slope)) + geom_boxplot() + ylab("slope") + xlab("")
# bdim
l0_agg <- aggregate(bdim ~ nb, l0, mean)
l0_agg <- merge(l0_agg, l0[, c("nb", "train_test")], all.x = TRUE)
p3 <- ggplot(l0_agg, aes(factor(train_test), bdim)) + geom_boxplot() + ylab("soil depth") + xlab("")
# sand
l0_agg <- aggregate(sand ~ nb, l0, mean)
l0_agg <- merge(l0_agg, l0[, c("nb", "train_test")], all.x = TRUE)
p4 <- ggplot(l0_agg, aes(factor(train_test), sand)) + geom_boxplot() + ylab("sand") + xlab("")
# evi
l0_agg <- aggregate(evi ~ nb, l0, mean)
l0_agg <- merge(l0_agg, l0[, c("nb", "train_test")], all.x = TRUE)
p5 <- ggplot(l0_agg, aes(factor(train_test), evi)) + geom_boxplot() + ylab("EVI") + xlab("")
# hand
l0_agg <- aggregate(hand ~ nb, l0, mean)
l0_agg <- merge(l0_agg, l0[, c("nb", "train_test")], all.x = TRUE)
p6 <- ggplot(l0_agg, aes(factor(train_test), hand)) + geom_boxplot() + ylab("hand") + xlab("")

ggsave(file = "Basin infos/train_test_boxplots.png",
       arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2), 
       height = 18, width = 30, unit = "cm")
