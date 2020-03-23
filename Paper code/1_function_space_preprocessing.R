#
# Moritz Feigl, Jan 2019
#
# 
#

setwd("FSO_paper")
library(feather)
library(tidyverse)
source("Functions/CFG_functions")


# 1. Create Grammar --------------------------------------------------------------------
g <- grammar(tf = "numeric * <eq> + numeric, <eq> + numeric, <eq> ",
             eq = "<fs>, <eq><op><fs>, <eq><op>numeric, <fs><op>numeric",
             fs = "<sp>,  <f>(<sp>), <sp><op><sp>, numeric",
             f = "exp, log",
             op = "+, -, *, /",
             sp = "slope, evi, sand, clay, elevation, hand, noise, bdim")
spatial_predictor_variables <- c("slope", "evi", "sand", "clay", "elevation", "hand", "noise", "bdim")

# 2. Random sample grammar & simplify --------------------------------------------------
function_df2 <- par.grammar.sampler(n = 5000000,
                                    cfgram = g,
                                    max.depth = 3,
                                    save_feather = FALSE,
                                    parallel = TRUE,
                                    no_cores = 18)
# take only unique functions
function_df2 <- function_df2[!duplicated(function_df2$Transfer_Function), ]
# save function_list as feather
write_feather(function_df, "functions_V2.feather")
function_df2$Transfer_Function_simple <- parlapply_simplify(funs = function_df2$Transfer_Function,
                                                            function_variables = spatial_predictor_variables)
# remove NA functions
function_df2 <- function_df2[!is.na(function_df2$Transfer_Function_simple), ]
write_feather(function_df2, "Data/functions_simple_onlyfunctions.feather")

# 3. Put in numeric values from [-1.5, 1.5] ----------------------------------------------
numbers <- c(-1*seq(0.1, 1.5, 0.1), seq(0.1, 1.5, 0.1))
# Functions for randomly putting in numeric values from given range
num_input <- function(fun){
  while(length(grep("numeric", fun)) != 0){
    fun <- sub("numeric", sample(numbers, 1), fun)
  }
  return(fun)
}
nums_input <- function(fun){
  fun <- gsub(" ", "", fun)
  funs <- replicate(10, num_input(fun))
  funs <- gsub("+-", "-", funs, fixed = TRUE)
  funs <- gsub("--", "+", funs, fixed = TRUE)
  funs <- gsub("++", "+", funs, fixed = TRUE)
  funs <- gsub("-+", "-", funs, fixed = TRUE)
  return(funs)
}
num_df_input <- function(functions, numbers){
  number_list <- lapply(functions, nums_input)
  return(do.call(c, number_list))
}
# Create new function vector with 10 numeric inputs each
functions_simple_10_numerics <- num_df_input(function_df2$Transfer_Function_simple, numbers)
# take only unique functions
functions_simple_10_numerics_unique <- functions_simple_10_numerics$TF_simple_10numerics[
  !duplicated(functions_simple_10_numerics$TF_simple_10numerics)
  ]

# remove exotic functions ----------------------------------------------------------------
# * 6
all_id6 <- grep("6*", functions_simple_10_numerics_unique, fixed = TRUE)
id6 <- all_id6[!(all_id6 %in% grep(".6*", functions_simple_10_numerics_unique, fixed = TRUE))]
functions_simple_10_numerics_unique <- functions_simple_10_numerics_unique[-id6]
# * 5
all_id5 <- grep("5*", functions_simple_10_numerics_unique, fixed = TRUE)
id5 <- all_id5[!(all_id5 %in% grep(".5*", functions_simple_10_numerics_unique, fixed = TRUE))]
functions_simple_10_numerics_unique <- functions_simple_10_numerics_unique[-id5]
# ^5
idp5 <- grep("^5", functions_simple_10_numerics_unique, fixed = TRUE)
functions_simple_10_numerics_unique <- functions_simple_10_numerics_unique[-idp5]
# ^4
idp4 <- grep("^4", functions_simple_10_numerics_unique, fixed = TRUE)
functions_simple_10_numerics_unique <- functions_simple_10_numerics_unique[-idp4]
# ^3
idp3 <- grep("^3", functions_simple_10_numerics_unique, fixed = TRUE)
functions_simple_10_numerics_unique <- functions_simple_10_numerics_unique[-idp3]
# * 3
all_id3 <- grep("3*", functions_simple_10_numerics_unique, fixed = TRUE)
id3 <- all_id3[!(all_id3 %in% grep(".3*", functions_simple_10_numerics_unique, fixed = TRUE))]
functions_simple_10_numerics_unique <- functions_simple_10_numerics_unique[-id3]
# * 4
all_id4 <- grep("4/*", functions_simple_10_numerics_unique, fixed = FALSE)
id4 <- all_id4[!(all_id4 %in% grep(".4", functions_simple_10_numerics_unique, fixed = TRUE))]
functions_simple_10_numerics_unique <- functions_simple_10_numerics_unique[-id4]
id42 <- grep("(4", functions_simple_10_numerics_unique, fixed = TRUE)
functions_simple_10_numerics_unique <- functions_simple_10_numerics_unique[-id42]
id43 <- grep("-4", functions_simple_10_numerics_unique, fixed = TRUE)
functions_simple_10_numerics_unique <- functions_simple_10_numerics_unique[-id43]
id44 <- grep("/4", functions_simple_10_numerics_unique, fixed = TRUE)
functions_simple_10_numerics_unique <- functions_simple_10_numerics_unique[-id43]

# save
 write_feather(data.frame("TF_simple_10numerics" = functions_simple_10_numerics_unique, stringsAsFactors = FALSE),
 "Data/functions_simple_10_numerics_pp.feather")

# 4. Distribution properties -------------------------------------------------------------
functions_simple_10_numerics <- read_feather("Data/functions_simple_10_numerics_pp.feather.feather")
# Load spatial predictors
source("Functions/FSO_functions.R")
spatial_predictors <- load_sp_mur(scale = TRUE, na.approx = FALSE, only_training_basins = FALSE)

evaluate_function_from_string <- function(string, l0){
  # evaluates a function given as a string for all relevant l0 layer
  # Input:
  #    string: transfer function as string
  #    l0: data frame with l0 layer
  # Output:
  #    vector with all evaluated function outputs
  tf <- unlist(strsplit(string, c("[/^()*+-]")))
  tf <- gsub(" ", "", tf, fixed = TRUE)
  relevant_predictors <- l0[which(names(l0) %in% tf)]
  if(ncol(relevant_predictors) == 0){
    args <- ""
    eval(parse(text = paste('f <- function(', args, ') { return(' , string , ')}', sep='')))
    f_evaluated <- rep(f(), nrow(l0))
  } else {
    args <- paste(names(relevant_predictors), collapse = ', ')

    eval(parse(text = paste('f <- function(', args, ') { return(' , string , ')}', sep='')))
    f_evaluated <- eval(parse(text = paste('mapply(f, ',
                                           paste0('relevant_predictors$',
                                                  names(relevant_predictors), collapse = ', '),
                                           ')')))
  }
  return(f_evaluated)
}
# function to create distribution parameters of a given TF
distribution_values_from_tf <- function(tf, spatial_predictors, cut_off = NULL){
  tf_evaluated <- evaluate_function_from_string(tf, l0 = spatial_predictors)
  if(!is.null(cut_off)){
    tf_evaluated[tf_evaluated < cut_off[1]] <- cut_off[1]
    tf_evaluated[tf_evaluated > cut_off[2]] <- cut_off[2]
  }
  distribution_values <- round(c(
    min(tf_evaluated, na.rm = TRUE),
    quantile(tf_evaluated, probs = c(0.1, 0.2, 0.3, 0.4), na.rm = TRUE),
    mean(tf_evaluated, na.rm = TRUE),
    quantile(tf_evaluated, probs = c(0.6, 0.7, 0.8, 0.9), na.rm = TRUE),
    max(tf_evaluated, na.rm = TRUE)
    ), 4)
  names(distribution_values) <- c("min", "10%", "20%", "30%", "40%", "mean",
                                  "60%", "70%", "80%", "90%", "max")
  return(distribution_values)
}

tf_distributions <- parallel::mclapply(functions_simple_10_numerics$TF_simple_10numerics,
                           distribution_values_from_tf, spatial_predictors = spatial_predictors,
                           mc.cores = 18,
                           mc.set.seed = TRUE)
tf_distributions_df <- do.call(rbind, tf_distributions)
tf_distributions_df <- data.frame(transfer_function = functions_simple_10_numerics$TF_simple_10numerics,
                                  tf_distributions_df,
                                  stringsAsFactors = FALSE)
write_feather(tf_distributions_df, "Data/functions_simple_10_numerics_Distribution_indiv_scale_allBasins.feather")

# Cut off distributions ------------------------------------------------------------------
generator_data <- read_feather("Data/generator_data_simple_10numerics.feather")
tf_dist <- read_feather("Data/functions_simple_10_numerics_Distribution_indiv_scale_allBasins.feather")
# remove functions with NA mean
na_ids <- is.na(tf_dist$mean)
generator_data <- generator_data[!na_ids, ]
tf_dist <- tf_dist[!na_ids, ]
# remove functions with Nan or Inf
inf_ids <- is.infinite(tf_dist$mean)
generator_data <- generator_data[!inf_ids, ]
tf_dist <- tf_dist[!inf_ids, ]
# Check in which interval 80% of the functions are located
threshold <- 11
or_ind <- abs(tf_dist$min) > threshold | abs(tf_dist$max) > threshold
sum(or_ind)/nrow(tf_dist)
# create dfs with all out of range distribution -> recalculate distributions with threshold
tf_dist_recalc <- tf_dist[or_ind, ]
generator_data_recalc <- generator_data[or_ind, ]
# remove out of range distributions from 
tf_dist <- tf_dist[!or_ind, ]
generator_data <- generator_data[!or_ind, ]
# remove distributions where all values are either below or above the range
out_of_range <- tf_dist_recalc$max > threshold & tf_dist_recalc$min > threshold | 
  tf_dist_recalc$max < -threshold & tf_dist_recalc$min < -threshold
tf_dist_recalc <- tf_dist_recalc[!out_of_range, ]
generator_data_recalc <- generator_data_recalc[!out_of_range, ]
# recalculate distributions values with threshold
tf_distributions_recalc1 <- parallel::mclapply(tf_dist_recalc$transfer_function,
                                    distribution_values_from_tf, 
                                    spatial_predictors = spatial_predictors, 
                                    cut_off = c(-threshold, threshold),
                                    mc.cores = 18)
# fix null values created by mclapply
if(sum(sapply(tf_distributions_recalc1, is.null)) > 0){
  id_nulls <- which(sapply(tf_distributions_recalc1, is.null))
  tf_nulls <- tf_dist_recalc$transfer_function[id_nulls]
  library(parallel)
  cl <- makeCluster(mc <- getOption("cl.cores", 18))
  clusterExport(cl=cl, varlist=c("evaluate_function_from_string"))
  dist_nulls <- parLapply(cl, tf_nulls, distribution_values_from_tf,
                    spatial_predictors = spatial_predictors, 
                    cut_off = c(-threshold, threshold))
}
# fix them
for(i in seq_along(id_nulls)){
  tf_distributions_recalc1[id_nulls[i]] <- dist_nulls[i]
}
tf_distributions_recalc <- do.call(rbind, tf_distributions_recalc1)
tf_distributions_recalc <- data.frame(transfer_function = tf_dist_recalc$transfer_function,
                                      tf_distributions_recalc,
                                  stringsAsFactors = FALSE)
write_feather(rbind(tf_dist, tf_distributions_recalc), 
              "Data/functions_simple_10_numerics_Distribution_indiv_scale_wrecalc_allBasins.feather")
write_feather(rbind(generator_data, generator_data_recalc), 
              "Data/generator_data_simple_10numerics_wrecalc_allBasins.feather")
# remove all functions that are mainly produce the max and min values
tf_dist <- read_feather("Data/functions_simple_10_numerics_Distribution_indiv_scale_wrecalc_allBasins.feather")
generator_data <- read_feather("Data/generator_data_simple_10numerics_wrecalc_allBasins.feather")
extreme_tfs_ind <- tf_dist$min == -11 & tf_dist$X10. == -11 & tf_dist$X20. == -11 |
  tf_dist$X80. == 11 & tf_dist$X90. == 11 & tf_dist$max == 11
tf_dist <- tf_dist[!extreme_tfs_ind, ]
generator_data <- generator_data[!extreme_tfs_ind, ]
write_feather(tf_dist, 
              "Data/functions_simple_10_numerics_Distribution_indiv_scale_wrecalc_allBasins_no_extremes.feather")
write_feather(generator_data,  
              "Data/generator_data_simple_10numerics_wrecalc_allBasins_no_extremes.feather")




