#
# Optimization in function space
#

# 1. Prepare data and parameters -------------------------------------------------------
if(!exists("spatial_predictors")) {
  spatial_predictors <- load_sp_mur(scale = TRUE, 
                                    na.approx = TRUE, 
                                    only_training_basins = TRUE,
                                    training_basins = training_basins)
}
# get parameter table and bounds
source("Functions/parameter_and_para_bounds_GR4J.R")
# Get true parameter Field
true_para_field_df <- true_para_field(Test_number)
# 4. Objective Function ------------------------------------------------------------------
objective_function <- function(point, Test_number, 
                               true_para_field_df, spatial_predictors, parameter_bounds,
                               para, para_1km, training_basins){
  latent_dim <- 6
  # Make a counter for the optimization
  counter <- get("counter", envir = .GlobalEnv)
  counter <- counter + 1
  assign("counter", counter, envir = .GlobalEnv)
  cat("\n****** Point nr.", paste0(counter, "/", iterations), "******\n")
  
  # calculate the transfer function of the point
  point_tf <- vector(mode = "list")
  point_tf[[1]] <- tf_generator(matrix(as.numeric(point[1:6]), ncol = latent_dim))
  point_tf[[2]] <- tf_generator(matrix(as.numeric(point[7:12]), ncol = latent_dim))
  point_tf[[3]] <- tf_generator(matrix(as.numeric(point[13:18]), ncol = latent_dim))
  names(point_tf) <- c("GR4Jx1", "GR4Jx3", "GR4Jx4")
  cat("x1 = ", point_tf[[1]], "\n")
  cat("x2 = 0\n")
  cat("x3 = ", point_tf[[2]], "\n")
  cat("x4 = ", point_tf[[3]], "\n")
  function_eval <- try({
    lapply(point_tf, tf_evaluation)
  }, silent = TRUE)
  # Stop at empty or wrong functions
  if("" %in% unlist(point_tf) | "try-error" %in% sapply(function_eval, class)) {
    cat("No valid function found!\nResulting overall loss: NA\n")
    result_tracker_df <- get("result_tracker_df", envir = .GlobalEnv)
    assign("result_tracker_df",
           rbind(result_tracker_df,
                 data.frame(best_x1 = result_tracker_df[nrow(result_tracker_df), "best_x1"],
                            best_x3 = result_tracker_df[nrow(result_tracker_df), "best_x3"],
                            best_x4 = result_tracker_df[nrow(result_tracker_df), "best_x4"],
                            full_loss = result_tracker_df[nrow(result_tracker_df), "full_loss"],
                            NSE = result_tracker_df[nrow(result_tracker_df), "NSE"],
                            wNSE = result_tracker_df[nrow(result_tracker_df), "wNSE"],
                            x1 = point_tf[[1]],
                            x2 = point_tf[[2]],
                            x3 = point_tf[[3]],
                            n_iteration_used = 0,
                            n_iterations_since_BF_change = result_tracker_df[nrow(result_tracker_df), "n_iterations_since_BF_change"] + 1,
                            stringsAsFactors = FALSE)),
           envir = .GlobalEnv)
    return(-9999)
  }
  # check if always a spatial predictor is in any function
  functions_splitted <- lapply(point_tf, function_splitter)
  spatial_predictor_check <- sum(sapply(functions_splitted, function(x) sum(x %in% names(spatial_predictors)) < 1))
  
  if(spatial_predictor_check > 0){
    cat("No valid function found!\n")
    result_tracker_df <- get("result_tracker_df", envir = .GlobalEnv)
    assign("result_tracker_df",
           rbind(result_tracker_df,
                 data.frame(best_x1 = result_tracker_df[nrow(result_tracker_df), "best_x1"],
                            best_x3 = result_tracker_df[nrow(result_tracker_df), "best_x3"],
                            best_x4 = result_tracker_df[nrow(result_tracker_df), "best_x4"],
                            full_loss = result_tracker_df[nrow(result_tracker_df), "full_loss"],
                            NSE = result_tracker_df[nrow(result_tracker_df), "NSE"],
                            wNSE = result_tracker_df[nrow(result_tracker_df), "wNSE"],
                            x1 = point_tf[[1]],
                            x2 = point_tf[[2]],
                            x3 = point_tf[[3]],
                            n_iteration_used = 0,
                            n_iterations_since_BF_change = result_tracker_df[nrow(result_tracker_df), "n_iterations_since_BF_change"] + 1,
                            stringsAsFactors = FALSE)),
           envir = .GlobalEnv)
    return(-9999)
  }
  # model size loss
  # 0.01 divided by number of parameters
  model_size_loss <- size_loss(functions_splitted)
  # Function evaluation
  point_tf$GR4Jx2 <- "0"
  new_gr4j_para <- try(create_GR4J_para(transfer_functions = point_tf,
                                        l0 = spatial_predictors,
                                        parameter_bounds = parameter_bounds,
                                        gmean_parameter = gmean_parameter),
                       silent = TRUE)
  
  if(class(new_gr4j_para) == "try-error" ) {
    cat("No valid function found!\n")
    result_tracker_df <- get("result_tracker_df", envir = .GlobalEnv)
    assign("result_tracker_df",
           rbind(result_tracker_df,
                 data.frame(best_x1 = result_tracker_df[nrow(result_tracker_df), "best_x1"],
                            best_x3 = result_tracker_df[nrow(result_tracker_df), "best_x3"],
                            best_x4 = result_tracker_df[nrow(result_tracker_df), "best_x4"],
                            full_loss = result_tracker_df[nrow(result_tracker_df), "full_loss"],
                            NSE = result_tracker_df[nrow(result_tracker_df), "NSE"],
                            wNSE = result_tracker_df[nrow(result_tracker_df), "wNSE"],
                            x1 = point_tf[[1]],
                            x2 = point_tf[[2]],
                            x3 = point_tf[[3]],
                            n_iteration_used = 0,
                            n_iterations_since_BF_change = result_tracker_df[nrow(result_tracker_df), "n_iterations_since_BF_change"] + 1,
                            stringsAsFactors = FALSE)),
           envir = .GlobalEnv)
    return(-9999)
  }
  # Catch if NAs are produced
  if(sum(is.na(new_gr4j_para)) > 0){
    cat("Function produces too many NA NaN values!\n")
    result_tracker_df <- get("result_tracker_df", envir = .GlobalEnv)
    assign("result_tracker_df",
           rbind(result_tracker_df,
                 data.frame(best_x1 = result_tracker_df[nrow(result_tracker_df), "best_x1"],
                            best_x3 = result_tracker_df[nrow(result_tracker_df), "best_x3"],
                            best_x4 = result_tracker_df[nrow(result_tracker_df), "best_x4"],
                            full_loss = result_tracker_df[nrow(result_tracker_df), "full_loss"],
                            NSE = result_tracker_df[nrow(result_tracker_df), "NSE"],
                            wNSE = result_tracker_df[nrow(result_tracker_df), "wNSE"],
                            x1 = point_tf[[1]],
                            x2 = point_tf[[2]],
                            x3 = point_tf[[3]],
                            n_iteration_used = 0,
                            n_iterations_since_BF_change = result_tracker_df[nrow(result_tracker_df), "n_iterations_since_BF_change"] + 1,
                            stringsAsFactors = FALSE)),
           envir = .GlobalEnv)
    return(-9999)
  }
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
  cat("\n", file = "GR4J_distributed/input/para_Mur_GR4J_fsOptim.txt")
  suppressWarnings(write.table(para_new, "GR4J_distributed/input/para_Mur_GR4J_fsOptim.txt",
              append = TRUE, row.names = FALSE, quote = FALSE))
  # start model run
  setwd("GR4J_distributed/")
  sys::exec_wait("start_GR4J_case_study.bat",
                 std_out = "GR4J_output.txt")
  setwd("..")
  # Get statistics and calculate losses
  statistics <- read.table("GR4J_distributed/output/statistics_gr4j_Mur.txt",
                           skip = 21, header = TRUE)
  if(is.factor(statistics$NSE)){
    cat("Function produces no valid model outputs\n")
    result_tracker_df <- get("result_tracker_df", envir = .GlobalEnv)
    assign("result_tracker_df",
           rbind(result_tracker_df,
                 data.frame(best_x1 = result_tracker_df[nrow(result_tracker_df), "best_x1"],
                            best_x3 = result_tracker_df[nrow(result_tracker_df), "best_x3"],
                            best_x4 = result_tracker_df[nrow(result_tracker_df), "best_x4"],
                            full_loss = result_tracker_df[nrow(result_tracker_df), "full_loss"],
                            NSE = result_tracker_df[nrow(result_tracker_df), "NSE"],
                            wNSE = result_tracker_df[nrow(result_tracker_df), "wNSE"],
                            x1 = point_tf[[1]],
                            x2 = point_tf[[2]],
                            x3 = point_tf[[3]],
                            n_iteration_used = 0,
                            n_iterations_since_BF_change = result_tracker_df[nrow(result_tracker_df), "n_iterations_since_BF_change"] + 1,
                            stringsAsFactors = FALSE)),
           envir = .GlobalEnv)
    return(-9999)
  }
  if(Test_number %in% c(3.1, 3.2, 3.3)){
    # 1 KM Model parameters
    new_gr4j_para_1km <- create_GR4J_para(transfer_functions = point_tf, 
                                          l0 = spatial_predictors, 
                                          parameter_bounds = parameter_bounds, 
                                          gmean_parameter = gmean_parameter, 
                                          km1 = TRUE)
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
    cat("\n", file = "GR4J_distributed_1km/input/para_Mur_GR4J_fsOptim.txt")
    suppressWarnings(write.table(para_new, "GR4J_distributed_1km/input/para_Mur_GR4J_fsOptim.txt",
                append = TRUE, row.names = FALSE, quote = FALSE))
    # 1 KM Model run
    setwd("GR4J_distributed_1km/")
    sys::exec_wait("start_GR4J_case_study.bat",
                   std_out = "GR4J_output.txt")
    setwd("..")
    # Get statistics and calculate losses
    statistics_1km <- read.table("GR4J_distributed_1km/output/statistics_gr4j_Mur.txt",
                                 skip = 21, header = TRUE)
  } else {
    statistics_1km <- NULL
  }
  # calculate mean loss for all basins
  evaluation <- GR4J_model_quality(statistics = statistics,
                                   Test_number = Test_number,
                                   true_para_field_df = true_para_field_df, 
                                   model_size_loss = model_size_loss,
                                   new_gr4j_para = new_gr4j_para,
                                   statistics_1km = statistics_1km,
                                   relevant_basins = training_basins)
  mean_NSE <- evaluation[["mean_NSE"]]
  wmean_NSE <- evaluation[["wmean_NSE"]]
  full_loss <- evaluation[["full_loss"]]
  model_loss <- evaluation[["model_loss"]]
  if(is.na(full_loss)) full_loss <- -9999
  # Track results in an external data.frame
  result_tracker_df <- get("result_tracker_df", envir = .GlobalEnv)
  old_best <- result_tracker_df$full_loss[nrow(result_tracker_df)]
  if(!is.na(full_loss)){
    if(full_loss > old_best) {
      assign("result_tracker_df", rbind(result_tracker_df,
                                        data.frame(best_x1 = point_tf[[1]],
                                                   best_x3 = point_tf[[2]],
                                                   best_x4 = point_tf[[3]],
                                                   full_loss = full_loss,
                                                   NSE = mean_NSE,
                                                   wNSE = wmean_NSE,
                                                   x1 = point_tf[[1]],
                                                   x2 = point_tf[[2]],
                                                   x3 = point_tf[[3]],
                                                   n_iteration_used = 1,
                                                   n_iterations_since_BF_change = 1,
                                                   stringsAsFactors = FALSE)),
             envir = .GlobalEnv)
    } else {
      assign("result_tracker_df",
             rbind(result_tracker_df,
                   data.frame(best_x1 = result_tracker_df[nrow(result_tracker_df), "best_x1"],
                              best_x3 = result_tracker_df[nrow(result_tracker_df), "best_x3"],
                              best_x4 = result_tracker_df[nrow(result_tracker_df), "best_x4"],
                              full_loss = result_tracker_df[nrow(result_tracker_df), "full_loss"],
                              NSE = result_tracker_df[nrow(result_tracker_df), "NSE"],
                              wNSE = result_tracker_df[nrow(result_tracker_df), "wNSE"],
                              x1 = point_tf[[1]],
                              x2 = point_tf[[2]],
                              x3 = point_tf[[3]],
                              n_iteration_used = 1,
                              n_iterations_since_BF_change = result_tracker_df[nrow(result_tracker_df), "n_iterations_since_BF_change"] + 1,
                              stringsAsFactors = FALSE)),
             envir = .GlobalEnv)
    }
  } else {
    assign("result_tracker_df",
           rbind(result_tracker_df,
                 data.frame(best_x1 = result_tracker_df[nrow(result_tracker_df), "best_x1"],
                            best_x3 = result_tracker_df[nrow(result_tracker_df), "best_x3"],
                            best_x4 = result_tracker_df[nrow(result_tracker_df), "best_x4"],
                            full_loss = result_tracker_df[nrow(result_tracker_df), "full_loss"],
                            NSE = result_tracker_df[nrow(result_tracker_df), "NSE"],
                            wNSE = result_tracker_df[nrow(result_tracker_df), "wNSE"],
                            x1 = point_tf[[1]],
                            x2 = point_tf[[2]],
                            x3 = point_tf[[3]],
                            n_iteration_used = 1,
                            n_iterations_since_BF_change = result_tracker_df[nrow(result_tracker_df), "n_iterations_since_BF_change"] + 1,
                            stringsAsFactors = FALSE)),
           envir = .GlobalEnv)
  }
  cat("\nPSO optimization for test nr", Test_number, "results:\n")
  cat("mean NSE:", mean_NSE, "\n")
  if(Test_number %in% c(3.1, 3.2, 3.3)){
    cat("wmean NSE 2km:", weighted.mean(statistics$NSE, w = 1-statistics$NSE), "\n")
    cat("wmean NSE 1km:", weighted.mean(statistics_1km$NSE, w = 1-statistics_1km$NSE), "\n")
  } else {
    cat("wmean NSE:", wmean_NSE, "\n")
  }
  cat("overall loss:", full_loss, "\n")
  result_tracker_df <- get("result_tracker_df", envir = .GlobalEnv)
  cat("\nThe best functions are:\n")
  cat("x1 = ", result_tracker_df[nrow(result_tracker_df), "best_x1"], "\n")
  cat("x2 = 0\n")
  cat("x3 = ", result_tracker_df[nrow(result_tracker_df), "best_x3"], "\n")
  cat("x4 = ", result_tracker_df[nrow(result_tracker_df), "best_x4"], "\n")
  return(full_loss)
}

# 3. PSO Optimization ------------------------------------------------------------------
assign("result_tracker_df", data.frame(best_x1 = "init",
                                       best_x3 = "init",
                                       best_x4 = "init",
                                       full_loss = -9999.0,
                                       NSE = -9999.,
                                       wNSE = -9999.,
                                       x1 = "init",
                                       x2 = "init",
                                       x3 = "init",
                                       n_iteration_used = 0,
                                       n_iterations_since_BF_change = 0,
                                       stringsAsFactors = FALSE), envir = .GlobalEnv)
assign("counter", 0, envir = .GlobalEnv)
iterations <- iterations # now part of the FSO functions
optim <- try(pso::psoptim(rep(NA, 18), fn = {function(x) {
  -objective_function(x, Test_number = Test_number, 
                      true_para_field_df = true_para_field_df,
                      spatial_predictors = spatial_predictors, 
                      parameter_bounds = parameter_bounds,
                      para = para, para_1km = para_1km, 
                      training_basins = training_basins)}},
  lower = rep(-10, 18),
  upper = rep(10, 18),
  control = list(maxit = Inf, maxf = iterations)))
while(class(optim) == "try-error"){
  assign("result_tracker_df", data.frame(best_x1 = "init",
                                         best_x3 = "init",
                                         best_x4 = "init",
                                         full_loss = -9999.0,
                                         NSE = -9999.,
                                         wNSE = -9999.,
                                         x1 = "init",
                                         x2 = "init",
                                         x3 = "init",
                                         n_iteration_used = 0,
                                         n_iterations_since_BF_change = 0,
                                         stringsAsFactors = FALSE), envir = .GlobalEnv)
  assign("counter", 0, envir = .GlobalEnv)
  iterations <- iterations # now part of the FSO function
  optim <- try(pso::psoptim(rep(NA, 18), fn = {function(x) {
    -objective_function(x, Test_number = Test_number, 
                        true_para_field_df = true_para_field_df,
                        spatial_predictors = spatial_predictors, 
                        parameter_bounds = parameter_bounds,
                        para = para, para_1km = para_1km, 
                        training_basins = training_basins)}},
                            lower = rep(-10, 18),
                            upper = rep(10, 18),
                            control = list(maxit = Inf, maxf = iterations)))
}
runs <- sum(result_tracker_df$n_iteration_used)
feather::write_feather(data.frame(result_tracker_df, method = "PSO", n_runs = runs, stringsAsFactors = FALSE),
                       paste0("Test ", substr(Test_number, 1, 1),"/",
                              "Test ", Test_number, "/PSO_GR4J_optimization_", 
                              Test_number, "_run", run, ".feather"))
# Run the model one more time to get the training results for all basins
end_results <- tail(result_tracker_df, 1)
suppressWarnings(evaluate_training_basins(end_results,
                                        run = run, Test_number = Test_number, 
                                        spatial_predictors = spatial_predictors,
                                        Optimizer = Optimizer,
                                        para = para, para_1km = para_1km,
                                        training_basins = training_basins, 
                                        test_basins = test_basins)
)




