# Functions for d-GR4J Function Space optimization case study
# Moritz Feigl, 2019
#

# Function to load spatial predictors
load_sp_mur <- function(scale = TRUE,
                        na.approx = FALSE, only_training_basins = TRUE, 
                        full_dataset = FALSE, training_basins = NULL){
  path <- "Data/spatial_predictors_mur/"
  clay <- raster::raster(paste0(path, "l0_clay_mur.asc"))
  bdim <- raster::raster(paste0(path, "l0_bdim_mur.asc"))
  elevation <- raster::raster(paste0(path, "l0_dem_mur.asc"))
  evi <- raster::raster(paste0(path, "l0_evi_mur.asc"))
  hand <- raster::raster(paste0(path, "l0_hand_mur.asc"))
  noise <- raster::raster(paste0(path, "l0_noise_mur.asc"))
  sand <- raster::raster(paste0(path, "l0_sand_mur.asc"))
  slope <- raster::raster(paste0(path, "l0_slope_mur.asc"))
  nz <- raster::raster(paste0(path, "l0_nz2000_mur.asc"))
  hshade <- raster::raster(paste0(path, "l0_hshade_mur.asc"))
  nb <- raster::raster(paste0(path, "l0_nb2000_mur.asc"))
  n <- length(raster::values(nb))
  # spatial predictor list
  spatial_predictors <- data.frame(nb = raster::values(nb), nz = raster::values(nz), clay = raster::values(clay),
                                   elevation = raster::values(elevation), evi = raster::values(evi),
                                   hand = raster::values(hand), noise = raster::values(noise),
                                   sand = raster::values(sand), slope = raster::values(slope),
                                   bdim = raster::values(bdim))
  # use mask for NZ = subset
  if(!full_dataset){
    if(only_training_basins){
      nb_subset <- training_basins
      spatial_predictors <- spatial_predictors[spatial_predictors$nb %in% nb_subset, ]
    } else {
      spatial_predictors <- spatial_predictors[!is.na(spatial_predictors$nb), ]
    }
  }
  # scale data to [0, 1]
  range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
  
  if(na.approx){
    nas <- which(is.na(spatial_predictors), arr.ind = TRUE)
    for(i in 1:nrow(nas)){
      spatial_predictors[nas[i,1], nas[i, 2]] <- mean(c(
        spatial_predictors[nas[i,1]-1, nas[i,2]],
        spatial_predictors[nas[i,1]+1, nas[i,2]]), na.rm = TRUE)
    }
  }
  if(scale){
    # standardize spatial predictors, except for NB, NZ, noise, evi
    spatial_predictors[, c("clay", "sand")] <- range01(spatial_predictors[, c("clay", "sand")], c(0, 100), na.rm = TRUE)
    spatial_predictors$bdim <- range01(spatial_predictors$bdim, c(0, 2), na.rm = TRUE)
    spatial_predictors$slope <- range01(spatial_predictors$slope, c(0, 90), na.rm = TRUE)
    spatial_predictors$elevation <- range01(spatial_predictors$elevation, c(0, 6000), na.rm = TRUE)
    spatial_predictors$hand <- range01(spatial_predictors$hand, c(0, 6000), na.rm = TRUE)
  }
  spatial_predictors <- tibble::as.tibble(spatial_predictors)
  spatial_predictors[which(spatial_predictors == 0, arr.ind = TRUE)] <- 0.0001
  return(spatial_predictors)
}
# NSE
NSE <- function(observations, predictions){
  1- (sum((predictions - observations)^2, na.rm = TRUE) /
        sum((mean(observations, na.rm = TRUE) - observations)^2, na.rm = TRUE))
}

evaluate_function_from_string <- function(string, l0){
  # Evaluate a function from a string
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
  f_evaluated[is.infinite(f_evaluated)] <- NA
  return(f_evaluated)
}

# Rescale function
rescale <- function(x, to, from = c(-11, 11), ...) { #-11, 11
  if(sd(x, na.rm = TRUE) != 0) {
    return((x - from[1]) / diff(from) * diff(to) + to[1])
  } else {
    return(x)
  }
}

# Generalized mean function: assumes that data is always positive
gmean <- function(x, p, ...){
  if(p == 1){
    gmean <- mean(x, ...)
  } else {
    # check for negative values -> not defined in this version of the generalized mean
    if(sum(x < 0) > 0) stop("Generalized mean with p != 1 is only defined for positive values. Check parameter values!")
    if(abs(p) < 0.001){
      # computational stable formular for geometric mean
      logx <- log(x)
      gmean <- exp(mean(logx[is.finite(logx)]))
    } else {
      # computational stable formular of generalized mean for large p
      x_max <- max(x, ...)
      psum_wmax <- sum((x[x != x_max]/x_max)^p, ...)
      n <- length(x)
      # catch extreme values
      if(is.infinite(psum_wmax)){
        if(p < 0) gmean <- min(x, ...)
        if(p > 0) gmean <- max(x, ...)
      } else{
        gmean <- max(x, ...) * exp( 1/p * (log(1 + psum_wmax) -log(n)))
      }
    }
  }
  return(gmean)
}

# Create GR4J parameter from given tf and spatial predictors
create_GR4J_para <- function(transfer_functions, l0, parameter_bounds, gmean_parameter,
                             km1 = FALSE){
  # create GR4J parameter from a list with transferfunctions
  # Input:
  #    transfer_functions: named list with transfer strings for GR4Jx1, GR4Jx2, GR4Jx3, GR4Jx4
  #    l0: data frame with l0 layer
  #    parameter_bounds: named list with parameter bounds for GR4Jx1, GR4Jx2, GR4Jx3, GR4Jx4
  # Output:
  #    data frame with new GR4J parameters
  new_gr4j_para <- data.frame(NB_ = l0$nb, NZ_ = l0$nz,
                              GR4Jx1 = NA, GR4Jx2 = NA, GR4Jx3 = NA, GR4Jx4 = NA)
  new_gr4j_para$GR4Jx1 <- suppressWarnings(
    evaluate_function_from_string(transfer_functions$GR4Jx1, l0 = l0))
  new_gr4j_para$GR4Jx2 <- suppressWarnings(
    evaluate_function_from_string(transfer_functions$GR4Jx2, l0 = l0))
  new_gr4j_para$GR4Jx3 <- suppressWarnings(
    evaluate_function_from_string(transfer_functions$GR4Jx3, l0 = l0))
  new_gr4j_para$GR4Jx4 <- suppressWarnings(
    evaluate_function_from_string(transfer_functions$GR4Jx4, l0 = l0))
  # scale parameter to parameter bounds
  new_gr4j_para$GR4Jx1 <- round(rescale(new_gr4j_para$GR4Jx1, to = parameter_bounds$GR4Jx1), 2)
  new_gr4j_para$GR4Jx2 <- round(rescale(new_gr4j_para$GR4Jx2, to = parameter_bounds$GR4Jx2), 2)
  new_gr4j_para$GR4Jx3 <- round(rescale(new_gr4j_para$GR4Jx3, to = parameter_bounds$GR4Jx3), 2)
  new_gr4j_para$GR4Jx4 <- round(rescale(new_gr4j_para$GR4Jx4, to = parameter_bounds$GR4Jx4), 2)
  if(km1){
    NB_1km <- raster("Data/km1_NB.asc")
    NZ_1km <- raster("Data/km1_NZ.asc")
    values(NB_1km)[!values(NB_1km) %in% new_gr4j_para$NB_] <- NA
    km1_df <- data.frame(new_gr4j_para[, c("GR4Jx1", "GR4Jx2", "GR4Jx3", "GR4Jx4")],
                         "NB_" = values(NB_1km)[!is.na(values(NB_1km))], 
                         "NZ_" = values(NZ_1km)[!is.na(values(NB_1km))])
    new_gr4j_para <- aggregate(. ~ NZ_ + NB_, km1_df, 
                               function(x) round(gmean(x,  p = gmean_parameter, na.rm = TRUE), 2))
    new_gr4j_para <- new_gr4j_para[order(new_gr4j_para$NZ_), ]
  } else {
    # aggregate to 2km scale
    new_gr4j_para <- aggregate(. ~ NZ_ + NB_, new_gr4j_para, 
                               function(x) round(gmean(x,  p = gmean_parameter, na.rm = TRUE), 2))
    new_gr4j_para <- new_gr4j_para[order(new_gr4j_para$NZ_), ]
  }
  return(new_gr4j_para)
}

# string splitter function
function_splitter <- function(point_tf){
  function_splitted <- unlist(strsplit(point_tf, c("[/^()*+-]")))
  function_splitted <- gsub(" ", "", function_splitted)
  function_splitted <- function_splitted[function_splitted != ""]
  return(function_splitted)
}

# Model size loss
size_loss <- function(functions_splitted){
  length(unlist(functions_splitted)) * 0.001/3
}

# Load true parameter field depending on the Test_number
true_para_field <- function(Test_number){
  path <- "Data/spatial_predictors_mur/"
  library(raster)
  if(Test_number == 2.1){
    # Training and Test basin definitions
    true_para_raster <- raster::raster("True parameters/x1_1km.asc")
    nz <- raster(paste0(path, "l0_nz2000_mur.asc"))
    nb <- raster(paste0(path, "l0_nb2000_mur.asc"))
    true_para_field <- data.frame(nz = values(nz), nb = values(nb), true_para = values(true_para_raster))
    true_para_field <- true_para_field[!is.na(true_para_field$nz), ]
    true_para_field <- true_para_field[!duplicated(true_para_field), ]
  }
  if(Test_number == 2.2){
    # Training and Test basin definitions
    true_para_raster <- raster("True parameters/x3_1km.asc")
    nz <- raster(paste0(path, "l0_nz2000_mur.asc"))
    nb <- raster(paste0(path, "l0_nb2000_mur.asc"))
    true_para_field <- data.frame(nz = values(nz), nb = values(nb), true_para = values(true_para_raster))
    true_para_field <- true_para_field[!is.na(true_para_field$nz), ]
    true_para_field <- true_para_field[!duplicated(true_para_field), ]
  }
  if(Test_number == 2.3){
    # Training and Test basin definitions
    true_para_raster <- raster::raster("True parameters/x4_1km.asc")
    nz <- raster(paste0(path, "l0_nz2000_mur.asc"))
    nb <- raster(paste0(path, "l0_nb2000_mur.asc"))
    true_para_field <- data.frame(nz = values(nz), nb = values(nb), true_para = values(true_para_raster))
    true_para_field <- true_para_field[!is.na(true_para_field$nz), ]
    true_para_field <- true_para_field[!duplicated(true_para_field), ]
  }
  return(true_para_field)
}

# Create a raster object from a transfer function
raster_from_tf <- function(tf, tf_bounds, only_catchment = TRUE, 
                           aggregate = FALSE, gmean_parameter, km1 = FALSE){
  # Input: a transfer function as a string
  # output: a raster object with the parameter field
  path <- "Data/spatial_predictors_mur/"
  # get raster objects of catchment
  sand <- raster::raster(paste0(path, "l0_sand_mur.asc"))
  nb <- raster::raster(paste0(path, "l0_nb2000_mur.asc"))
  nz <- raster::raster(paste0(path, "l0_nz2000_mur.asc"))
  if(aggregate){
    # 0. Load necessary data
    l0_all <- load_sp_mur(na.approx = TRUE, scale = TRUE,
                          only_training_basins = FALSE,
                          full_dataset = FALSE)
    if(km1){
      # 1. Create 250m parameter field
      paraf <- sand
      raster::values(paraf)[is.na(raster::values(nb))] <- NA
      raster::values(paraf)[!is.na(raster::values(nb))] <- rescale(
        evaluate_function_from_string(tf, l0 = l0_all),  to = tf_bounds)
      # 2. get 1km NB and IZ information together with 1km values in a df
      NB_1km <- raster::raster("Data/km1_NB.asc")
      IZ_1km <- raster::raster("Data/km1_IZ.asc")
      NZ_1km <- raster::raster("Data/km1_NZ.asc")
      km1_df <- data.frame("para" = raster::values(paraf)[!is.na(raster::values(nb))],
                           "nb" = raster::values(NB_1km)[!is.na(raster::values(NB_1km))], 
                           "iz" = raster::values(IZ_1km)[!is.na(raster::values(IZ_1km))],
                           "nz" = raster::values(NZ_1km)[!is.na(raster::values(NZ_1km))])
      km1_df$unique_id <- 1:nrow(km1_df)
      # 3. Aggregate to 1km raster
      km1_values <- aggregate(para ~ nz, km1_df, FUN = gmean, 
                              p = gmean_parameter, na.rm = TRUE)
      km1_df_all <- merge(km1_df[, -1], km1_values, by = "nz", 
                          all.x = TRUE)
      km1_df_all <- km1_df_all[order(km1_df_all$unique_id), ]
      
      km1_paraf <- sand
      raster::values(km1_paraf)[is.na(raster::values(nb))] <- NA
      raster::values(km1_paraf)[!is.na(raster::values(nb))] <- km1_df_all$para
      tf_grid <- km1_paraf
    } else {
      # 4. Aggregate to 2km raster
      km2_paraf <- sand
      # define values depending on tf
      raster::values(km2_paraf)[is.na(raster::values(nb))] <- NA
      raster::values(km2_paraf)[!is.na(raster::values(nb))] <- rescale(
        evaluate_function_from_string(tf, l0 = l0_all), 
        to = tf_bounds)
      # create df with parameter values and nb & nz
      grid_values <- data.frame("para" = raster::values(km2_paraf)[!is.na(raster::values(nb))],
                                "nb" = raster::values(nb)[!is.na(raster::values(nb))], 
                                "nz" = raster::values(nz)[!is.na(raster::values(nb))])
      grid_values$unique_id <- 1:nrow(grid_values)
      # aggregate parameter field
      grid_values_agg <- aggregate(para ~ nz, grid_values, FUN = gmean, 
                                   p = gmean_parameter, na.rm = TRUE)
      # add unique id again
      grid_values_all <- merge(grid_values[, -1], grid_values_agg, by = "nz", 
                               all.x = TRUE)
      grid_values_all <- grid_values_all[order(grid_values_all$unique_id), ]
      raster::values(km2_paraf)[!is.na(raster::values(nb))] <- grid_values_all$para
      tf_grid <- km2_paraf
    }
  } else {
    l0_all <- load_sp_mur(na.approx = FALSE, scale = TRUE, 
                          full_dataset = TRUE)
    tf_grid <- sand
    # define values depending on tf
    raster::values(tf_grid) <- rescale(
      evaluate_function_from_string(tf, l0 = l0_all), 
      to = tf_bounds)
    # subset for the catchment area if wanted
    if(only_catchment){
      raster::values(tf_grid)[is.na(values(nb))] <- NA
    }
  }
  return(tf_grid)
}

# DDS Function
dds_fs<- function(xBounds.df, numIter, OBJFUN, search_dim, Test_number, 
                  true_para_field_df, spatial_predictors, parameter_bounds,
                  para, para_1km, training_basins){
  # INPUTS:
  # xBounds.df must be a dataframe with 1st column as minimum, 2nd column as maximum
  # numIter is an integer
  # OBJFUN is a function which returns a scalar value, for which we are trying to minimize.
  #
  # OUTPUTS:
  # outputs.df is a two entry list, containing x_best and y_best, as they evolve over numIter iterations.
  
  # Format xBounds.df colnames
  colnames(xBounds.df) <- c("min", "max")
  # Generate initial first guess
  x_init <- rnorm(search_dim)
  # Evaluate first cost function
  x_evaluated <- OBJFUN(x_init, Test_number = Test_number, 
                        true_para_field_df = true_para_field_df,
                        spatial_predictors = spatial_predictors, 
                        parameter_bounds = parameter_bounds,
                        para = para, para_1km = para_1km, 
                        training_basins = training_basins)
  x_ini <- x_evaluated$`Current point in function space`
  x_best <- matrix(x_init, nrow = 1)
  if(!is.na( x_evaluated$loss)){
    y_init <- x_evaluated$loss
  } else {
    y_init <- -999
  }
  y_best <- y_init
  r = 0.2
  # Select which entry to peturb at each iteration
  peturbIdx <- probPeturb(xBounds.df, numIter)
  # Peturb each entry by N(0,1)*r(x_max - x_min) reflecting if @ boundaries
  sigma <- xBounds.df$max - xBounds.df$min
  for (i in 2:numIter){
    # Set up test x
    x_test <- x_best[i-1, ]
    # Get entries we will peturb
    idx <- peturbIdx[[i]]
    # Initialize vector of peturbations initially zeros with same length of x so we will add this vector to peturb x
    peturbVec <- rep(0, length(x_test))
    # Generate the required number of random normal variables
    N <- rnorm(length(x_test), mean=0, sd=1)
    # Set up vector of peturbations
    peturbVec[idx] <- r*N[idx]*sigma[idx]
    # Temporary resulting x value if we peturbed it
    testPeturb <- x_test + peturbVec
    # Find the values in testPeturb OBJFUN <- wrapper_ofthat have boundary violations.  Store the indices in boundaryViolationsIdx
    boundaryViolationIdx <- which(testPeturb<xBounds.df$min | testPeturb > xBounds.df$max)
    # Reset those violated indices to the opposite peturbation direction
    peturbVec[boundaryViolationIdx]<-(-1*r*N[boundaryViolationIdx]*sigma[boundaryViolationIdx])
    # Find values still at violations of min or max and set them to the minimum or maximum values
    testPeturb<-x_test + peturbVec
    minViolationIdx<-which(testPeturb<xBounds.df$min)
    maxViolationIdx<-which(testPeturb>xBounds.df$max)
    testPeturb[minViolationIdx]<-xBounds.df$min[minViolationIdx]
    testPeturb[maxViolationIdx]<-xBounds.df$max[maxViolationIdx]
    # Peturb the test vector
    x_test <- x_test + peturbVec
    # Evaluate objective function #§ a bit sloppy.. but never mind...
    x_evaluated <- OBJFUN(x_test, Test_number = Test_number, 
                          true_para_field_df = true_para_field_df,
                          spatial_predictors = spatial_predictors, 
                          parameter_bounds = parameter_bounds,
                          para = para, para_1km = para_1km, 
                          training_basins = training_basins)
    x_test <- x_evaluated$`Current point in function space`
    y_test <- x_evaluated$loss
    if(!is.na(y_test)) {
      y_best[i] <- max(c(y_test, y_best[i-1]))
      bestIdx <- which.max(c(y_test, y_best[i-1]))
    } else {
      y_best[i] <- y_best[i-1]
      bestIdx <- 2
    }
    x_choices <- cbind(x_test, x_best[i-1, ])
    x_best <- rbind(x_best, x_choices[,bestIdx])
  }
  output.list <- list(t(x_best), y_best)
  return(output.list)
}


probPeturb <- function(x, numIter){
  # perturber function for DDS
  # Returns numIter length list of entries to be peturbed
  # Input is xBounds & numIter.
  # Returns numIter entry list with the indices which will be peturbed
  xDims <- nrow(x)
  probabilityVector<- 1-log(1:numIter)/log(numIter)
  peturbIdx <- apply(matrix(unlist(lapply(probabilityVector, function(x) as.logical(rbinom(xDims, 1, x)))), byrow=TRUE, ncol=xDims), 1, which)
  return(peturbIdx)
}

SPAEF <- function(observations, predictions){
  # SPAtial EFficiency (SPAEF) metric
  # Inputs: numeric vectors of observawtions and predictions
  # Output: numeric
  spaef_try <- try({
    alpha <- cor(predictions, observations)
    beta <- (sd(observations)/mean(observations)) / (sd(predictions)/mean(predictions))
    
    # scale for histogram distance
    observations <- scale(observations)
    predictions <- scale(predictions)
    
    range_true <- max(observations) - min(observations)
    breaks <- seq(min(observations), max(observations), range_true/100)
    c1 <- as.integer(table(cut(observations, breaks = breaks)))
    c2 <- as.integer(table(cut(predictions, breaks = breaks)))
    
    c_min <- numeric()
    for(i in seq_along(c1)) c_min <- c(c_min, min(c1[i], c2[i]))
    gamma <- sum(c_min)/sum(c1)
    spaef <- 1 - sqrt((alpha-1)^2 + (beta-1)^2 + (gamma-1)^2)
  }, silent = TRUE)
  if(class(spaef_try) == "try-error") spaef <- NA
  return(spaef)
}

# GR4J model quality
GR4J_model_quality <- function(statistics, Test_number, true_para_field_df, 
                               model_size_loss, new_gr4j_para, relevant_basins = NULL,
                               statistics_1km = NULL){
  # Test 1.1 mean NSE
  if(Test_number == 1.1){
    output <- list()
    # subset relevant basins
    if(!is.null(relevant_basins)) statistics <- statistics[statistics$sb %in% relevant_basins, ]    # calculate mean loss for all basins
    mean_NSE <- mean(statistics$NSE)
    wmean_NSE <- weighted.mean(statistics$NSE, w = 1.01-statistics$NSE)
    if(is.infinite(mean_NSE)) mean_NSE <- NA
    if(is.nan(mean_NSE)) mean_NSE <- NA
    # calculate overall loss
    full_loss <- mean_NSE - model_size_loss
    model_loss <- mean_NSE - model_size_loss
    output["mean_NSE"] <- mean_NSE
    output["wmean_NSE"] <- wmean_NSE
    output["model_loss"] <- model_loss
    output["full_loss"] <- full_loss
  }
  
  # test 1.2/1.3 weighted mean NSE
  if(Test_number == 1.2){
    output <- list()
    # subset relevant basins
    if(!is.null(relevant_basins)) statistics <- statistics[statistics$sb %in% relevant_basins, ]
    # calculate mean loss for all basins
    mean_NSE <- mean(statistics$NSE)
    wmean_NSE <- weighted.mean(statistics$NSE, w = 1.01-statistics$NSE)
    if(is.infinite(mean_NSE)) mean_NSE <- NA
    if(is.nan(mean_NSE)) mean_NSE <- NA
    # calculate overall loss
    full_loss <- wmean_NSE - model_size_loss
    model_loss <- wmean_NSE - model_size_loss
    output["mean_NSE"] <- mean_NSE
    output["wmean_NSE"] <- wmean_NSE
    output["model_loss"] <- model_loss
    output["full_loss"] <- full_loss
  }
  
  # Test 2.1-2.3 multi-objective weighted mean NSE using parameter fields
  if(Test_number %in% c(2.1, 2.2, 2.3)){
    output <- list()
    # subset relevant basins
    if(!is.null(relevant_basins)) statistics <- statistics[statistics$sb %in% relevant_basins, ]
    # calculate mean loss for all basins
    mean_NSE <- mean(statistics$NSE)
    wmean_NSE <- weighted.mean(statistics$NSE, w = 1.01-statistics$NSE)
    # Calculate loss by comparison with observed storage parameter field
    if(Test_number == 2.1) estimated_para_field <- new_gr4j_para[, c("NZ_", "NB_", "GR4Jx1")]
    if(Test_number == 2.2) estimated_para_field <- new_gr4j_para[, c("NZ_", "NB_", "GR4Jx3")]
    if(Test_number == 2.3) estimated_para_field <- new_gr4j_para[, c("NZ_", "NB_", "GR4Jx4")]
    
    # in case relevant basins were chosen (testing)
    if(!is.null(relevant_basins)){
      true_para_field_df <- true_para_field_df[true_para_field_df$nb %in% relevant_basins, ]
      estimated_para_field <- estimated_para_field[estimated_para_field$NB_ %in% relevant_basins, ]
    }
    
    merged_paras <- merge(estimated_para_field, true_para_field_df, by.x = c("NZ_", "NB_"), 
                          by.y = c("nz", "nb"))
    para_field_loss <- NSE(observations = merged_paras$true_para, predictions = merged_paras$GR4Jx3)
    # calculate overall loss
    model_loss <- wmean_NSE - model_size_loss
    full_loss <- mean(c(para_field_loss, model_loss))
    output["mean_NSE"] <- mean_NSE
    output["wmean_NSE"] <- wmean_NSE
    output["model_loss"] <- model_loss
    output["full_loss"] <- full_loss
  }
  
  # Test 2.4-2.6 multi-objective weighted mean NSE using states
  if(Test_number %in% c(2.4, 2.5)){
    output <- list()
    # subset relevant basins
    if(!is.null(relevant_basins)) statistics <- statistics[statistics$sb %in% relevant_basins, ]
    
    if(Test_number == 2.4) {
      gr4j_parameter <- "BW0GEB"
      gr4j_state <- "GR4JState1"
      true_states <- feather::read_feather("True parameters/GR4J_state1.feather")
    }
    if(Test_number == 2.5) {
      gr4j_parameter <- "BW3GEB"
      gr4j_state <- "GR4JState2"
      true_states <- feather::read_feather("True parameters/GR4J_state2.feather")
    }
    # Calculate loss by comparison with observed storage states
    # read and format state results
    state_names <- read.table("GR4J_distributed/output/COSERO.plus1", 
                              header=FALSE, fill=TRUE, skip = 1, nrows = 1, 
                              stringsAsFactors = FALSE)
    state_ind <- grep(state_names, pattern = gr4j_parameter)
    state_classes <- rep("NULL", length(state_names))
    state_classes[state_ind] <- "numeric"
    state_classes[1:3] <- "integer"
    states <- read.table("GR4J_distributed/output/COSERO.plus1", 
                         header = FALSE, fill = FALSE, skip = 243,
                         colClasses = state_classes)
    state_classes[c(1:3, state_ind)] <- "character"
    
    names(states) <- read.table("GR4J_distributed/output/COSERO.plus1", 
                                header=FALSE, fill=FALSE, skip = 1,nrows = 1,
                                colClasses = state_classes)
    # get correct names
    state_ind2 <- grep(names(states), pattern = gr4j_parameter)
    names(states)[state_ind2] <- gsub(names(states)[state_ind2], 
                                      pattern = gr4j_parameter, 
                                      replacement = gr4j_state)
    # in case relevant basins were chosen (testing)
    if(!is.null(relevant_basins)){
      relevant_columns <- integer()
      for(basins in seq_along(relevant_basins)){
        current_basin <- paste0("000", relevant_basins[basins])
        current_basin <- substring(current_basin,
                                   first = nchar(current_basin) - 3)
        relevant_columns <- c(relevant_columns, grep(names(states), pattern = current_basin))
      }
      relevant_columns <- relevant_columns[relevant_columns != 0]
      states <- states[, c(1:3, relevant_columns)]
    }
    # get only the modelled catchments
    true_states <- true_states[, names(true_states) %in% names(states)]
    # get only the modelled timesteps
    true_states <- as.data.frame(true_states[
      paste0(true_states$yyyy, "-", true_states$mm, "-", true_states$dd) %in% 
        paste0(states$yyyy, "-" ,states$mm, "-", states$dd), ])
    # create a list for each catchment and calculate state loss
    state_loss <- vector(mode = "list")
    for(i in 4:ncol(states)){
      state_loss[[i-3]] <- data.frame("true" = true_states[, i], "model" = states[, i])
    }
    state_loss_nse <- sapply(state_loss, 
                             function(x) {NSE(observations = x$true, predictions = x$model)})
    # calculate mean total loss for all basins
    mean_NSE <- mean(statistics$NSE)
    wmean_NSE <- weighted.mean(statistics$NSE, w = 1.01-statistics$NSE)
    full_loss <- numeric(length(state_loss_nse))
    for(i in 1:length(state_loss_nse)){
      full_loss[i] <- mean(c(statistics$NSE[i] - model_size_loss, state_loss_nse[i]))
    }
    full_loss <- weighted.mean(full_loss, w = 1.01-full_loss)
    # The mean NSE minus the model size loss -> to be consistent with the other tests
    model_loss <- wmean_NSE - model_size_loss
    # calculate overall loss
    output["mean_NSE"] <- mean_NSE
    output["wmean_NSE"] <- wmean_NSE
    output["model_loss"] <- model_loss
    output["full_loss"] <- full_loss
  }
  if(Test_number == 2.6){
    output <- list()
    # subset relevant basins
    if(!is.null(relevant_basins)) statistics <- statistics[statistics$sb %in% relevant_basins, ]
    state_loss_nse <- list()
    for(states_qualtiy in 1:2){
      if(states_qualtiy == 1) {
        gr4j_parameter <- "BW0GEB"
        gr4j_state <- "GR4JState1"
        true_states <- feather::read_feather("True parameters/GR4J_state1.feather")
      }
      if(states_qualtiy == 2) {
        gr4j_parameter <- "BW3GEB"
        gr4j_state <- "GR4JState2"
        true_states <- feather::read_feather("True parameters/GR4J_state2.feather")
      }
      # Calculate loss by comparison with observed storage states
      # read and format state results
      state_names <- read.table("GR4J_distributed/output/COSERO.plus1", 
                                header=FALSE, fill=TRUE, skip = 1, nrows = 1, 
                                stringsAsFactors = FALSE)
      state_ind <- grep(state_names, pattern = gr4j_parameter)
      state_classes <- rep("NULL", length(state_names))
      state_classes[state_ind] <- "numeric"
      state_classes[1:3] <- "integer"
      states <- read.table("GR4J_distributed/output/COSERO.plus1", 
                           header = FALSE, fill = FALSE, skip = 243,
                           colClasses = state_classes)
      state_classes[c(1:3, state_ind)] <- "character"
      names(states) <- read.table("GR4J_distributed/output/COSERO.plus1", 
                                  header=FALSE, fill=FALSE, skip = 1,nrows = 1,
                                  colClasses = state_classes)
      # get correct names
      state_ind2 <- grep(names(states), pattern = gr4j_parameter)
      names(states)[state_ind2] <- gsub(names(states)[state_ind2], 
                                        pattern = gr4j_parameter, 
                                        replacement = gr4j_state)
      # in case relevant basins were chosen (testing)
      if(!is.null(relevant_basins)){
        relevant_columns <- integer()
        for(basins in seq_along(relevant_basins)){
          current_basin <- paste0("000", relevant_basins[basins])
          current_basin <- substring(current_basin,
                                     first = nchar(current_basin) - 3)
          relevant_columns <- c(relevant_columns, grep(names(states), pattern = current_basin))
        }
        relevant_columns <- relevant_columns[relevant_columns != 0]
        states <- states[, c(1:3, relevant_columns)]
      }
      # get only the modelled catchments
      true_states <- true_states[, names(true_states) %in% names(states)]
      # get only the modelled timesteps
      true_states <- as.data.frame(true_states[
        paste0(true_states$yyyy, "-", true_states$mm, "-", true_states$dd) %in% 
          paste0(states$yyyy, "-" ,states$mm, "-", states$dd), ])
      # create a list for each catchment and calculate state loss
      state_loss <- vector(mode = "list")
      for(i in 4:ncol(states)){
        state_loss[[i-3]] <- data.frame("true" = true_states[, i], "model" = states[, i])
      }
      state_loss_nse[[states_qualtiy]] <- sapply(state_loss, 
                                                 function(x) {
                                                   NSE(observations = x$true, 
                                                       predictions = x$model)})
    }
    # calculate mean total loss for all basins
    mean_NSE <- mean(statistics$NSE)
    wmean_NSE <- weighted.mean(statistics$NSE, w = 1.01-statistics$NSE)
    full_loss <- numeric(length(state_loss_nse[[1]]))
    for(i in 1:length(full_loss)){
      full_loss[i] <- mean(c(statistics$NSE[i] - model_size_loss, 
                             state_loss_nse[[1]][i],
                             state_loss_nse[[2]][i]))
    }
    full_loss <- weighted.mean(full_loss, w = 1.01-full_loss)
    # The mean NSE minus the model size loss -> to be consistent with the other tests
    model_loss <- wmean_NSE - model_size_loss
    # calculate overall loss
    output["mean_NSE"] <- mean_NSE
    output["wmean_NSE"] <- wmean_NSE
    output["model_loss"] <- model_loss
    output["full_loss"] <- full_loss
  }
  
  # Test 3.1-3.3 multi-objective weighted mean NSE using states for both 1km and 2km d-GR4J
  if(Test_number == 3.1){
    output <- list()
    # subset relevant basins
    if(!is.null(relevant_basins)) {
      statistics <- statistics[statistics$sb %in% relevant_basins, ]
      statistics_1km <- statistics_1km[statistics_1km$sb %in% relevant_basins, ]
    }
    # 2 KM NSE
    mean_NSE_2km <- mean(statistics$NSE)
    wmean_NSE_2km <- weighted.mean(statistics$NSE, w = 1.01-statistics$NSE)
    
    # 1 KM NSE
    mean_NSE_1km <- mean(statistics_1km$NSE)
    wmean_NSE_1km <- weighted.mean(statistics_1km$NSE, w = 1.01-statistics_1km$NSE)
    
    # overall NSE
    mean_NSE <- mean(c(mean_NSE_1km, mean_NSE_2km))
    wmean_NSE <- mean(c(wmean_NSE_1km, wmean_NSE_2km))
    if(is.infinite(mean_NSE_2km)) mean_NSE_2km <- NA
    if(is.nan(mean_NSE_2km)) mean_NSE_2km <- NA
    model_loss <- wmean_NSE - model_size_loss
    # calculate overall loss
    full_loss <- mean_NSE - model_size_loss
    output["mean_NSE"] <- mean_NSE_2km
    output["wmean_NSE"] <- wmean_NSE_2km
    output["model_loss"] <- model_loss
    output["full_loss"] <- full_loss
  }
  if(Test_number %in% c(3.2, 3.3)){
    output <- list()
    # subset relevant basins
    if(!is.null(relevant_basins)) {
      statistics <- statistics[statistics$sb %in% relevant_basins, ]
      statistics_1km <- statistics_1km[statistics_1km$sb %in% relevant_basins, ]
    }
    # 2 KM NSE
    mean_NSE_2km <- mean(statistics$NSE)
    wmean_NSE_2km <- weighted.mean(statistics$NSE, w = 1.01-statistics$NSE)
    
    # 1 KM NSE
    mean_NSE_1km <- mean(statistics_1km$NSE)
    wmean_NSE_1km <- weighted.mean(statistics_1km$NSE, w = 1.01-statistics_1km$NSE)
    
    # Losses from states
    if(Test_number == 3.2) {
      gr4j_parameter <- "BW0GEB"
      gr4j_state <- "GR4JState1"
      true_states <- feather::read_feather("True parameters/GR4J_state1.feather")
      true_states_1km <- feather::read_feather("True parameters/GR4J_state1_1km.feather")
    }
    if(Test_number == 3.3) {
      gr4j_parameter <- "BW3GEB"
      gr4j_state <- "GR4JState2"
      true_states <- feather::read_feather("True parameters/GR4J_state2.feather")
      true_states_1km <- feather::read_feather("True parameters/GR4J_state2_1km.feather")
    }
    # Calculate state loss
    state_loss_nse <- list("1km" = NA, "2km" = NA)
    scale_names <- list("1km" = "_1km", "2km" = "")
    for(scale in c("1km", "2km")){
      scale_name <- scale_names[[scale]]
      # Calculate loss by comparison with observed storage states
      # read and format state results
      state_names <- read.table(paste0("GR4J_distributed", scale_name, "/output/COSERO.plus1"), 
                                header = FALSE, fill = TRUE, skip = 1, nrows = 1, 
                                stringsAsFactors = FALSE)
      state_ind <- grep(state_names, pattern = gr4j_parameter)
      state_classes <- rep("NULL", length(state_names))
      state_classes[state_ind] <- "numeric"
      state_classes[1:3] <- "integer"
      states <- read.table(paste0("GR4J_distributed", scale_name, "/output/COSERO.plus1"), 
                           header = FALSE, fill = FALSE, skip = 243,
                           colClasses = state_classes)
      state_classes[c(1:3, state_ind)] <- "character"
      
      names(states) <- read.table(paste0("GR4J_distributed", scale_name, "/output/COSERO.plus1"), 
                                  header = FALSE, fill = FALSE, skip = 1,nrows = 1,
                                  colClasses = state_classes)
      # get correct names
      state_ind2 <- grep(names(states), pattern = gr4j_parameter)
      names(states)[state_ind2] <- gsub(names(states)[state_ind2], 
                                        pattern = gr4j_parameter, 
                                        replacement = gr4j_state)
      # in case relevant basins were chosen
      if(!is.null(relevant_basins)){
        relevant_columns <- integer()
        for(basins in seq_along(relevant_basins)){
          current_basin <- paste0("000", relevant_basins[basins])
          current_basin <- substring(current_basin,
                                     first = nchar(current_basin) - 3)
          relevant_columns <- c(relevant_columns, grep(names(states), pattern = current_basin))
        }
        relevant_columns <- relevant_columns[relevant_columns != 0]
        states <- states[, c(1:3, relevant_columns)]
      }
      true_states_scale <- get(paste0("true_states", scale_name))
      # get only the modelled catchments
      true_states_scale <- true_states_scale[, names(true_states) %in% names(states)]
      # get only the modelled timesteps
      true_states_scale <- as.data.frame(true_states_scale[
        paste0(true_states_scale$yyyy, "-", true_states_scale$mm, "-", true_states_scale$dd) %in% 
          paste0(states$yyyy, "-" ,states$mm, "-", states$dd), ])
      # create a list for each catchment and calculate state loss
      state_loss <- vector(mode = "list")
      for(i in 4:ncol(states)){
        state_loss[[i-3]] <- data.frame("true" = true_states_scale[, i], "model" = states[, i])
      }
      state_loss_nse[[scale]] <- sapply(state_loss, 
                                        function(x) {NSE(observations = x$true, predictions = x$model)})
    }
    # calculate overall mean NSE
    mean_NSE <- mean(c(mean_NSE_1km, mean_NSE_2km))
    wmean_NSE <- mean(c(wmean_NSE_1km, wmean_NSE_2km))
    if(is.infinite(mean_NSE_2km)) mean_NSE_2km <- NA
    if(is.nan(mean_NSE_2km)) mean_NSE_2km <- NA
    # calculate full loss 
    full_loss <- numeric(length(state_loss_nse))
    for(i in 1:length(state_loss_nse)){
      full_loss_i_2km <- mean(c(statistics$NSE[i] - model_size_loss, state_loss_nse[["2km"]][i]))
      full_loss_i_1km <- mean(c(statistics_1km$NSE[i] - model_size_loss, state_loss_nse[["1km"]][i]))
      full_loss[i] <- mean(full_loss_i_2km, full_loss_i_1km)
    }
    full_loss <- weighted.mean(full_loss, w = 1.01-full_loss)
    model_loss <- wmean_NSE - model_size_loss
    output["mean_NSE"] <- mean_NSE_2km
    output["wmean_NSE"] <- wmean_NSE_2km
    output["model_loss"] <- model_loss
    output["full_loss"] <- full_loss
  }
  
  # Test 4.1-4.6 multi-objective optimization with NSE for Q and SPAEF for states
  # Test 4.1-4.3 using SPAEF only for the last maps of states
  if(Test_number %in% c(4.1, 4.2)){
    output <- list()
    # subset relevant basins
    if(!is.null(relevant_basins)) statistics <- statistics[statistics$sb %in% relevant_basins, ]
    
    if(Test_number == 4.1) {
      gr4j_state <- "GR4JSt1"
      true_last_states <- as.data.frame(
        feather::read_feather(
          paste0("True parameters/", training, "_last_GR4J_state1.feather")))
    }
    if(Test_number == 4.2) {
      gr4j_state <- "GR4JSt2"
      true_last_states <- as.data.frame(
        feather::read_feather(
          paste0("True parameters/", training, "_last_GR4J_state2.feather")))
    }
    
    #SPAEF
    spaef_try <- try({
      # Calculate loss by comparison with observed storage states
      last_state_names <- read.table("GR4J_distributed/output/statevar_GR4J.dmp", 
                                     header=FALSE, fill=TRUE, skip = 2, nrows = 1, 
                                     stringsAsFactors = FALSE)
      
      last_states <- read.table("GR4J_distributed/output/statevar_GR4J.dmp", 
                                header = FALSE, fill = FALSE, skip = 3, nrows = 2859)
      colnames(last_states) <- last_state_names
      last_states <-last_states[last_states$NB %in% relevant_basins, 
                                c("NB", "IZ", gr4j_state)]
      true_last_states <- true_last_states[true_last_states$NB %in% relevant_basins, ]
      
      a <- true_last_states[, gr4j_state]
      a <- a[is.finite(a)]
      b <- last_states[, gr4j_state]
      b <- b[is.finite(b)]
      alpha <- cor(b, a)
      beta <- (sd(a)/mean(a)) / (sd(b)/mean(b))
      
      range_true <- max(a) - min(a)
      breaks <- seq(min(a), max(a), range_true/100)
      c1 <- as.integer(table(cut(a, breaks = breaks)))
      c2 <- as.integer(table(cut(b, breaks = breaks)))
      
      c_min <- numeric()
      for(i in seq_along(c1)) c_min <- c(c_min, min(c1[i], c2[i]))
      gamma <- sum(c_min)/sum(c1)
      spaef <- 1 - sqrt((alpha-1)^2 + (beta-1)^2 + (gamma-1)^2)
    })
    if(class(spaef_try) == "try-error") spaef <- NA
    
    # calculate mean total loss for all basins
    mean_NSE <- mean(statistics$NSE)
    wmean_NSE <- weighted.mean(statistics$NSE, w = 1.01-statistics$NSE)
    full_loss <- wmean_NSE + spaef - model_size_loss
    
    # The mean NSE minus the model size loss -> to be consistent with the other tests
    model_loss <- wmean_NSE - model_size_loss
    
    # calculate overall loss
    output["mean_NSE"] <- mean_NSE
    output["wmean_NSE"] <- wmean_NSE
    output["model_loss"] <- spaef
    output["full_loss"] <- full_loss
  }
  if(Test_number == 4.3){
    output <- list()
    # subset relevant basins
    if(!is.null(relevant_basins)) statistics <- statistics[statistics$sb %in% relevant_basins, ]
    
    # GR4J States 1
    gr4j_state <- "GR4JSt1"
    true_last_states <- as.data.frame(
      feather::read_feather(
        paste0("True parameters/", training, "_last_GR4J_state1.feather")))
    
    
    
    
    #SPAEF
    spaef_try1 <- try({
      # Calculate loss by comparison with observed storage states
      last_state_names <- read.table("GR4J_distributed/output/statevar_GR4J.dmp", 
                                     header=FALSE, fill=TRUE, skip = 2, nrows = 1, 
                                     stringsAsFactors = FALSE)
      
      last_states <- read.table("GR4J_distributed/output/statevar_GR4J.dmp", 
                                header = FALSE, fill = FALSE, skip = 3, nrows = 2859)
      colnames(last_states) <- last_state_names
      last_states <-last_states[last_states$NB %in% relevant_basins, 
                                c("NB", "IZ", gr4j_state)]
      true_last_states <- true_last_states[true_last_states$NB %in% relevant_basins, ]
      
      a <- true_last_states[, gr4j_state]
      a <- a[is.finite(a)]
      b <- last_states[, gr4j_state]
      b <- b[is.finite(b)]
      alpha <- cor(b, a)
      beta <- (sd(a)/mean(a)) / (sd(b)/mean(b))
      
      range_true <- max(a) - min(a)
      breaks <- seq(min(a), max(a), range_true/100)
      c1 <- as.integer(table(cut(a, breaks = breaks)))
      c2 <- as.integer(table(cut(b, breaks = breaks)))
      
      c_min <- numeric()
      for(i in seq_along(c1)) c_min <- c(c_min, min(c1[i], c2[i]))
      gamma <- sum(c_min)/sum(c1)
      spaef1 <- 1 - sqrt((alpha-1)^2 + (beta-1)^2 + (gamma-1)^2)
    })
    if(class(spaef_try1) == "try-error") spaef1 <- NA
    
    
    
    # GR4J States 2
    gr4j_state <- "GR4JSt2"
    true_last_states <- as.data.frame(
      feather::read_feather(
        paste0("True parameters/", training, "_last_GR4J_state2.feather")))
    
    # Calculate loss by comparison with observed storage states
    last_state_names <- read.table("GR4J_distributed/output/statevar_GR4J.dmp", 
                                   header=FALSE, fill=TRUE, skip = 2, nrows = 1, 
                                   stringsAsFactors = FALSE)
    
    last_states <- read.table("GR4J_distributed/output/statevar_GR4J.dmp", 
                              header = FALSE, fill = FALSE, skip = 3, nrows = 2859)
    colnames(last_states) <- last_state_names
    last_states <-last_states[last_states$NB %in% relevant_basins, 
                              c("NB", "IZ", gr4j_state)]
    true_last_states <- true_last_states[true_last_states$NB %in% relevant_basins, ]
    #SPAEF
    spaef_try2 <- try({
      # Calculate loss by comparison with observed storage states
      last_state_names <- read.table("GR4J_distributed/output/statevar_GR4J.dmp", 
                                     header=FALSE, fill=TRUE, skip = 2, nrows = 1, 
                                     stringsAsFactors = FALSE)
      
      last_states <- read.table("GR4J_distributed/output/statevar_GR4J.dmp", 
                                header = FALSE, fill = FALSE, skip = 3, nrows = 2859)
      colnames(last_states) <- last_state_names
      last_states <-last_states[last_states$NB %in% relevant_basins, 
                                c("NB", "IZ", gr4j_state)]
      true_last_states <- true_last_states[true_last_states$NB %in% relevant_basins, ]
      
      a <- true_last_states[, gr4j_state]
      a <- a[is.finite(a)]
      b <- last_states[, gr4j_state]
      b <- b[is.finite(b)]
      alpha <- cor(b, a)
      beta <- (sd(a)/mean(a)) / (sd(b)/mean(b))
      
      range_true <- max(a) - min(a)
      breaks <- seq(min(a), max(a), range_true/100)
      c1 <- as.integer(table(cut(a, breaks = breaks)))
      c2 <- as.integer(table(cut(b, breaks = breaks)))
      
      c_min <- numeric()
      for(i in seq_along(c1)) c_min <- c(c_min, min(c1[i], c2[i]))
      gamma <- sum(c_min)/sum(c1)
      spaef2 <- 1 - sqrt((alpha-1)^2 + (beta-1)^2 + (gamma-1)^2)
    })
    if(class(spaef_try2) == "try-error") spaef2 <- NA
    
    if(sum(is.na(c(spaef1, spaef2))) == 0){
      spaef <- mean(spaef1, spaef2)
    } else spaef <- NA
    
    # calculate mean total loss for all basins
    mean_NSE <- mean(statistics$NSE)
    wmean_NSE <- weighted.mean(statistics$NSE, w = 1.01-statistics$NSE)
    full_loss <- wmean_NSE + spaef - model_size_loss
    
    # The mean NSE minus the model size loss -> to be consistent with the other tests
    model_loss <- wmean_NSE - model_size_loss
    
    # calculate overall loss
    output["mean_NSE"] <- mean_NSE
    output["wmean_NSE"] <- wmean_NSE
    output["model_loss"] <- spaef
    output["full_loss"] <- full_loss
  }
  # Test 4.4-4.6 using SPAEF for the time series of state maps
  if(Test_number %in% c(4.4, 4.5)){
    output <- list()
    # subset relevant basins
    if(!is.null(relevant_basins)) statistics <- statistics[statistics$sb %in% relevant_basins, ]
    
    gr4j_state <- ifelse(Test_number == 4.4, "GR4JSt1", "GR4JSt2")
    # Calculate loss by comparison with observed storage states
    # read and format state results
    state_files <- list.files("GR4J_distributed/cdr/output", full.names = TRUE)
    state_list <- lapply(state_files, function(x) read.table(x, header=TRUE, sep = ","))
    
    if(!is.null(relevant_basins)){
      if(sum(relevant_basins %in% test_basins) == length(relevant_basins)){
        training_nz <- feather::read_feather("True parameters/training_basins_nz.feather")
        state_list <- lapply(state_list, function(x){
          x[!(x$NZ %in% training_nz$NZ), ]
        })
      }
      if(sum(relevant_basins %in% training_basins) == length(relevant_basins)){
        training_nz <- feather::read_feather("True parameters/training_basins_nz.feather")
        state_list <- lapply(state_list, function(x){
          x[x$NZ %in% training_nz$NZ, ]
        })
      }
      
    }
    
    # standardize states
    
    # get relevant state and standardize
    state_list <- lapply(state_list, function(x) {
      data.frame(NZ = x$NZ, gr4j_state = x[, gr4j_state])})
    all_state_list <- Map(merge, 
                          state_list, 
                          true_state_list[1:length(state_list)], 
                          by="NZ", all.y = FALSE)
    
    # For testing use only the testing time steps
    if(training == "testing") all_state_list <- all_state_list[2434:length(all_state_list)]
    if(training == "training") all_state_list <- all_state_list[241:length(all_state_list)]
    
    # compute spaef for every timestep
    all_spaef <- sapply(all_state_list, function(x) {
      SPAEF(observations = x[, 2], prediction = x[, 3])
    })
    spaef <- mean(all_spaef)
    
    # calculate loss
    mean_NSE <- mean(statistics$NSE)
    wmean_NSE <- weighted.mean(statistics$NSE, w = 1.01-statistics$NSE)
    if(is.na(wmean_NSE) | is.na(spaef)){
      full_loss <- NA
    } else full_loss <- mean(c(wmean_NSE, spaef)) - model_size_loss
    
    # output
    output["mean_NSE"] <- mean_NSE
    output["wmean_NSE"] <- wmean_NSE
    output["model_loss"] <- spaef
    output["full_loss"] <- full_loss
  }
  if(Test_number == 4.6){
    output <- list()
    # subset relevant basins
    if(!is.null(relevant_basins)) statistics <- statistics[statistics$sb %in% relevant_basins, ]
    
    # Calculate loss by comparison with observed storage states
    # read and format state results
    state_files <- list.files("GR4J_distributed/cdr/output", full.names = TRUE)
    state_list <- lapply(state_files, function(x) read.table(x, header=TRUE, sep = ","))
    
    if(!is.null(relevant_basins)){
      if(sum(relevant_basins %in% test_basins) == length(relevant_basins)){
        training_nz <- feather::read_feather("True parameters/training_basins_nz.feather")
        state_list <- lapply(state_list, function(x){
          x[!(x$NZ %in% training_nz$NZ), ]
        })
      }
      if(sum(relevant_basins %in% training_basins) == length(relevant_basins)){
        training_nz <- feather::read_feather("True parameters/training_basins_nz.feather")
        state_list <- lapply(state_list, function(x){
          x[x$NZ %in% training_nz$NZ, ]
        })
      }
    }
    
    spaef_list <- list("GR4JSt1" = NA,
                       "GR4JSt2" = NA)
    for(gr4j_state in c("GR4JSt1", "GR4JSt2")){
      
      # get relevant state and standardize
      state_list_sub <- lapply(state_list, function(x) {
        data.frame(NZ = x$NZ, gr4j_state = x[, gr4j_state])})
      true_state_list_sub <- lapply(true_state_list[1:length(state_list)], 
                                    function(x) {
                                      data.frame(NZ = x$NZ, gr4j_state = x[, gr4j_state])})
      all_state_list <- Map(merge, 
                            state_list_sub, 
                            true_state_list_sub, 
                            by="NZ", all.y = FALSE)
      
      # For testing use only the testing time steps
      if(training == "testing") all_state_list <- all_state_list[2434:length(all_state_list)]
      if(training == "training") all_state_list <- all_state_list[241:length(all_state_list)]
      
      # compute spaef for every timestep
      all_spaef <- sapply(all_state_list, function(x) {
        SPAEF(observations = x[, 2], prediction = x[, 3])
      })
      spaef_list[[gr4j_state]] <- mean(all_spaef)
    }
    spaef <- mean(unlist(spaef_list))
    # calculate loss
    mean_NSE <- mean(statistics$NSE)
    wmean_NSE <- weighted.mean(statistics$NSE, w = 1.01-statistics$NSE)
    if(is.na(wmean_NSE) | is.na(spaef)){
      full_loss <- NA
    } else full_loss <- mean(c(wmean_NSE, spaef)) - model_size_loss
    
    # output
    output["mean_NSE"] <- mean_NSE
    output["wmean_NSE"] <- wmean_NSE
    output["model_loss"] <- spaef
    output["full_loss"] <- full_loss
  }
  
  return(output)
}

# Evaluate test basin and test time period performance
evaluate_test_basins <- function(test_functions, Optimizer, Test_number, run,
                                 para, para_1km, training_basins, test_basins,
                                 true_state_list){
  true_para_field_df <- true_para_field(Test_number)
  if(Test_number %in% c(4.4, 4.5)){
    if(!exists("true_state_list")) {
      true_state_list <- readRDS("True parameters/true_states_list")
      gr4j_state <- ifelse(Test_number == 4.4, "GR4JSt1", "GR4JSt2")
      true_state_list <- lapply(true_state_list, function(x) {
        data.frame(NZ = x$NZ, gr4j_state = x[, gr4j_state])})
    }
  }
  if(Test_number == 4.6){
    if(!exists("true_state_list")) {
      true_state_list <- readRDS("True parameters/true_states_list")
    }
  }
  
  l0 <- load_sp_mur(na.approx = TRUE, scale = TRUE,
                    only_training_basins = FALSE,
                    full_dataset = FALSE)
  test_functions$GR4Jx2 <- "0"
  new_gr4j_para <- create_GR4J_para(transfer_functions = test_functions,
                                    l0 = l0,
                                    parameter_bounds = parameter_bounds,
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
  cat("\n", file = "GR4J_distributed/input/para_Mur_GR4J_fsOptim.txt")
  suppressWarnings(write.table(para_new, "GR4J_distributed/input/para_Mur_GR4J_fsOptim.txt",
                               append = TRUE, row.names = FALSE, quote = FALSE))
  
  setwd("GR4J_distributed/")
  sys::exec_wait("start_GR4J_case_study_full_run.bat",
                 std_out = "GR4J_output.txt")
  setwd("..")
  # 2 KM statistics
  statistics <- read.table("GR4J_distributed/output/statistics_gr4j_Mur.txt",
                           skip = 21, header = TRUE)
  if(Test_number %in% c(3.1, 3.2, 3.3)){
    # 1 KM Model parameters
    new_gr4j_para_1km <- create_GR4J_para(transfer_functions = test_functions, 
                                          l0 = l0, 
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
    sys::exec_wait("start_GR4J_case_study_full_run.bat",
                   std_out = "GR4J_output.txt")
    setwd("..")
    # Get statistics and calculate losses
    statistics_1km <- read.table("GR4J_distributed_1km/output/statistics_gr4j_Mur.txt",
                                 skip = 21, header = TRUE)
  } else {
    statistics_1km <- NULL
  }
  
  # calculate mean loss for all basins
  mean_NSE <- mean(statistics$NSE)
  mean_train_NSE <-  mean(statistics$NSE[statistics$sb %in% training_basins])
  mean_test_NSE <-  mean(statistics$NSE[statistics$sb %in% test_basins])
  # Model size loss
  functions_splitted <- lapply(test_functions, function_splitter)
  model_size_loss <- size_loss(functions_splitted)
  # All losses for all basins
  evaluation <- GR4J_model_quality(statistics = statistics, Test_number = Test_number, 
                                   true_para_field_df = true_para_field_df, 
                                   model_size_loss = model_size_loss,
                                   new_gr4j_para = new_gr4j_para, 
                                   statistics_1km = statistics_1km,
                                   relevant_basins = c(training_basins, test_basins),
                                   training = FALSE, true_state_list = true_state_list)
  # test losses
  test_evaluation <- GR4J_model_quality(statistics = statistics, Test_number = Test_number, 
                                        true_para_field_df = true_para_field_df, 
                                        model_size_loss = model_size_loss,
                                        new_gr4j_para = new_gr4j_para,
                                        relevant_basins = test_basins,
                                        statistics_1km = statistics_1km,
                                        training = FALSE, true_state_list = true_state_list)
  # train losses
  train_evaluation <- GR4J_model_quality(statistics = statistics, Test_number = Test_number, 
                                         true_para_field_df = true_para_field_df, 
                                         model_size_loss = model_size_loss,
                                         new_gr4j_para = new_gr4j_para,
                                         relevant_basins = training_basins,
                                         statistics_1km = statistics_1km,
                                         training = FALSE, true_state_list = true_state_list)
  
  names(statistics)[1] <- "Basin"
  results <- list(data.frame("mean NSE" = round(mean_NSE, 3),
                             "mean training NSE" = round(mean_train_NSE, 3),
                             "mean test NSE" = round(mean_test_NSE, 3), check.names = FALSE),
                  # overall model results
                  data.frame(
                    "mean NSE" = round(evaluation[["mean_NSE"]], 3),
                    "weighted mean NSE" = round(evaluation[["wmean_NSE"]], 3),
                    "SPAEF or model loss" = round(evaluation[["model_loss"]], 3), 
                    "full loss" = round(evaluation[["full_loss"]], 3), check.names = FALSE),
                  # training model results
                  data.frame(
                    "mean NSE" = round(train_evaluation[["mean_NSE"]], 3),
                    "weighted mean NSE" = round(train_evaluation[["wmean_NSE"]], 3),
                    "SPAEF or model loss" = round(train_evaluation[["model_loss"]], 3), 
                    "full loss" = round(train_evaluation[["full_loss"]], 3), check.names = FALSE),
                  
                  # test model results
                  data.frame(
                    "mean NSE" = round(test_evaluation[["mean_NSE"]], 3),
                    "weighted mean NSE" = round(test_evaluation[["wmean_NSE"]], 3),
                    "SPAEF or model loss" = round(test_evaluation[["model_loss"]], 3), 
                    "full loss" = round(test_evaluation[["full_loss"]], 3), check.names = FALSE),
                  statistics[, 1:2])
  
  cat("\nSaved testing results of test nr.", Test_number, "in corresponding folder\n")
  file <- paste0("Test ", substr(Test_number, 1, 1),"/",
                 "Test ", Test_number, "/testing/",
                 Optimizer, "_testing_", Test_number,
                 "_run", run, ".txt")
  cat("Testing results for test number", Test_number, ":", Optimizer, 
      "- run", run, "\n\n", file = file)
  cat("General results:\n", file = file, append = TRUE)
  suppressWarnings(write.table(t(results[[1]]), file = file, append = TRUE, 
                               col.names = FALSE, quote = FALSE))
  cat("\n", file = file, append = TRUE)
  # all basins
  cat("All basins:\n", file = file, append = TRUE)
  suppressWarnings(write.table(t(results[[2]]), file = file, append = TRUE, 
                               col.names = FALSE, quote = FALSE))
  cat("\n", file = file, append = TRUE)
  # training basins
  cat("Training basins:\n", file = file, append = TRUE)
  suppressWarnings(write.table(t(results[[3]]), file = file, append = TRUE, 
                               col.names = FALSE, quote = FALSE))
  cat("\n", file = file, append = TRUE)
  # Test basins
  cat("Test basins:\n", file = file, append = TRUE)
  suppressWarnings(write.table(t(results[[4]]), file = file, append = TRUE, 
                               col.names = FALSE, quote = FALSE))
  cat("\n", file = file, append = TRUE)
  
  cat("All basins NSE:\n", file = file, append = TRUE)
  suppressWarnings(write.table(results[[5]], file = file, append = TRUE, 
                               row.names = FALSE, quote = FALSE))
  cat("\n\nThe tested functions are:\n", file = file, append = TRUE)
  cat("x1 = ", test_functions[["GR4Jx1"]], "\n", file = file, append = TRUE)
  cat("x2 = 0\n", file = file, append = TRUE)
  cat("x3 = ", test_functions[["GR4Jx3"]], "\n", file = file, append = TRUE)
  cat("x4 = ", test_functions[["GR4Jx4"]], "\n", file = file, append = TRUE)
  return(results)
}
# Evaluate training basins and training time period performance
evaluate_training_basins <- function(end_results, run, Test_number, 
                                     spatial_predictors, Optimizer, para, para_1km = NULL,
                                     training_basins, test_basins){
  best_tfs <- list(GR4Jx1 = end_results$best_x1,
                   GR4Jx2 = "0",          
                   GR4Jx3 = end_results$best_x3,
                   GR4Jx4 = end_results$best_x4)
  new_gr4j_para <- try(create_GR4J_para(transfer_functions = best_tfs,
                                        l0 = spatial_predictors,
                                        parameter_bounds = parameter_bounds,
                                        gmean_parameter = gmean_parameter),
                       silent = TRUE)
  if(class(new_gr4j_para) == "try-error") {
    stop("Failed evaluating training results. No valid transfer function was found.")
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
  
  if(Test_number %in% c(3.1, 3.2, 3.3)){
    # 1 KM Model parameters
    new_gr4j_para_1km <- create_GR4J_para(transfer_functions = best_tfs, 
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
  
  # Get statistics and calculate losses
  statistics <- read.table("GR4J_distributed/output/statistics_gr4j_Mur.txt",
                           skip = 21, header = TRUE)
  
  # model size loss
  functions_splitted <- lapply(best_tfs, function_splitter)
  model_size_loss <- size_loss(functions_splitted)
  true_para_field_df <- true_para_field(Test_number)
  # Evaluate model quality
  if(Test_number %in% c(4.4, 4.5)){
    if(!exists("true_state_list")) {
      true_state_list <- readRDS("True parameters/true_states_list")
      gr4j_state <- ifelse(Test_number == 4.4, "GR4JSt1", "GR4JSt2")
      true_state_list <- lapply(true_state_list, function(x) {
        data.frame(NZ = x$NZ, gr4j_state = x[, gr4j_state])})
    }
  }
  if(Test_number == 4.6){
    if(!exists("true_state_list")) {
      true_state_list <- readRDS("True parameters/true_states_list")
    }
  }  
  evaluation <- GR4J_model_quality(statistics = statistics, Test_number = Test_number, 
                                   true_para_field_df = true_para_field_df, 
                                   model_size_loss = model_size_loss,
                                   new_gr4j_para = new_gr4j_para, 
                                   statistics_1km = statistics_1km,
                                   relevant_basins = training_basins,
                                   true_state_list = true_state_list)
  # define training result df
  names(statistics)[1] <- "Basin"
  train_results <- list(data.frame("mean_NSE" = round(evaluation[["mean_NSE"]], 3), 
                                   "weighted_mean_NSE" = round(evaluation[["wmean_NSE"]], 3),
                                   "SPAEF/model_loss" = round(evaluation[["model_loss"]], 3), 
                                   "full_loss" = round(evaluation[["full_loss"]], 3)),
                        statistics[, 1:2])
  # save and cat
  cat("\nSaved best training results of test nr.", Test_number, "in corresponding folder\n")
  file <- paste0("Test ", substr(Test_number, 1, 1),"/",
                 "Test ", Test_number, "/training/",
                 Optimizer, "_training_", Test_number,
                 "_run", run, ".txt")
  cat("Training results for test number", Test_number, ":", Optimizer, 
      "- run", run, "\n\n", file = file)
  suppressWarnings(write.table(t(train_results[[1]]), file = file, append = TRUE, 
                               col.names = FALSE, quote = FALSE))
  cat("\n", file = file, append = TRUE)
  suppressWarnings(write.table(train_results[[2]], file = file, 
                               append = TRUE, row.names = FALSE, quote = FALSE))
  cat("\n\nThe optimized functions are:\n", file = file, append = TRUE)
  cat("x1 = ", best_tfs[["GR4Jx1"]], "\n", file = file, append = TRUE)
  cat("x2 = 0\n", file = file, append = TRUE)
  cat("x3 = ", best_tfs[["GR4Jx3"]], "\n", file = file, append = TRUE)
  cat("x4 = ", best_tfs[["GR4Jx4"]], "\n", file = file, append = TRUE)
}

# Function Space Optimization
FSO <- function(Optimizer, Test_number, run, iterations,
                training_basins, test_basins){
  # 1. Setup
  cat("\n***", "Test number", Test_number, "-", Optimizer, "optimization run", run,  "***\n")
  if(Test_number == 1.1) cat("Test info: optimization of mean NSE\n")
  if(Test_number == 1.2) cat("Test info: optimization of weighted mean NSE\n")
  if(Test_number == 2.1) cat("Test info: optimization of weighted mean NSE and X1 parameter field\n")
  if(Test_number == 2.2) cat("Test info: optimization of weighted mean NSE and X3 parameter field\n")
  if(Test_number == 2.3) cat("Test info: optimization of weighted mean NSE and X4 parameter field\n")
  if(Test_number == 2.4) cat("Test info: optimization of weighted mean NSE and GR4J state S NSE\n")
  if(Test_number == 2.5) cat("Test info: optimization of weighted mean NSE and GR4J state R NSE\n")
  if(Test_number == 2.5) cat("Test info: optimization of weighted mean NSE and GR4J states S & R NSE\n")
  if(Test_number == 3.1) cat("Test info: optimization of weighted mean NSE of 2 km and 1 km model")
  if(Test_number == 3.2) cat("Test info: optimization of weighted mean NSE of 2 km and 1 km model 
                             and GR4J state S NSE\n")
  if(Test_number == 3.3) cat("Test info: optimization of weighted mean NSE of 2 km and 1 km model 
                             and GR4J state R NSE\n")
  if(Test_number == 4.1) cat("Test info: optimization of weighted mean NSE and the last time step GR4J state S SPAEF\n")
  if(Test_number == 4.2) cat("Test info: optimization of weighted mean NSE and the last time step GR4J state R SPAEF\n")
  if(Test_number == 4.3) cat("Test info: optimization of weighted mean NSE and the last time step GR4J states S & R SPAEF\n")
  if(Test_number == 4.4) cat("Test info: optimization of weighted mean NSE and GR4J state S SPAEF\n")
  if(Test_number == 4.4) cat("Test info: optimization of weighted mean NSE and GR4J state R SPAEF\n")
  if(Test_number == 4.4) cat("Test info: optimization of weighted mean NSE and GR4J states S & R SPAEF\n")
  
  # Generate Test specific folders in directory
  general_test_folder <- paste0("Test ", substr(Test_number, 1, 1))
  subtest_folder <- paste0("Test ", substr(Test_number, 1, 1),"/",
                           "Test ", Test_number)
  training_folder <- paste0("Test ", substr(Test_number, 1, 1),"/",
                            "Test ", Test_number, "/training")
  testing_folder <- paste0("Test ", substr(Test_number, 1, 1),"/",
                           "Test ", Test_number, "/testing")
  
  if (!dir.exists(general_test_folder)){
    dir.create(general_test_folder)
  }
  if (!dir.exists(subtest_folder)){
    dir.create(subtest_folder)
  }
  if (!dir.exists(training_folder)){
    dir.create(training_folder)
  }
  if (!dir.exists(testing_folder)){
    dir.create(testing_folder)
  }
  # Plot paths
  # path for saving rasters
  para_fields <- paste0("Test ", substr(Test_number, 1, 1),"/",
                        "Test ", Test_number, "/parameter fields/")
  if (!dir.exists(para_fields)){
    dir.create(para_fields)
  }
  # paths for run specific para fields
  para_fields2 <- paste0(para_fields, "run_", run, "/")
  if (!dir.exists(para_fields2)){
    dir.create(para_fields2)
  }
  if (!dir.exists(paste0(para_fields2, "Plots/"))){
    dir.create(paste0(para_fields2, "Plots/"))
  }
  # Path for plots
  diag_path1 <- paste0(para_fields, "../diagnostic plots/")
  diag_path2 <- paste0(para_fields, "../diagnostic plots/run_", run, "/")
  if (!dir.exists(diag_path1)){
    dir.create(diag_path1)
  }
  if (!dir.exists(diag_path2)){
    dir.create(diag_path2)
  }
  
  diag_path3 <- paste0(para_fields, "../diagnostic plots/run_", run, "/", Optimizer, "/")
  if (!dir.exists(diag_path3)){
    dir.create(diag_path3)
  }
  
  
  # Load spatial information about storage parameter
  true_para_field_df <- true_para_field(Test_number)
  # remove old states if existent
  tmp <- do.call(file.remove, list(
    list.files("GR4J_distributed/cdr/output", full.names = TRUE)))
  # 2. Training
  
  # CHANGE 2 KM DEFAULT FILE:
  default <- readLines("GR4J_distributed/input/defaults.txt")
  # set para file
  default[14] <- "para_Mur_GR4J_fsOptim.txt"
  default[15] <- "para_Mur_GR4J_true.txt"
  # set dates
  default[32] <- "2003 01 02 00 00" # start date
  default[33] <- "2009 08 31 00 00" # backup
  default[36] <- "2009 08 31 00 00"# end date
  default[37] <- "2012 12 31 00 00" # backup
  # set spin-up
  default[40] <- "241"
  default[41] <- "2433"
  # set Datafile
  default[7] <- "QObs_24h_synth_GR4J.txt" # used
  default[8] <- "QObs_24h.txt" # backup
  # set output type
  default[26] <- ifelse(Test_number %in% c(4.4, 4.5), "1", "0")
  writeLines(default, con = "GR4J_distributed/input/defaults.txt")
  
  # CHANGE 1 KM DEFAULT FILE:
  default <- readLines("GR4J_distributed_1km/input/defaults.txt")
  # set para file
  default[16] <- "para_Mur_GR4J_fsOptim.txt"
  default[17] <- "para_Mur_GR4J_true.txt"
  # set dates
  default[31] <- "2003 01 02 00 00" # start date
  default[32] <- "2009 08 31 00 00" # backup
  default[35] <- "2009 08 31 00 00" # end date
  default[36] <- "2012 12 31 00 00" # backup
  # set spin-up time
  default[39] <- "241"
  default[40] <- "2433"
  # set Datafile
  default[7] <- "QObs_24h_synth_GR4J.txt" # used
  default[8] <- "QObs_24h.txt" # backup
  # set output type
  default[25] <- "1"
  # write new default file
  writeLines(default, con = "GR4J_distributed_1km/input/defaults.txt")
  
  
  # FSO optimization
  source(paste0("Functions/", Optimizer, "_optimization_GR4J.R"), local = TRUE)
  
  results <- feather::read_feather(paste0("Test ", substr(Test_number, 1, 1),"/",
                                          "Test ", Test_number, "/", Optimizer,
                                          "_GR4J_optimization_", 
                                          Test_number, "_run", run, ".feather"))
  # define optimized functions 
  result_functions <- list(GR4Jx1 = tail(results$best_x1, 1),
                           GR4Jx3 = tail(results$best_x3, 1),
                           GR4Jx4 = tail(results$best_x4, 1))
  c("\n------------------------------------------\n")
  cat("Finished Function Space Optimization\nOptimized function:\n")
  for(i in 1:3) cat(names(result_functions)[i], ":", result_functions[[i]], "\n")
  cat("\nCreating and saving diagnostic plots in",
      paste0("Test ", substr(Test_number, 1, 1),"/",
             "Test ", Test_number, "/run_", run, "/", Optimizer))
  # 3. Testing
  
  # CHANGE 2 KM DEFAULT FILE:
  default <- readLines("GR4J_distributed/input/defaults.txt")
  # set para file
  default[14] <- "para_Mur_GR4J_fsOptim.txt"
  default[15] <- "para_Mur_GR4J_true.txt"
  # set dates
  default[32] <- "2003 01 02 00 00" # start date
  default[33] <- "2009 08 31 00 00" # nackup -> never used because spin-up is covering that
  default[36] <- "2012 12 31 00 00" # end date
  default[37] <- "2009 08 31 00 00" # backup
  # set spin-up
  default[40] <- "2433"
  default[41] <- "241"
  # set Datafile
  default[7] <- "QObs_24h_synth_GR4J.txt" # used
  default[8] <- "QObs_24h.txt" # backup
  writeLines(default, con = "GR4J_distributed/input/defaults.txt")
  
  # CHANGE 1 KM DEFAULT FILE:
  default <- readLines("GR4J_distributed_1km/input/defaults.txt")
  # set para file
  default[16] <- "para_Mur_GR4J_fsOptim.txt"
  default[17] <- "para_Mur_GR4J_true.txt"
  # set dates
  default[31] <- "2003 01 02 00 00" # start date
  default[32] <- "2009 08 31 00 00" # backup
  default[35] <- "2012 12 31 00 00" # end date
  default[36] <- "2009 08 31 00 00" # backup
  # set spin-up time
  default[39] <- "2433"
  default[40] <- "241"
  # set Datafile
  default[7] <- "QObs_24h_synth_GR4J.txt" # used
  default[8] <- "QObs_24h.txt" # backup
  # set output type
  default[25] <- "1"
  # write new default file
  writeLines(default, con = "GR4J_distributed_1km/input/defaults.txt")
  
  # Start test evaluation
  test_results <- evaluate_test_basins(test_functions = result_functions, 
                                       Optimizer = Optimizer,
                                       Test_number = Test_number,
                                       run = run,
                                       training_basins = training_basins, 
                                       test_basins = test_basins, 
                                       para = para, para_1km = para_1km,
                                       true_state_list = true_state_list)
  
  # 4. Diagnostic Plots & rasters
  library(ggplot2, quietly = TRUE)
  library(ggpubr, quietly = TRUE)
  for(parameter in c("x1", "x3", "x4")){
    plot_parameter <- list("x1" = "X1", "x3" = "X3", "x4" = "X4")
    raster_250m <- raster_from_tf(tf = result_functions[[eval(paste0("GR4J", parameter))]], 
                                  tf_bounds = parameter_bounds[[eval(paste0("GR4J", parameter))]])
    raster_1km <- raster_from_tf(tf = result_functions[[eval(paste0("GR4J", parameter))]], 
                                 tf_bounds = parameter_bounds[[eval(paste0("GR4J", parameter))]],
                                 aggregate = TRUE, km1 = TRUE, gmean_parameter = gmean_parameter)
    raster_2km <- raster_from_tf(tf = result_functions[[eval(paste0("GR4J", parameter))]], 
                                 tf_bounds = parameter_bounds[[eval(paste0("GR4J", parameter))]],
                                 aggregate = TRUE, km1 = FALSE, gmean_parameter = gmean_parameter)
    for(rasters in c("250m", "1km", "2km")){
      writeRaster(get(paste0("raster_", rasters)), 
                  paste0(para_fields2, Optimizer, "_", Test_number, "_", parameter, "_", rasters, ".asc"), 
                  format = "ascii", overwrite = TRUE)
      png(paste0(para_fields2,"/Plots/", Optimizer, "_", Test_number, "_", parameter, "_", rasters, ".png"),
          width = 1200, height = 800)
      plot(get(paste0("raster_", rasters)), 
           main = paste0("Parameter", plot_parameter[[parameter]], " - ", Optimizer))
      dev.off()
      
    }
    
    # get true parameter field
    assign("true_parameter", raster::raster(paste0("True parameters/", parameter, "_2km.asc")))
    
    # Plot True vs. predicted parameters
    
    plot_df <- data.frame("Observation" = values(true_parameter), 
                          "Prediction" = values(get(paste0("raster_2km"))))
    max_val <- max(plot_df, na.rm = TRUE)
    min_val <- min(plot_df, na.rm = TRUE)
    plot_df <- plot_df[!is.na(plot_df$Observation), ]
    correlation <- round(cor(plot_df, use = "pairwise.complete.obs")[2, 1], 2)
    
    plot_df <- plot_df[sample(nrow(plot_df), 500), ]
    ggplot(plot_df, aes(Observation, Prediction)) + geom_point(col = "cornsilk4") + 
      geom_smooth(method='lm', col = "darkgoldenrod2") + 
      labs(x = paste0("True ",  plot_parameter[[parameter]]),
           y = paste0("Predicted ",  plot_parameter[[parameter]])
      ) + annotate(geom = "text", 
                   x =  max_val, y = max_val,
                   label = paste0("R = ", correlation), 
                   size = 4, hjust = 1, vjust = 2) +
      ylim(min_val, max_val) + xlim(min_val, max_val) +
      theme_bw()
    
    ggsave(paste0(diag_path2, Optimizer, "/4_", Optimizer, "_", parameter, "_vs_true.png"),
           width = 7, height = 7, units = "in")
    
    plot_df_melt <- suppressWarnings(reshape2::melt(plot_df))
    ggplot(plot_df_melt, aes(value, fill = variable)) +
      geom_density(alpha = 0.4) +
      labs(
        x = plot_parameter[[parameter]],
        title = paste0("GR4J parameter ", plot_parameter[[parameter]]),
        subtitle = paste0(Optimizer, " optimization, ", "density estimation")
      ) + scale_fill_discrete(labels = c("true values",
                                         "predicted values"), 
                              name= "")
    ggsave(paste0(diag_path2, Optimizer, "/4_", Optimizer, "_", parameter, "_vs_true_density.png"),
           width = 7, height = 7, units = "in")
    
  }
  
  # Plot map with good and bad results
  if(Sys.info()[1] == "Linux"){
    path <- "/media/cfgrammar/Data/Dropbox/Diss/CF_Grammar/Data/spatial_predictors_mur/"
  } else {
    path <- "D:/Dropbox/Diss/CF_Grammar/Data/spatial_predictors_mur/"
  }
  # get raster objects of catchment
  nb <- raster::raster(paste0(path, "l0_nb2000_mur.asc"))
  
  test_result <- read.table(paste0(para_fields, "../testing/", 
                                   Optimizer, "_testing_", Test_number, "_run", run,  ".txt"),
                            skip = 26, header = TRUE, nrow = 112)
  # 1. Plot: Good vs. Bad
  good_results <- test_result[test_result$NSE >= 0.8,]
  bad_results <- test_result[test_result$NSE < 0.8,]
  # define good/bad/training as factor values
  good_bad <- nb
  values(good_bad)[values(nb) %in% good_results$Basin] <- 1
  values(good_bad)[values(nb) %in% bad_results$Basin] <- 2
  values(good_bad)[values(nb) %in% training_basins] <- 3
  values(good_bad)[is.na(values(nb))] <- NA
  
  # Make plotting data frame
  good_bad.p <- raster::rasterToPoints(good_bad)
  df <- data.frame(good_bad.p)
  # Make appropriate column headings
  colnames(df) <- c("Longitude", "Latitude", "Quality")
  df$Quality <- factor(df$Quality, levels = c(1, 2, 3), labels = c("good", "bad", "training"))
  # Now make the map
  ggplot(data = df, aes(y = Latitude, x = Longitude)) +
    geom_raster(aes(fill = Quality)) +
    theme_bw() +
    coord_equal() +
    scale_fill_manual(values = c("chartreuse2", "brown1", "grey"), 
                      labels = c("NSE >= 0.8", "NSE < 0.8", "Training"),
                      name = "", drop = FALSE) +
    labs(
      title = paste0(Optimizer, " optimization"),
      subtitle = "Testing time period 2009 - 2012"
    ) +
    ggsave(paste0(diag_path2, Optimizer, "/1_result_map_", Optimizer, ".png"),
           width = 7, height = 7, units = "in")
  
  
  # 2. plot: NSE map
  nse_basins <- nb
  for (i in unique(nb)){
    values(nse_basins)[values(nb) == i] <- test_result$NSE[test_result$Basin == i]
  }
  nse_basins.p <- raster::rasterToPoints(nse_basins)
  df2 <- data.frame(nse_basins.p)
  #Make appropriate column headings
  colnames(df2) <- c("Longitude", "Latitude", "NSE")
  df2$NSE <- round(df2$NSE, 4)
  #Now make the map
  min_NSE <- min(df2$NSE)
  if(min_NSE < -1){
    breaks <- c(1, 0, -1, min_NSE)
    cols <- c("black", "brown1", "chartreuse3")
  } else {
    breaks <- c(1, 0, -1)
    cols <- c("brown1", "chartreuse3")
  }
  if(min_NSE > 0){
    nse_lim <- c(0, 1)
  } else {
    nse_lim <- c(min_NSE, 1)
  }
  ggplot(data = df2, aes(y = Latitude, x = Longitude)) +
    geom_raster(aes(fill=NSE)) +
    theme_bw() +
    coord_equal() +
    scale_fill_gradientn("NSE", limits = nse_lim, breaks = breaks,
                         colors = cols) +
    labs(
      title = paste0(Optimizer, " optimization"),
      subtitle = "Testing time period 2009 - 2012"
    ) +
    ggsave(paste0(diag_path2, Optimizer, "/2_NSE_map_", Optimizer, ".png"),
           width = 7, height = 7, units = "in")
  
  # 3. plot: Elevation vs. good and bad predictions
  compare_results_df <- load_sp_mur(scale = FALSE, na.approx = FALSE, 
                                    only_training_basins = FALSE, full_dataset = FALSE)
  #compare_results_df <- aggregate(elevation ~ nb, compare_results_df, mean)
  compare_results_df <- merge(compare_results_df, test_result, by.x = "nb", by.y= "Basin")
  compare_results_df$good_bad <- "NSE > 0.8"
  compare_results_df$good_bad[compare_results_df$nb %in% bad_results$Basin] <- "NSE <= 0.8"
  compare_results_df$good_bad[compare_results_df$nb %in% training_basins] <- "Training"
  compare_results_df$good_bad <- factor(compare_results_df$good_bad, 
                                        levels = c("NSE <= 0.8", "NSE > 0.8", "Training"))
  
  ggplot(compare_results_df, aes(good_bad, elevation)) + geom_boxplot() +
    scale_x_discrete("good_bad", drop = FALSE) +
    labs(
      x = "", y = "mean zone elevation",
      title = paste0(Optimizer, " optimization"),
      subtitle = "Testing time period 2009 - 2012"
    ) +
    ggsave(paste0(diag_path2, Optimizer, "/3_elevation_vs_prediction_", 
                  Optimizer, ".png"),
           width = 7, height = 7, units = "in")
  
  # 4. Compare Training/test NSE for test time period
  test_result$train_test <- "Test"
  test_result$train_test[test_result$Basin %in% training_basins] <- "Training"
  
  
  ggplot(test_result, aes(factor(train_test), NSE)) + geom_boxplot() +
    labs(
      x = "", y = "NSE of testing time period",
      title = paste0(Optimizer, " optimization NSE"),
      subtitle = "Testing time period 2009 - 2012"
    ) +
    ggsave(paste0(diag_path2, Optimizer, "/4_NSE_distributions_", 
                  Optimizer, ".png"),
           width = 7, height = 7, units = "in")
  
  cat("\nDone!\n")
}

# Load FSO setup
FSO_setup <- function(){
  cat("\nCase study is setup as defined in Functions/case_study_setup.\n")
  # Load setup
  source("Functions/case_study_setup.R")
  # Load function space VAE 
  source("Functions/FSO_VAE_generator.R")
  cat("\nFSO setup complete!\n")
}




