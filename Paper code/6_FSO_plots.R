# FSO performance plots
# Moritz Feigl, 2019
#


setwd("FSO_paper")


source("Functions/FSO_functions.R")
transfer_functions <- list(GR4Jx1 = "0.5 + evi*1.5 + exp(bdim)*0.9", # production store
                           GR4Jx2 = "0",          # interaction direct runoff and routing store
                           GR4Jx3 = "-1.3  - log(slope)",  # routing store
                           GR4Jx4 = "-log(hand)*0.2 - slope*1.5 - 1.5") # length of UH
training_basins <- c(1, 3, 4, 8, 27, 30, 35, 38, 46, 47, 50, 51, 57, 60, 69, 72, 75, 84, 90, 
                     93, 95, 96, 98, 100, 101, 103, 111) # 27 training basins
test_basins <- c(2, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                 21, 22, 23, 24, 25, 26, 28, 29, 31, 32, 33, 34, 36, 37, 39, 40,
                 41,  42, 43, 44, 45, 48, 49, 52, 53, 54, 55, 56, 58, 59, 61, 62,
                 63, 64, 65, 66, 67, 68, 70, 71, 73, 74, 76, 77, 78, 79, 80, 81,
                 82, 83, 85, 86, 87, 88, 89, 91, 92, 94, 97, 99, 102, 104, 105, 106,
                 107, 108, 109, 110, 112)
bounds <- read.table("GR4J_distributed/parameter.txt",
                     skip = 1, sep = "\t", header = FALSE)
parameter_bounds <- list(GR4Jx1 = as.numeric(bounds[7, c(3, 4)]), # production store
                         GR4Jx2 = as.numeric(bounds[8, c(3, 4)]),          # transboundary flow
                         GR4Jx3 = as.numeric(bounds[9, c(3, 4)]),  # routing store
                         GR4Jx4 = as.numeric(bounds[10, c(3, 4)])) #length of the UH


FSO_plot <- function(Test_number){
  # Plot all test results
  
  # load spatial predictors
  l0 <- load_sp_mur(scale = TRUE, na.approx = FALSE, 
                    only_training_basins = FALSE)
  # define test path
  test_path <- paste0("Test ", substr(Test_number, 1, 1),"/",
                      "Test ", Test_number, "/")
  library(tidyverse)
  library(feather)
  # # find relevant files
  files <- list.files(test_path)
  files <- files[grep(".feather",files)]
  files <- files[grep(Test_number,files)]
  
  # load first result
  data_df <- feather::read_feather(paste0(test_path, files[1]))
  data_df$run <- substr(files[1], nchar(files[1]) - 8, nchar(files[1])- 8)
  data_df <- data_df[data_df$full_loss != -9999, ]
  data_df$cum_it <- cumsum(data_df$n_iteration_used)
  data_df$method <- strsplit(files[1], "_")[[1]][1]
  
  for(i in 2:length(files)) {
    data <- feather::read_feather(paste0(test_path, files[i]))
    data <- data[data$full_loss != -9999, ]
    data$run <- substr(files[i], nchar(files[i]) - 8, nchar(files[i])- 8)
    data$cum_it <- cumsum(data$n_iteration_used)
    data$method <- strsplit(files[i], "_")[[1]][1]
    data_df <- rbind(data_df, data)
  }
  
  names(data_df)[12] <- "Method"
  
  # Check if really all optimizer results are present
  if(sum(unique(data_df$Method) %in% c("GA", "DDS", "PSO")) != 3) {
    stop("Not all Optimizer results are found. Check if test ", Test_number, 
         " was run with GA, DDS and PSO.")
  }
  
  max_Loss <- c(
    max(data_df[data_df$Method == "GA" , "full_loss"]),
    max(data_df[data_df$Method == "DDS" , "full_loss"]),
    max(data_df[data_df$Method == "PSO" , "full_loss"]))
  max_NSEs <- unlist(c(
    data_df[data_df$Method == "GA" & data_df$full_loss == max_Loss[1], "NSE"][1, 1],
    data_df[data_df$Method == "DDS" & data_df$full_loss == max_Loss[2], "NSE"][1, 1],
    data_df[data_df$Method == "PSO" & data_df$full_loss == max_Loss[3], "NSE"][1, 1]))
  
  optimized_funs <- list(data.frame(best_x1 = transfer_functions[["GR4Jx1"]],
                                    best_x3 = transfer_functions[["GR4Jx3"]],
                                    best_x4 = transfer_functions[["GR4Jx4"]], stringsAsFactors = FALSE),
                         as.data.frame(data_df[data_df$Method == "GA", ][
                           which.max(as.matrix(data_df[data_df$Method == "GA", "full_loss"])), 1:3], 
                           stringsAsFactors = FALSE),
                         as.data.frame(data_df[data_df$Method == "DDS", ][
                           which.max(as.matrix(data_df[data_df$Method == "DDS", "full_loss"])), 1:3], 
                           stringsAsFactors = FALSE),
                         as.data.frame(data_df[data_df$Method == "PSO", ][
                           which.max(as.matrix(data_df[data_df$Method == "PSO", "full_loss"])), 1:3], 
                           stringsAsFactors = FALSE))
  optimized_funs <- do.call(rbind, optimized_funs)
  names(optimized_funs) <- c("X1", "X3", "X4")
  optimized_funs <- cbind(Method = c("True functions", "GA", "DDS", "PSO"), optimized_funs, stringsAsFactors = FALSE)
  
  # best run
  example_runs <- c(
    as.character(data_df[data_df$Method == "GA", ][which.max(as.matrix(data_df[data_df$Method == "GA", "full_loss"])), "run"][1,]),
    as.character(data_df[data_df$Method == "DDS", ][which.max(as.matrix(data_df[data_df$Method == "DDS", "full_loss"])), "run"][1,]),
    as.character(data_df[data_df$Method == "PSO", ][which.max(as.matrix(data_df[data_df$Method == "PSO", "full_loss"])), "run"][1,]))
  
  # iteration needed until best result
  nr_it_needed <- c(
    as.character(data_df[data_df$Method == "GA", ][which.max(as.matrix(data_df[data_df$Method == "GA", "full_loss"])), "cum_it"][1,]),
    as.character(data_df[data_df$Method == "DDS", ][which.max(as.matrix(data_df[data_df$Method == "DDS", "full_loss"])), "cum_it"][1,]),
    as.character(data_df[data_df$Method == "PSO", ][which.max(as.matrix(data_df[data_df$Method == "PSO", "full_loss"])), "cum_it"][1,]))
  
  # best functions simplified
  optimized_funs[, 2] <- c(optimized_funs[1, 2],
                           Deriv::Simplify(optimized_funs[2, 2]),
                           Deriv::Simplify(optimized_funs[3, 2]),
                           Deriv::Simplify(optimized_funs[4, 2]))
  
  optimized_funs[, 3] <- c(optimized_funs[1, 3],
                           Deriv::Simplify(optimized_funs[2, 3]),
                           Deriv::Simplify(optimized_funs[3, 3]),
                           Deriv::Simplify(optimized_funs[4, 3]))
  optimized_funs[, 4] <- c(optimized_funs[1, 4],
                           Deriv::Simplify(optimized_funs[2, 4]),
                           Deriv::Simplify(optimized_funs[3, 4]),
                           Deriv::Simplify(optimized_funs[4, 4]))
  
  # iterations needed to get NSE > 0.95
  it_until_larger_nse <- unlist(c(
    min(data_df[data_df$Method == "GA" & data_df$NSE > 0.95, "cum_it"][, 1]),
    min(data_df[data_df$Method == "DDS" & data_df$NSE > 0.95, "cum_it"][, 1]),
    min(data_df[data_df$Method == "PSO" & data_df$NSE > 0.95, "cum_it"][, 1])))
  
  # test losses
  ga_test <- read.table(paste0(test_path, "/testing/GA_testing_", Test_number, 
                               "_run", example_runs[1],  ".txt"),
                        skip = 26, header = TRUE, nrow = 112)
  dds_test <- read.table(paste0(test_path, "/testing/DDS_testing_", Test_number, 
                                "_run", example_runs[2],  ".txt"),
                         skip = 26, header = TRUE, nrow = 112)
  pso_test <- read.table(paste0(test_path, "/testing/PSO_testing_", Test_number, 
                                "_run", example_runs[3],  ".txt"),
                         skip = 26, header = TRUE, nrow = 112)
  
  ga_test_nse <- mean(ga_test$NSE[ga_test$Basin %in% test_basins])
  dds_test_nse <- mean(dds_test$NSE[dds_test$Basin %in% test_basins])
  pso_test_nse <- mean(pso_test$NSE[pso_test$Basin %in% test_basins])
  
  # number of correct l0 of all possible correct l0
  get_l0s <- function(tfs, l0){
    l0s <- lapply(tfs, function_splitter)
    l0s <- lapply(l0s, function(x) x[(x %in% names(l0))])
    return(l0s)
  }
  l0s_true <- get_l0s(as.list(optimized_funs[1, -1]), l0)
  l0s_GA <- get_l0s(as.list(optimized_funs[2, -1]), l0)
  l0s_DDS <- get_l0s(as.list(optimized_funs[3, -1]), l0)
  l0s_PSO <- get_l0s(as.list(optimized_funs[4, -1]), l0)
  
  total_number_of_l0s <- length(unlist(l0s_true))
  how_many_l0_true <- function(l0s_predict, l0s_true){
    hits <- integer()
    for(i in 1:3){
      hits[i] <- sum(unique(l0s_predict[[i]]) %in% unique(l0s_true[[i]]))
    }
    return(sum(hits))
  }
  l0_hits <- c(how_many_l0_true(l0s_GA, l0s_true), 
               how_many_l0_true(l0s_DDS, l0s_true), 
               how_many_l0_true(l0s_PSO, l0s_true))
  
  # prepare results table
  results <- data.frame(optimized_funs,
                        "True l0" = c("", 
                                      paste0(l0_hits, "/", total_number_of_l0s)),
                        "Iterations" = c(" ", nr_it_needed),
                        #"Iterations until NSE > 0.95" = c(" ", it_until_larger_nse),
                        "Train NSE" = c(" ", round(max_NSEs, 3)),
                        #"Train Loss" = c(" ", round(max_Loss, 3)),
                        "Test NSE" = c(" ", as.character(round(ga_test_nse, 3)), 
                                       as.character(round(dds_test_nse, 3)), 
                                       as.character(round(pso_test_nse, 3))), check.names = FALSE)
  
  # 1. Plot: NSE vs. n_iter
  library(ggplot2)
  library(gtable)
  library(grid)
  library(gridExtra)
  # melt for plotting
  sub_df <- data_df[, c("cum_it", "NSE", "Method", "run")]
  m_sub_df <- suppressWarnings(reshape2::melt(sub_df, id.vars = c("Method", "run", "cum_it")))
  m_sub_df$cum_it <- as.integer(m_sub_df$cum_it)
  names(m_sub_df)[1] <- "Optimizer"
  if(Test_number == 1.2){
    ylim_low <- 0.98
  } else {
    ylim_low <- 0.8
  }
  p1 <-  ggplot(m_sub_df[m_sub_df$run == 0, ], aes(x = cum_it, y = value, col = Optimizer)) + 
    geom_line(alpha = 0.8) +
    geom_line(data = m_sub_df[m_sub_df$run == 1, ], 
              aes(x = cum_it, y = value, col = Optimizer), alpha = 0.8) +
    geom_line(data = m_sub_df[m_sub_df$run == 2, ], 
              aes(x = cum_it, y = value, col = Optimizer), alpha = 0.8) +
    geom_line(data = m_sub_df[m_sub_df$run == 3, ], 
              aes(x = cum_it, y = value, col = Optimizer), alpha = 0.8) +
    geom_line(data = m_sub_df[m_sub_df$run == 4, ], 
              aes(x = cum_it, y = value, col = Optimizer), alpha = 0.8) +
    geom_line(data = m_sub_df[m_sub_df$run == 5, ], 
              aes(x = cum_it, y = value, col = Optimizer), alpha = 0.8) +
    geom_hline(aes(yintercept = 1)) +
    xlab("Number of iterations") + ylab("NSE") +
    coord_cartesian(ylim = c(ylim_low, 1)) +
    scale_x_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000),
                       labels = c("0", "500", "1000", "1500", "2000", "2500", "3000"),
                       limits = c(0, 3000)) +
    theme(legend.text = element_text(size = 11),
          legend.title=element_text(size=12),
          axis.text=element_text(size=10),
          axis.title=element_text(size=13,face="bold")) + theme_minimal()
  
  mytheme <- ttheme_default(base_size=10, core=list(fg_params=list(parse=TRUE)))
  
  
  for(i in 2:4) {
    results[, i] <- gsub("*", "%*%", results[, i], fixed = TRUE)
    results[, i] <- gsub("/", "%/%", results[, i], fixed = TRUE)
  }
  g <- tableGrob(results, rows = NULL, theme = mytheme)
  g <- gtable_add_grob(g,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 2, b = nrow(g), l = 1, r = ncol(g))
  g <- gtable_add_grob(g,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 1, l = 1, r = ncol(g))
  hlay <- rbind(c(1, 1), c(2, 2))
  
  cat("Plot is saved under:\n", paste0(Test_number,"_Optimization_results.png"), "\n")
  ggsave(file = paste0(Test_number,"_Optimization_results.png"),
         arrangeGrob(p1, g, #g_real, g,
                     heights = c(8, 8),
                     layout_matrix = hlay,
                     top = textGrob("", gp=gpar(fontsize=18,font=8))),
         height = 18, width = 30, unit = "cm") #saving plot
  
  # write out which runs were the best
  cat("Overview of best runs per optimizer are saved under:\n",  paste0(test_path, "example_runs.txt"), "\n")
  write.table(data.frame(Optimizer = c("GA", "DDS", "PSO"),
                         example_run = example_runs),
              paste0(test_path, "example_runs.txt"))
  
  img <- png::readPNG(paste0(Test_number,"_Optimization_results.png"), native = FALSE, info = FALSE)
  grid::grid.raster(img)
}
sapply(c(1.2, 2.4, 2.5, 2.6, 4.4, 4.5, 4.6), FSO_plot)



single_FSO_plot <- function(Test_number, Optimizer = "DDS", ylim_low){
  # Plot for a single FSO result
  
  # load spatial predictors
  l0 <- load_sp_mur(scale = TRUE, na.approx = FALSE, 
                    only_training_basins = FALSE)
  # define test path
  test_path <- paste0("Test ", substr(Test_number, 1, 1),"/",
                      "Test ", Test_number, "/")
  
  library(tidyverse)
  library(feather)
  # # find relevant files
  files <- list.files(test_path)
  files <- files[grep(".feather", files)]
  files <- files[grep(Optimizer, files)]
  files <- files[grep(Test_number,files)]
  
  # load first result
  data_df <- feather::read_feather(paste0(test_path, files[1]))
  data_df$run <- substr(files[1], nchar(files[1]) - 8, nchar(files[1])- 8)
  data_df <- data_df[data_df$full_loss != -9999, ]
  data_df$cum_it <- cumsum(data_df$n_iteration_used)
  data_df$method <- strsplit(files[1], "_")[[1]][1]
  
  if(length(files) > 1){
    for(i in 2:length(files)) {
      data <- feather::read_feather(paste0(test_path, files[i]))
      data <- data[data$full_loss != -9999, ]
      data$run <- substr(files[i], nchar(files[i]) - 8, nchar(files[i])- 8)
      data$cum_it <- cumsum(data$n_iteration_used)
      data$method <- strsplit(files[i], "_")[[1]][1]
      data_df <- rbind(data_df, data)
    }
  }
  names(data_df)[12] <- "Method"
  ## Timeseries plot
  library(ggplot2)
  library(gtable)
  library(grid)
  library(gridExtra)
  # melt for plottin
  sub_df <- data_df[, c("cum_it", "NSE", "full_loss", "Method", "run")]
  m_sub_df <- suppressWarnings(reshape2::melt(sub_df, id.vars = c("Method", "run", "cum_it")))
  m_sub_df$cum_it <- as.integer(m_sub_df$cum_it)
  names(m_sub_df)[1:2] <- c("Optimizer", "Run")
  
  p1 <-  ggplot(m_sub_df, aes(x = cum_it, y = value, col = Run, linetype = variable)) +
    geom_line() +
    geom_hline(yintercept = 1) +
    coord_cartesian(ylim = c(ylim_low, 1)) +
    xlab("\nNumber of iterations") + ylab("NSE\n") +
    scale_x_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000),
                       labels = c("0", "500", "1000", "1500", "2000", "2500", "3000"),
                       limits = c(0, 3000)) +
    scale_linetype_manual(name = "Loss type", values = c(1, 2), labels = c("Mean basin\nNSE", "Optimization\nNSE")) +
    theme(axis.text=element_text(size=10),
      axis.title=element_text(size=13,face="bold"))+ theme_minimal() +
    ggsave(file = paste0("Result Plots/", Test_number, "/", Test_number, "_FSO_training.png"), width = 9, height = 7, units = "in")
  print(p1)
  
  # Plot with only mean NSE
  p2 <-  ggplot(m_sub_df[m_sub_df$variable == "NSE", ], aes(x = cum_it, y = value, col = Run)) +
    geom_line() +
    geom_hline(yintercept = 1) +
    coord_cartesian(ylim = c(ylim_low, 1)) +
    xlab("\nNumber of iterations") + ylab("NSE\n") +
    scale_x_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000),
                       labels = c("0", "500", "1000", "1500", "2000", "2500", "3000"),
                       limits = c(0, 3000)) +
    theme(legend.text = element_text(size = 100),
      legend.title=element_text(size=17),
      axis.text=element_text(size=25),
      axis.title=element_text(size=13,face="bold"))+ theme_minimal(base_size = 20) +
    ggsave(file = paste0("Result Plots/", Test_number, "/", Test_number, "_FSO_training_simple.png"), width = 9, height = 7, units = "in")
}

FSO_result_table <- function(Test_number, example_run){
  # load spatial predictors
  l0 <- load_sp_mur(scale = TRUE, na.approx = FALSE, 
                    only_training_basins = FALSE)
  # define test path
  test_path <- paste0("Test ", substr(Test_number, 1, 1),"/",
                      "Test ", Test_number, "/")
  
  library(tidyverse)
  library(feather)
  # # find relevant files
  files <- list.files(test_path)
  files <- files[grep(".feather", files)]
  files <- files[grep("DDS", files)]
  files <- files[grep(Test_number,files)]
  
  # load first result
  data_df <- feather::read_feather(paste0(test_path, files[1]))
  data_df$run <- substr(files[1], nchar(files[1]) - 8, nchar(files[1])- 8)
  data_df <- data_df[data_df$full_loss != -9999, ]
  data_df$cum_it <- cumsum(data_df$n_iteration_used)
  data_df$method <- strsplit(files[1], "_")[[1]][1]
  
  if(length(files) > 1){
    for(i in 2:length(files)) {
      data <- feather::read_feather(paste0(test_path, files[i]))
      data <- data[data$full_loss != -9999, ]
      data$run <- substr(files[i], nchar(files[i]) - 8, nchar(files[i])- 8)
      data$cum_it <- cumsum(data$n_iteration_used)
      data$method <- strsplit(files[i], "_")[[1]][1]
      data_df <- rbind(data_df, data)
    }
  }
  names(data_df)[12] <- "Method"
  max_Loss <- max(data_df[data_df$Method == "DDS" & data_df$run == example_run, "full_loss"])
  max_NSEs <- unlist(data_df[data_df$Method == "DDS" & data_df$run == example_run 
                             & data_df$full_loss == max_Loss, "NSE"][1, 1])
  
  optimized_funs <- list(data.frame(best_x1 = transfer_functions[["GR4Jx1"]],
                                    best_x3 = transfer_functions[["GR4Jx3"]],
                                    best_x4 = transfer_functions[["GR4Jx4"]], stringsAsFactors = FALSE),
                         as.data.frame(data_df[data_df$Method == "DDS", ][
                           which.max(as.matrix(data_df[data_df$Method == "DDS", "full_loss"])), 1:3], 
                           stringsAsFactors = FALSE))
  optimized_funs <- do.call(rbind, optimized_funs)
  names(optimized_funs) <- c("X1", "X3", "X4")
  optimized_funs <- cbind(Method = c("True functions", "DDS"), optimized_funs, stringsAsFactors = FALSE)
  
  # iteration needed until best result
  nr_it_needed <- as.character(data_df[data_df$Method == "DDS", ][which.max(as.matrix(data_df[data_df$Method == "DDS", "full_loss"])), "cum_it"][1,])
  
  # best functions simplified
  optimized_funs[, 2] <- c(optimized_funs[1, 2],
                           Deriv::Simplify(optimized_funs[2, 2]))
  
  optimized_funs[, 3] <- c(optimized_funs[1, 3],
                           Deriv::Simplify(optimized_funs[2, 3]))
  optimized_funs[, 4] <- c(optimized_funs[1, 4],
                           Deriv::Simplify(optimized_funs[2, 4]))
  
  # iterations needed to get NSE > 0.95
  it_until_larger_nse <- unlist(min(data_df[data_df$Method == "DDS" & data_df$NSE > 0.95 & data_df$run == example_run, "cum_it"][, 1]))
  
  # test losses
  dds_test <- read.table(paste0(test_path, "/testing/DDS_testing_", Test_number, 
                                "_run", example_run,  ".txt"),
                         skip = 26, header = TRUE, nrow = 112)
  
  dds_test_nse <- mean(dds_test$NSE[dds_test$Basin %in% test_basins])
  
  # number of correct l0 of all possible correct l0
  get_l0s <- function(tfs, l0){
    l0s <- lapply(tfs, function_splitter)
    l0s <- lapply(l0s, function(x) x[(x %in% names(l0))])
    return(l0s)
  }
  l0s_true <- get_l0s(as.list(optimized_funs[1, -1]), l0)
  l0s_DDS <- get_l0s(as.list(optimized_funs[2, -1]), l0)
  
  total_number_of_l0s <- length(unlist(l0s_true))
  how_many_l0_true <- function(l0s_predict, l0s_true){
    hits <- integer()
    for(i in 1:3){
      hits[i] <- sum(unique(l0s_predict[[i]]) %in% unique(l0s_true[[i]]))
    }
    return(sum(hits))
  }
  l0_hits <- how_many_l0_true(l0s_DDS, l0s_true)
  
  # prepare results table
  results <- data.frame(optimized_funs,
                        "True l0" = c("", 
                                      paste0(l0_hits, "/", total_number_of_l0s)),
                        "Iterations" = c(" ", nr_it_needed),
                        #"Iterations until NSE > 0.95" = c(" ", it_until_larger_nse),
                        "Train NSE" = c(" ", round(max_NSEs, 3)),
                        #"Train Loss" = c(" ", round(max_Loss, 3)),
                        "Test NSE" = c(" ", as.character(round(dds_test_nse, 3))), check.names = FALSE)
  
  for(i in 2:4) {
    results[, i] <- gsub("*", "%*%", results[, i], fixed = TRUE)
    results[, i] <- gsub("/", "%/%", results[, i], fixed = TRUE)
  }
  results[2, 1] <- "FSO"
  names(results)[1] <- ""
  # Table grob
  mytheme <- ttheme_default(base_size=10, core=list(fg_params=list(parse=TRUE)))
  g <- tableGrob(results, rows = NULL, theme = mytheme)
  g <- gtable_add_grob(g,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 2, b = nrow(g), l = 1, r = ncol(g))
  g <- gtable_add_grob(g,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 1, l = 1, r = ncol(g))
  hlay <- rbind(c(1, 1))
  
  ggsave(file = paste0("Result Plots/", Test_number, "/", Test_number, "_FSO_table.png"), 
         width = 15, height = 7, units = "in",
         plot = g)
  cat("\nsaved table in ",  paste0("Result Plots/", Test_number, "/", Test_number, "_FSO_table.png"))
  plot(g)
  
}

single_FSO_dist <- function(Test_number, example_run){
  # Diagnostic Plots & rasters
  library(ggplot2, quietly = TRUE)
  library(ggpubr, quietly = TRUE)
  test_path <- paste0("Test ", substr(Test_number, 1, 1),"/",
                      "Test ", Test_number, "/")
  
  for(parameter in c("x1", "x3", "x4")){
    plot_parameter <- list("x1" = "X1", "x3" = "X3", "x4" = "X4")
    # get true parameter field
    assign("true_parameter", raster::raster(paste0("../True parameters/", parameter, "_2km.asc")))
    assign("pred_parameter", 
           raster::raster(
             paste0(test_path, "parameter fields/run_", example_run, "/DDS_", Test_number, "_", parameter, "_2km.asc")))
    
    # Plot True vs. predicted parameters
    plot_df <- data.frame("Observation" = values(true_parameter), 
                          "Prediction" = values(pred_parameter))
    max_val <- max(plot_df, na.rm = TRUE)
    min_val <- min(plot_df, na.rm = TRUE)
    plot_df <- plot_df[!is.na(plot_df$Observation), ]
    correlation <- round(cor(plot_df, use = "pairwise.complete.obs")[2, 1], 2)
    
    plot_df <- plot_df[sample(nrow(plot_df), 500), ]
    ggplot(plot_df, aes(Observation, Prediction)) + geom_point(col = "cornsilk4") + 
      geom_smooth(method='lm', col = "darkorange1", se = FALSE) + 
      labs(
        x = paste0("True ",  plot_parameter[[parameter]]),
        y = paste0("Predicted ",  plot_parameter[[parameter]])
      ) + annotate(geom = "text", 
                   x =  (max_val - min_val)/2 + min_val, y = max_val,
                   label = paste0("R = ", correlation), 
                   size = 5, hjust = 0.7, vjust = 1) +
      ylim(min_val, max_val) + xlim(min_val, max_val) +
      theme_minimal() +
      ggsave(paste0("Result Plots/", Test_number, "/", Test_number, "_FSO_", parameter, "scatter.png"),
             width = 7, height = 7, units = "in")
    
    plot_df_melt <- suppressWarnings(reshape2::melt(plot_df))
    p_dens <- ggplot(plot_df_melt, aes(value, fill = variable)) +
      geom_density(alpha = 0.4) +
      labs(
        x = paste0("\n", plot_parameter[[parameter]]),
        y = "Density\n"
      ) + scale_fill_manual(labels = c("True values",
                                       "Predicted values"), 
                            name= "", values = c("deepskyblue2", "darkorange1")) + theme_minimal()
    means <- aggregate(value ~ variable, plot_df_melt, mean)
    dens <- ggplot_build(ggplot(plot_df_melt, aes(value, fill = variable)) +
                           geom_density(alpha = 0.4))$data[[1]] 
    dens_obs <- dens[dens$group == 1, ]
    mean_obs_height <- dens_obs[which.min(abs(dens_obs$x - means[1, 2])), "density"]
    dens_pred <- dens[dens$group == 2, ]
    mean_pred_height <- dens_pred[which.min(abs(dens_pred$x - means[2, 2])), "density"]
    mean_df <- cbind(means, height = c(mean_obs_height, mean_pred_height))
    p_dens +
      geom_segment(mean_df, mapping = aes(x = value, xend = value, y = 0, yend = height,
                                          col = variable),
                   linetype = "solid", size = 1) +
      scale_color_manual(name = "", labels = c("True mean",
                                               "Predicted mean"),
                         values = c("deepskyblue4", "darkorange3")) +
      ggsave(paste0("Result Plots/", Test_number, "/", Test_number, "_FSO_", parameter, "density.png"),
             width = 7, height = 7, units = "in")
  }
}

ylims_low <- c(0.98, rep(0.8, 6))
tests <- c(1.2, 2.4, 2.5, 2.6, 4.4, 4.5, 4.6)
# Find best run in testing

for(test in 1:7){
  if(tests[test] == 1.2) example_run <- 5
  if(tests[test] == 2.4) example_run <- 1
  if(tests[test] == 2.5) example_run <- 4
  if(tests[test] == 2.6) example_run <- 5
  if(tests[test] == 4.4) example_run <- 2
  if(tests[test] == 4.5) example_run <- 5
  if(tests[test] == 4.6) example_run <- 2
  single_FSO_plot(tests[test], ylim_low = ylims_low[test])
  FSO_result_table(tests[test], example_run = example_run)
  single_FSO_dist(tests[test], example_run = example_run)
}

