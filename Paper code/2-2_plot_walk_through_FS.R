#
# Script for plotting walk through function space
#


# Setup ----------------------------------------------------------------------------------
setwd("FSO_paper")
source("Functions/FSO_VAE_generator.R")
source("Functions/FSO_functions.R")

generator_data <- read_feather(
  "Data/generator_data_simple_10numerics_wrecalc_allBasins_no_extremes.feather")
index2word <- as.character(read_feather("Data/index2word_10numerics.feather"))
tf_dist <- read_feather(
  "Data/functions_simple_10_numerics_Distribution_indiv_scale_wrecalc_allBasins_no_extremes.feather")
# scale distribution values
min_max_scale <- function(x) round((x-min(x))/(max(x)-min(x)), 11)
dist_scaled <- apply(tf_dist[,-c(1, 2, 12)], 2, min_max_scale)
# Train/val split
train_ind <- readRDS("Data/train_ind_no_extremes.rds")
x_train_tf <- as.matrix(generator_data[train_ind, ])
x_val_tf <- as.matrix(generator_data[-train_ind, ])
y_train_tf <- array(x_train_tf, dim = c(nrow(x_train_tf),
                                        ncol(x_train_tf),
                                        1))
y_val_tf <- array(x_val_tf, dim = c(nrow(x_val_tf),
                                    ncol(x_val_tf),
                                    1))
train_dist <- as.matrix(dist_scaled[train_ind, ])
val_dist <- as.matrix(dist_scaled[-train_ind, ])

# Walk through Function Space ------------------------------------------------------------
weights <- c(seq(0, 1, 0.1))
set.seed(999)
function_idx <- sample(nrow(x_val_tf), 2)
proto_dist <- val_dist[function_idx, ]
proto_tfs <- x_val_tf[function_idx, ]
proto_tfs_words <- apply(proto_tfs, 1, ind2word, index2word)
proto_tfs_encoded <- predict(encoder, list(proto_tfs, proto_dist))
# depending on weight w, make function predictions
new_tfs <- matrix(NA, ncol = ncol(proto_tfs_encoded), nrow = length(weights))
for(i in seq_along(weights)){
  new_tfs[i, ] <- proto_tfs_encoded[1,] * (1-weights[i]) + proto_tfs_encoded[2,] * weights[i]
}
new_tfs_pred <- predict(generator_tf, new_tfs, batch_size = 1)
new_tfs_function <- apply(new_tfs_pred, 1, tf_prediction)
data.frame(weight_f2 = weights, new_tfs_function)
write.table(data.frame(weight_f2 = weights, new_tfs_function), "walk_through_function_space.txt",
            row.names = FALSE)
# Plot distribution change
sent_encoded <- predict(encoder, list(proto_tfs, proto_dist))

dist_prediction <- as.data.frame(predict(generator_dist, new_tfs, batch_size = 1))
dist_prediction <- as.data.frame(apply(dist_prediction, 2, rescale, from = c(0, 1), to = c(-11, 11)))
dist_prediction <- cbind(names = c("start", 
                                   paste0("step ", 1:(nrow(dist_prediction)-2)),
                                   "end"), 
                         dist_prediction, check_names = FALSE)
colnames(dist_prediction) <- c("name", 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
dist_prediction$name <- factor(dist_prediction$name, 
                               levels = c("start", 
                                          paste0("step ", 1:(nrow(dist_prediction)-2)),
                                          "end"))
dist_prediction[, -c(1, 11)] <- min_max_scale(dist_prediction[, -c(1, 11)])
plot_data <- reshape2:: melt(dist_prediction)
plot_data$variable <- as.numeric(as.character(plot_data$variable))
library(ggplot2)
cc <- scales::seq_gradient_pal("green", "red", "Lab")(seq(0,1,length.out = nrow(dist_prediction)))
ggplot(plot_data, aes(x = value, y = variable, col = name)) + geom_line() + 
  scale_colour_manual(values = cc, name = "") + xlab("scaled parameter values") + 
  ylim(0, 1) + xlim(0, 1) +
  ylab("cumulative probability") + theme_classic() +
  ggsave("walk_through_FS.png", width = 7, height = 7, units = "in")
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
