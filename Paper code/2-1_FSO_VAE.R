# Training FSO Variational Autoencoder
# Moritz Feigl, 2019
#
# Load preprocessed data -----------------------------------------------------------------
setwd("FSO_paper")
library(feather)
library(keras)
library(tensorflow)
generator_data <- read_feather(
  "Data/generator_data_simple_10numerics_wrecalc_allBasins_no_extremes.feather"
  )
index2word <- as.character(read_feather("Data/index2word_10numerics.feather"))
tf_dist <- read_feather(
  "Data/functions_simple_10_numerics_Distribution_indiv_scale_wrecalc_allBasins_no_extremes.feather"
  )
# scale distribution values
min_max_scale <- function(x) round((x-min(x))/(max(x)-min(x)), 11)
dist_scaled <- apply(tf_dist[,-c(1, 2, 12)], 2, min_max_scale)
# Train/val split
#train_ind <- sample(nrow(generator_data), floor(0.8*nrow(generator_data)))
#saveRDS(train_ind, "train_ind_no_extremes.rds")
train_ind <- readRDS("train_ind_no_extremes.rds")
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
# function from index to words to sentences
ind2word <- function(x, index2word){
  x <- x[x != 0]
  tf <- index2word[x]
  return(paste(tf, collapse = ""))
}

# Model variables ------------------------------------------------------------------------
# variables for model architecture
max_sent_length <- ncol(generator_data)
embedding_length <- 5
emb_input_dim <- length(index2word) + 1 
distibution_dim <- ncol(train_dist)
input_cols <- embedding_length
input_rows <- max_sent_length
intermediate_dim <- 156 
latent_dim <- 6 
decoder_dense_dim <- 20
epsilon_std <- 1.0
# variables for training
epochs <- 1 # custom loop
batch_size <- 1000
kl_weight <- 100
dist_weight <- 1000

# CCN encoder ----------------------------------------------------------------------------
x <- layer_input(shape = c(max_sent_length, NULL), name = "tf_input")
# Embedding
embedding <- x %>% layer_embedding(input_dim = emb_input_dim,
                                   output_dim = embedding_length,
                                   input_length = max_sent_length,
                                   trainable = TRUE) %>%
  layer_reshape(target_shape = c(input_rows, input_cols, 1))
# 1. column 3x3 Filter
encoder_conv1 <- layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation='tanh')
encoder_pool1 <- layer_max_pooling_2d(pool_size = c(1, input_cols - 2))
column1 <- embedding %>%
  encoder_conv1 %>%
  encoder_pool1 %>%
  layer_flatten
# 2. column 4x4 filter
encoder_conv2 <- layer_conv_2d(filters = 16, kernel_size = c(4, 4), activation='tanh')
encoder_pool2 <- layer_max_pooling_2d(pool_size = c(1, input_cols - 3))
column2 <- embedding %>%
  encoder_conv2 %>%
  encoder_pool2 %>%
  layer_flatten
# 3. column 5x5 filter
encoder_conv3 <- layer_conv_2d(filters = 16, kernel_size = c(5, 5), activation='tanh')
encoder_pool3 <- layer_max_pooling_2d(pool_size = c(1, input_cols - 4))
column3 <- embedding %>%
  encoder_conv3 %>%
  encoder_pool3 %>%
  layer_flatten
# Latent space
aggregated_space <- layer_concatenate(list(column1, column2, column3)) %>%
  layer_dense(units = 453, activation = "tanh")
cnn_encoding <- aggregated_space %>%
  layer_dense(units = intermediate_dim, activation = "relu")

# Distribution properties encoding -------------------------------------------------------
d <- layer_input(shape = c(distibution_dim), name = "dist_input")
dist_encoder <- d %>% layer_dense(units = 80, activation = "selu") %>%
  layer_dense(units = 40, activation = "selu") %>%
  layer_dense(units = 20, activation = "selu") 
cnn_encoding <-  layer_concatenate(list(cnn_encoding, dist_encoder))

# VAE sampling encoded space -------------------------------------------------------------
z_mean <-  cnn_encoding %>% layer_dense(units = latent_dim, activation = "linear")
z_log_var <-  cnn_encoding %>% layer_dense(units = latent_dim, activation = "linear")

# Sampling function
sampling <- function(args){
  z_mean <- args[[1]]
  z_log_var <- args[[2]]
  epsilon <- k_random_normal(
    shape = latent_dim,
    mean = 0.,
    stddev = epsilon_std
  )
  z_mean + k_exp(z_log_var/2)*epsilon
}
# Latent space
z <- layer_lambda(list(z_mean, z_log_var), sampling)

# LSTM decoding TFs ----------------------------------------------------------------------
for(i in 1:input_rows) {
  assign(paste0("decoder0", i), layer_dense(units = latent_dim, activation = "selu"))
}
for(i in 1:input_rows) {
  decoder0 <- get(paste0("decoder0", i))
  tmp_layer <- z %>% decoder0 %>% layer_reshape(target_shape = c(NULL, 1, latent_dim))
  assign(paste0("decoder0Z", i), tmp_layer)
}
decoder2 <- bidirectional(layer = layer_lstm(units = 200, return_sequences = TRUE), 
                          merge_mode="sum")
decoder3 <- layer_dense(units = emb_input_dim, activation = "softmax", name = "TF_output")
decoder_list <- vector(mode = "list")
for(i in 1:input_rows) decoder_list[[i]] <- get(paste0("decoder0Z", i))
decoded <- layer_concatenate(decoder_list, axis = 1)  %>%
  decoder2 %>% decoder3

# Decoding distribution properties -------------------------------------------------------
dist_decoder1 <- layer_dense(units = 80, activation = "selu")
dist_decoder2 <- layer_dense(units = 40, activation = "selu")
dist_decoder3 <- layer_dense(units = 20, activation = "selu")
dist_decoder4 <-layer_dense(units = distibution_dim, name = "dist_output")
dist_decoded <- z %>% dist_decoder1 %>%
  dist_decoder2 %>%
  dist_decoder3 %>%
  dist_decoder4

# Full model -----------------------------------------------------------------------------
autoencoder <- keras_model(c(x, d), c(decoded, dist_decoded))
summary(autoencoder)

# Loss function --------------------------------------------------------------------------
vae_loss <- function(y_true, y_pred){
  xent_loss <- k_sum(k_sparse_categorical_crossentropy(y_true, y_pred))
  kl_loss <- -0.5*k_sum(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), axis = -1L) 
  xent_loss <- k_mean(xent_loss)
  kl_loss <- k_mean(kl_loss)
  return(k_mean(xent_loss + kl_weight * kl_loss))
}
dist_MSE <- function(y_true, y_pred){
  k_mean(k_square(y_pred - y_true), axis=-1) * dist_weight
}

# Model definition -----------------------------------------------------------------------
autoencoder %>% compile(optimizer = "adam",
                        loss = list("TF_output" = vae_loss,
                                    "dist_output" = dist_MSE))

# Train Autoencoder ----------------------------------------------------------------------
if (!dir.exists("training")){
  dir.create("training")
}
for(epoch in 1:200){
  history <- autoencoder %>% fit(x = list("tf_input" = x_train_tf,
                                          "dist_input" = train_dist[, 1:distibution_dim]),
                                 y = list("TF_output" = y_train_tf,
                                          "dist_output" = train_dist[, 1:distibution_dim]),
                                 epochs = 1,
                                 batch_size = batch_size,
                                 validation_data = list(list("tf_input" = x_val_tf,
                                                             "dist_input" = val_dist[, 1:distibution_dim]),
                                                        list("TF_output" = y_val_tf,
                                                             "dist_output" = val_dist[, 1:distibution_dim])))

  val_loss <- as.character(round(autoencoder$history$history$val_loss[[1]], 0))
  autoencoder %>% save_model_weights_hdf5(
    paste0("training/FSO_VAE-ep", epoch,
           "-bs",  batch_size,
           "loss_", val_loss, ".h5"))
}

plot(autoencoder$history, xlab = "epochs", main = "Model loss", 
     ylim = c(0, max(history$metrics$loss)),
     ylab = "model loss", type="l", col="blue")
plot(history$metrics$loss, xlab = "epochs", main = "Model loss", 
     ylim = c(0, max(history$metrics$loss)),
     ylab = "model loss", type="l", col="blue")
lines(history$metrics$val_loss, col = "darkgreen")
legend("topright", c("training","validation"), 
       col=c("blue", "darkgreen"), lty=c(1,1), bty = "n")
