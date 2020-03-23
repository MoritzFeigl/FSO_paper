#
# Function Space Variational Autoencoder
#

library(feather)
library(keras)
library(tensorflow)
cat("Load functions for function space optimization\n-------------------------------------------------\n")
cat("1. Load autoencoder, generator_tf, generator_dist and encoder model.\n")
# Model definition -----------------------------------------------------------------------
# variables for model architecture
max_sent_length <- 32
embedding_length <- 5
emb_input_dim <- 35
distibution_dim <- 9 
input_cols <- embedding_length
input_rows <- max_sent_length
intermediate_dim <- 156 
latent_dim <- 6 
decoder_dense_dim <- 20
epsilon_std <- 1.0
# variables for training
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
for(i in 1:input_rows) assign(paste0("decoder0", i), layer_dense(units = latent_dim, activation = "selu"))
for(i in 1:input_rows) {
  decoder0 <- get(paste0("decoder0", i))
  tmp_layer <- z %>% decoder0 %>% layer_reshape(target_shape = c(NULL, 1, latent_dim))
  assign(paste0("decoder0Z", i), tmp_layer)
}

decoder2 <- bidirectional(layer = layer_lstm(units = 200, return_sequences = TRUE), merge_mode="sum")
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
# compile
autoencoder %>% compile(optimizer="adam",
                        loss = list("TF_output" = vae_loss,
                                    "dist_output" = dist_MSE))

# trained model weights ------------------------------------------------------------------
autoencoder %>% load_model_weights_hdf5(paste0(path, 
                                               "models/current_model/vae_current_ls6.h5"))

# encoder decoder def --------------------------------------------------------------------
encoder <- keras_model(c(x, d), z)
# Generator (decoder)
decoder_input <- layer_input(shape = latent_dim)
for(i in 1:input_rows) {
  decoder0 <- get(paste0("decoder0", i))
  tmp_layer <- decoder_input %>% decoder0 %>% layer_reshape(target_shape = c(NULL, 1, latent_dim))
  assign(paste0("decoder0Z", i), tmp_layer)
}
decoder_list <- vector(mode = "list")
for(i in 1:input_rows) decoder_list[[i]] <- get(paste0("decoder0Z", i))
decoded_g <-  layer_concatenate(decoder_list, axis = 1)  %>%
  decoder2 %>% decoder3
#decoded_g <- decoder_input %>% decoder1 %>% decoder2 %>% time_distributed(layer = decoder3, name = "TF_output")
dist_decoded_g <- decoder_input %>% dist_decoder1 %>%
  dist_decoder2 %>%
  dist_decoder3 %>%
  dist_decoder4
generator_tf <- keras_model(decoder_input, decoded_g)
generator_dist <- keras_model(decoder_input, dist_decoded_g)

# index2word and functions ---------------------------------------------------------------
cat("2. Load index2word vector for function generation from indices.\n")
index2word <- as.character(
  read_feather(paste0(path,"Data/index2word_10numerics.feather")))
cat("3. Load tf_prediction function for predicting indices from softmax.\n")
# index to word translation
ind2word <- function(x, index2word){
  x <- x[x != 0]
  tf <- index2word[x]
  return(paste(tf, collapse = ""))
}
index_reconstructor <- function(pred_matrix){
  index <- integer(32)
  for(i in 1:nrow(pred_matrix)) index[i] <- which.max(pred_matrix[i, ])
  return(index)
}
# function for sampling from softmax
index_sampler <- function(pred_matrix){
  index <- integer(nrow(pred_matrix))
  for(i in 1:nrow(pred_matrix)) index[i] <- sample(ncol(pred_matrix), 1, prob = pred_matrix[i, ])
  return(index)
}
# helper function that makes the function evalutaion
tf_evaluation <- function(predicted_tf){
  predicted_tf_num <- predicted_tf
  for(i in c("slope", "evi", "sand", "clay", "elevation", "hand", "noise", "bdim")){
    predicted_tf_num <- gsub(i, "1", predicted_tf_num)
  }
  # if there is by chance a 11 or 10. or so from setting the tfs to 1 -> error
  if(length(grep("11", predicted_tf_num)) > 0) predicted_tf_num <- "ERROR"#gsub("11", "ERROR", predicted_tf_num)
  if(length(grep("10.", predicted_tf_num)) > 0) predicted_tf_num <- "ERROR"#gsub("10.", "ERROR", predicted_tf_num)
  if(length(grep("12", predicted_tf_num)) > 0) predicted_tf_num <- "ERROR"#gsub("12", "ERROR", predicted_tf_num)
  if(length(grep("13", predicted_tf_num)) > 0) predicted_tf_num <- "ERROR" #gsub("13", "ERROR", predicted_tf_num)
  for(i in 1:9) {
    if(length(grep(paste0(i, "1"), predicted_tf_num)) > 0) predicted_tf_num <- gsub(paste0(i, "1"), "ERROR", predicted_tf_num)
  }
  tf_eval <- try({
    eval(parse(text = paste('f_test <- function() {' ,  predicted_tf_num , '}', sep='')))
    f_test()
  }, silent = TRUE)
  return(tf_eval)
}
# main function for index prediction
tf_prediction <- function(index_pred){
  # get indices for softmax
  # evaluate resulting function
  # if error -> sample new function until no error
  index_prediction <- index_reconstructor(index_pred)
  index_prediction <- index_prediction - 1
  predicted_tf <- ind2word(index_prediction, index2word = index2word)
  tf_eval <- tf_evaluation(predicted_tf)
  fail_count <- 0
  # if we have a try error sample until 100 valid function are sampled and take the one
  # with the highest probability
  sample_df <- character()
  if(class(tf_eval) == "try-error" | is.null(tf_eval)){

    while(length(sample_df) <= 200 & fail_count < 2000){#
      index_prediction <- index_sampler(index_pred)
      index_prediction <- index_prediction - 1
      predicted_tf <- ind2word(index_prediction, index2word = index2word)
      tf_eval <- tf_evaluation(predicted_tf)
      if(class(tf_eval) == "try-error" | is.null(tf_eval)){
        fail_count <- fail_count + 1
        next
      }
      sample_df <- c(sample_df, predicted_tf)
    }
    max_id <- which.max(table(sample_df))
    predictedf_tf <- names(table(sample_df)[max_id])
  }
  return(predicted_tf)
}
tf_generator <- function(point){
  index_prob_prediction <- predict(generator_tf, point, batch_size = 1)
  point_tf <- apply(index_prob_prediction, 1, tf_prediction)
  for(i in 1:5) point_tf <- gsub("--", "-", point_tf, fixed = TRUE)
  for(i in 1:5) point_tf <- gsub("++", "-", point_tf, fixed = TRUE)
  return(point_tf)
}
cat("-------------------------------------------------\nFinished!\n")
