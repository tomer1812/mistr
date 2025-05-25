library(grf)

base_dir = "./example_data/"

n_patients = 5000
n_test = 5000
p = 3
n_htes = 200
datatype <- "type200"
z_loc <- 5
z_thresh <- 0.5

seed <- 1
test_seed <- 10001

train_df <- read.csv(paste0(base_dir, "output/", datatype, "/n_", n_patients, "_seed_", seed, "_imputed.csv"))
test_df <- read.csv(paste0(base_dir, "data/", datatype, "/n_", n_test, "_seed_", test_seed, "_test.csv"),
                    row.names = 1)

Z.train <- (train_df[, paste0("X.", z_loc)] > z_thresh)

Y.max <- train_df[1, "Y.max"]

for (col in 1:n_htes) {
  if (col %% 50 == 0) {
    print(paste0("Updating Ti", col-1, " and Di", col-1))
  }

  ti_col_name <- paste0("Ti", col-1)
  di_col_name <- paste0("Di", col-1)

  # Find indices where Ti column value is greater than or equal to Y.max
  adms_censoring_index <- which(train_df[[ti_col_name]] >= Y.max)

  # For those indices, set the Ti column value to Y.max
  train_df[adms_censoring_index, ti_col_name] <- Y.max

  # For those indices, set the Di column value to 1
  train_df[adms_censoring_index, di_col_name] <- 1
}

X.train <- train_df[, paste0("X.", 1:p)]
W.train <- as.numeric(as.logical(train_df$W))
D.train <- as.numeric(as.logical(train_df$D))
Y.train <- train_df[, "Y"]

x_columns <- paste0("X.", 1:p)


for (i in 1:n_htes) {
  if (i %% 50 == 0) {
    print(paste("Iteration", i))
  }

  if (sum(train_df[[paste0("Di", i-1)]]) != nrow(train_df)) {
    print(paste("Check column", col, "got", sum(train_df[[paste0("Di", i-1)]]), "ones instead of", nrow(train_df)))
  }

  Y.train <- train_df[[paste0("Ti", i-1)]]
  c.forest <- instrumental_forest(X.train, Y.train, W.train, Z.train, ci.group.size = 8, min.node.size = 25)
  c.pred <- predict(c.forest, X.test, estimate.variance = TRUE)
  full_test_df[[paste0("HTE_", i-1)]] <- c.pred$predictions
  full_test_df[[paste0("VAR_HTE_", i-1)]] <- c.pred$variance.estimates
}

full_test_df[["MEAN_HTE"]] <- rowMeans(full_test_df[, paste0("HTE_", 0:(n_htes-1))])
full_test_df[["MEAN_VAR_HTE"]] <- rowMeans(full_test_df[, paste0("VAR_HTE_", 0:(n_htes-1))])

full_test_df[["IMPUTATION_VAR_HTE"]] <- apply(full_test_df[, paste0("HTE_", 0:(n_htes-1))], 1, function(x) {
  mean_x <- mean(x)  
  sum((x - mean_x)^2) / (length(x) - 1)  
})

full_test_df$TOTAL_VAR_MISTRIV <- full_test_df$MEAN_VAR_HTE + (1 + 1 / n_htes) * full_test_df$IMPUTATION_VAR_HTE

write.csv(full_test_df, paste0(base_dir, "output/", datatype, "/n_", n_patients, "_seed_", seed, "_test_with_htes_iv.csv"))
