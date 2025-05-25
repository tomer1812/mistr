library(grf)

n_patients = 5000
n_test = 5000
p = 5
n_htes = 200

datatype <- "type8"


base_dir = "example_data/"
train_filename = ""
test_filename = ""

seed = 1
test_seed = 10001

train_df <- read.csv(paste0(base_dir, "output/", datatype, "/n_", n_patients, "_seed_", seed, "_imputed.csv"))
test_df <- read.csv(paste0(base_dir, "data/", datatype, "/n_", n_test, "_seed_", test_seed, "_test.csv"),
                    row.names = 1)

Y.max <- train_df[1, "Y.max"]

for (col in 1:n_htes) {
  ti_col_name <- paste0("Ti", col-1)
  di_col_name <- paste0("Di", col-1)

  adms_censoring_index <- which(train_df[[ti_col_name]] >= Y.max)
  train_df[adms_censoring_index, ti_col_name] <- Y.max
  train_df[adms_censoring_index, di_col_name] <- 1
}

X.train <- train_df[, paste0("X.", 1:p)]
W.train <- train_df$W

X.test <- test_df[, paste0("X.", 1:p)]

for (i in 1:n_htes) {
  if (i %% 50 == 0) {
    print(paste("Iteration", i))
  }

  if (sum(train_df[[paste0("Di", i-1)]]) != nrow(train_df)) {
    print(paste("Check column", col, "got", sum(train_df[[paste0("Di", i-1)]]), "ones instead of", nrow(train_df)))
  }

  Y.train <- train_df[[paste0("Ti", i-1)]]
  c.forest <- causal_forest(X.train, Y.train, W.train, ci.group.size = 8, min.node.size = 25)
  c.pred <- predict(c.forest, X.test, estimate.variance = TRUE)
  test_df[[paste0("HTE_", i-1)]] <- c.pred$predictions
  test_df[[paste0("VAR_HTE_", i-1)]] <- c.pred$variance.estimates
}

test_df[["MEAN_HTE"]] <- rowMeans(test_df[, paste0("HTE_", 0:(n_htes-1))])
test_df[["MEAN_VAR_HTE"]] <- rowMeans(test_df[, paste0("VAR_HTE_", 0:(n_htes-1))])

test_df[["IMPUTATION_VAR_HTE"]] <- apply(test_df[, paste0("HTE_", 0:(n_htes-1))], 1, function(x) {
  mean_x <- mean(x)  
  sum((x - mean_x)^2) / (length(x) - 1)  
})

test_df$TOTAL_VAR_MISTR <- test_df$MEAN_VAR_HTE + (1 + 1 / n_htes) * test_df$IMPUTATION_VAR_HTE

write.csv(test_df, paste0(base_dir, "output/", datatype, "/n_", n_patients, "_seed_", seed, "_test_with_htes.csv"))
