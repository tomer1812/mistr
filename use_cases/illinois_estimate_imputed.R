library(grf)

options(warn = 1)

covariates <- c('age', 'male', 'white', 'black')

base_output_dir = 'example_data/illinois/'

group = 'hie'
#group = 'jsie'


print(paste('Running', group, 'group with test'))

rep <- 1
seed <- 1


age <- seq(20, 55, by = 1)        # Age values from 20 to 55
gender <- c(0, 1)                 # Gender values: 0 and 1
black <- c(0, 1)                  # Black values: 0 and 1
white <- c(0, 1)                  # White values: 0 and 1


n_imputations <- 2500
min_observed_leaf <- 18
min_node_size <- 18
ci_group_size <- 18


set.seed(seed)
dfx <- read.csv(paste0(base_output_dir, group, "_imputed_mol_", min_observed_leaf, "_nimp_", n_imputations, "_rep_", rep, ".csv"))
train_df <- na.omit(dfx)

n_htes = 200

Y.max <- (7*25)
num_trees <- 2000

duration_col = 'Y'
event_col = 'D'
treatment_col = 'W'
intrumental_col = 'Z'

print(paste0("Starting seed: ", seed))
start_time <- Sys.time()

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

for (i in 1:n_htes) {
  if (i %% 50 == 0) {
    print(paste("Iteration", i))
  }

  if (sum(train_df[[paste0("Di", i-1)]]) != nrow(train_df)) {
    print(paste("Check column", col, "got", sum(train_df[[paste0("Di", i-1)]]), "ones instead of", nrow(train_df)))
  }
  X.train <- as.matrix(train_df[, covariates])

  W.train <- train_df$W
  Z.train <- train_df$Z
  Y.train <- train_df[[paste0("Ti", i-1)]]
  c.forest <- instrumental_forest(X.train, Y.train, W.train, Z.train, ci.group.size = ci_group_size,
                                  min.node.size = min_node_size, num.trees = num_trees)
  c.pred <- predict(c.forest, X.train, estimate.variance = TRUE)
  train_df[[paste0("IV_HTE_", i-1)]] <- c.pred$predictions
  train_df[[paste0("IV_VAR_HTE_", i-1)]] <- c.pred$variance.estimates

  c.forest <- causal_forest(X.train, Y.train, W.train, ci.group.size = ci_group_size, min.node.size = min_node_size, num.trees = num_trees)
  c.pred <- predict(c.forest, X.train, estimate.variance = TRUE)
  train_df[[paste0("HTE_", i-1)]] <- c.pred$predictions
  train_df[[paste0("VAR_HTE_", i-1)]] <- c.pred$variance.estimates
}

train_df[["IV_MEAN_HTE"]] <- rowMeans(train_df[, paste0("IV_HTE_", 0:(n_htes-1))])
train_df[["IV_MEAN_VAR_HTE"]] <- rowMeans(train_df[, paste0("IV_VAR_HTE_", 0:(n_htes-1))])

train_df[["IV_IMPUTATION_VAR_HTE"]] <- apply(train_df[, paste0("IV_HTE_", 0:(n_htes-1))], 1, function(x) {
  mean_x <- mean(x)
  sum((x - mean_x)^2) / (length(x) - 1)
})

train_df$IV_TOTAL_VAR_OURS <- train_df$IV_MEAN_VAR_HTE + (1 + 1 / n_htes) * train_df$IV_IMPUTATION_VAR_HTE


train_df[["MEAN_HTE"]] <- rowMeans(train_df[, paste0("HTE_", 0:(n_htes-1))])
train_df[["MEAN_VAR_HTE"]] <- rowMeans(train_df[, paste0("VAR_HTE_", 0:(n_htes-1))])

train_df[["IMPUTATION_VAR_HTE"]] <- apply(train_df[, paste0("HTE_", 0:(n_htes-1))], 1, function(x) {
  mean_x <- mean(x)  # Calculate the mean of the row if not using MEAN_X column
  sum((x - mean_x)^2) / (length(x) - 1)  # Using the MEAN_X column directly
})

train_df$TOTAL_VAR_OURS <- train_df$MEAN_VAR_HTE + (1 + 1 / n_htes) * train_df$IMPUTATION_VAR_HTE



# Find indices where Y column value is greater than or equal to Y.max
adms_censoring_index <- which(train_df[['Y']] >= Y.max)

# For those indices, set the Y column value to Y.max
train_df[adms_censoring_index, 'Y'] <- Y.max

# For those indices, set the D column value to 1
train_df[adms_censoring_index, 'D'] <- 1

Y.train <- train_df$Y
D.train <- train_df$D

cs.forest = causal_survival_forest(X.train, Y.train, W.train, D.train, horizon = Y.max,
                                   num.trees = num_trees, ci.group.size = ci_group_size, min.node.size = min_node_size)
cs.pred <- predict(cs.forest, X.train, estimate.variance = TRUE)
train_df[["CSF"]] <- cs.pred$predictions
train_df[["TOTAL_VAR_CSF"]] <- cs.pred$variance.estimates

end_time <- Sys.time()
print(paste("Iteration", seed, "took", as.integer(as.numeric(end_time) - as.numeric(start_time)), "seconds"))

write.csv(train_df, paste0(base_output_dir, group, "_train_imputed_with_htes_", "ci_", ci_group_size, "_node_", min_node_size, "_mol_", min_observed_leaf, "_nimp_", n_imputations, "_nt_", num_trees, "_rep_", rep, ".csv"))
