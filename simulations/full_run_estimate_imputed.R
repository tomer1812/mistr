library(grf)

n_patients = 5000
n_test = 5000
p = 5
max_drop_covariates = 20  # (p+1):(p+max_drop_covariates) will be dropped if exists
reps = 100
n_htes = 200

n_percs = 20
datatypes <- c("type1")

base_dir = "examples_data/"


for (seed in 1:(reps)) {
  for (datatype in datatypes) {
    print(paste0("Starting seed: ", seed, " datatype: ", datatype))
    start_time <- Sys.time()

    full_test_df <- data.frame()

    for (test_seed in 0:n_percs) {
      test_seed <- 5000 + test_seed
      test_df <- read.csv(paste0(base_dir, "data/", datatype, "/n_", n_test, "_seed_", test_seed, "_test.csv"),
                          row.names = 1)

      for (jj in (p+1):(p+max_drop_covariates)) {
        if (paste0('X.', jj) %in% names(test_df)) {
          test_df[[paste0('X.', jj)]] <- NULL
        }
      }

      test_df$test_seed <- test_seed
      full_test_df <- rbind(full_test_df, test_df[1 ,])
    }

    test_seed <- 10000 + seed
    test_df <- read.csv(paste0(base_dir, "data/", datatype, "/n_", n_test, "_seed_", test_seed, "_test.csv"),
                        row.names = 1)

    for (jj in (p+1):(p+max_drop_covariates)) {
      if (paste0('X.', jj) %in% names(test_df)) {
        test_df[[paste0('X.', jj)]] <- NULL
      }
    }

    test_df$test_seed <- test_seed
    full_test_df <- rbind(full_test_df, test_df)

    train_df <- read.csv(paste0(base_dir, "output/", datatype, "/n_", n_patients, "_seed_", seed, "_imputed.csv"))

    for (jj in (p+1):(p+max_drop_covariates)) {
      if (paste0('X.', jj) %in% names(train_df)) {
        train_df[[paste0('X.', jj)]] <- NULL
      }
    }

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
    W.train <- train_df$W

    X.test <- full_test_df[, paste0("X.", 1:p)]

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
      full_test_df[[paste0("HTE_", i-1)]] <- c.pred$predictions
      full_test_df[[paste0("VAR_HTE_", i-1)]] <- c.pred$variance.estimates
    }

    full_test_df[["MEAN_HTE"]] <- rowMeans(full_test_df[, paste0("HTE_", 0:(n_htes-1))])
    full_test_df[["MEAN_VAR_HTE"]] <- rowMeans(full_test_df[, paste0("VAR_HTE_", 0:(n_htes-1))])

    full_test_df[["IMPUTATION_VAR_HTE"]] <- apply(full_test_df[, paste0("HTE_", 0:(n_htes-1))], 1, function(x) {
      mean_x <- mean(x)  # Calculate the mean of the row if not using MEAN_X column
      sum((x - mean_x)^2) / (length(x) - 1)  # Using the MEAN_X column directly
    })

    full_test_df$TOTAL_VAR_OURS <- full_test_df$MEAN_VAR_HTE + (1 + 1 / n_htes) * full_test_df$IMPUTATION_VAR_HTE

    # Find indices where Y column value is greater than or equal to Y.max
    adms_censoring_index <- which(train_df[['Y']] >= Y.max)

    # For those indices, set the Y column value to Y.max
    train_df[adms_censoring_index, 'Y'] <- Y.max

    # For those indices, set the D column value to 1
    train_df[adms_censoring_index, 'D'] <- 1

    Y.train <- train_df$Y
    D.train <- train_df$D


    cs.forest = causal_survival_forest(X.train, Y.train, W.train, D.train, horizon = Y.max,
                                       num.trees = 2000, ci.group.size = 8, min.node.size = 25)
    cs.pred <- predict(cs.forest, X.test, estimate.variance = TRUE)
    full_test_df[["CSF"]] <- cs.pred$predictions
    full_test_df[["TOTAL_VAR_CSF"]] <- cs.pred$variance.estimates

    write.csv(full_test_df, paste0(base_dir, "output/", datatype, "/n_", n_patients, "_seed_", seed, "_test_with_htes.csv"))

    end_time <- Sys.time()
    print(paste("Iteration", seed, "datatype", datatype, "took", as.integer(as.numeric(end_time) - as.numeric(start_time)), "seconds"))
  }
}



n_patients = 5000
n_test = 5000
n_causal_trees = 2000
p = 36
reps = 5
n_htes = 200
datatypes <- c('type140', 'type141')


for (seed in 1:reps) {
  for (datatype in datatypes) {
    print(paste0("Starting seed: ", seed, " datatype: ", datatype))
    start_time <- Sys.time()

    test_seed <- 10000 + seed
    full_test_df <- read.csv(paste0(base_dir, "data/", datatype, "/n_", n_test, "_seed_", test_seed, "_test.csv"),
                        row.names = 1)
    full_test_df$test_seed <- test_seed

    train_df <- read.csv(paste0(base_dir, "output/", datatype, "/n_", n_patients, "_seed_", seed, "_imputed.csv"))

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
    W.train <- train_df$W

    X.test <- full_test_df[, 1:p]

    for (i in 1:n_htes) {
      if (i %% 50 == 0) {
        print(paste("Iteration", i))
      }

      if (sum(train_df[[paste0("Di", i-1)]]) != nrow(train_df)) {
        print(paste("Check column", col, "got", sum(train_df[[paste0("Di", i-1)]]), "ones instead of", nrow(train_df)))
      }

      Y.train <- train_df[[paste0("Ti", i-1)]]
      c.forest <- causal_forest(X.train, Y.train, W.train, ci.group.size = 8, min.node.size = 25, num.trees = n_causal_trees)
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

    full_test_df$TOTAL_VAR_OURS <- full_test_df$MEAN_VAR_HTE + (1 + 1 / n_htes) * full_test_df$IMPUTATION_VAR_HTE

    # Find indices where Y column value is greater than or equal to Y.max
    adms_censoring_index <- which(train_df[['Y']] >= Y.max)

    # For those indices, set the Y column value to Y.max
    train_df[adms_censoring_index, 'Y'] <- Y.max

    # For those indices, set the D column value to 1
    train_df[adms_censoring_index, 'D'] <- 1

    Y.train <- train_df$Y
    D.train <- train_df$D

    cs.forest = causal_survival_forest(X.train, Y.train, W.train, D.train, horizon = Y.max,
                                       num.trees = n_causal_trees, ci.group.size = 8, min.node.size = 25)
    cs.pred <- predict(cs.forest, X.test, estimate.variance = TRUE)
    full_test_df[["CSF"]] <- cs.pred$predictions
    full_test_df[["TOTAL_VAR_CSF"]] <- cs.pred$variance.estimates

    write.csv(full_test_df, paste0(base_dir, "output/", datatype, "/n_", n_patients, "_seed_", seed, "_test_with_htes.csv"))

    end_time <- Sys.time()
    print(paste("Iteration", seed, "datatype", datatype, "took", as.integer(as.numeric(end_time) - as.numeric(start_time)), "seconds"))
  }
}










##################### Instrumental Variable


n_patients = 5000
n_test = 5000
p = 3
start_reps = 1
reps = 100
n_htes = 200
datatypes <- c("type212")
z_loc <- 5
z_thresh <- 0.5


for (seed in start_reps:(reps)) {
  for (datatype in datatypes) {
    print(paste0("Starting seed: ", seed, " datatype: ", datatype))
    start_time <- Sys.time()

    full_test_df <- data.frame()

    for (test_seed in 0:n_percs) {
      test_seed <- 5000 + test_seed
      test_df <- read.csv(paste0(base_dir, "data/", datatype, "/n_", n_test, "_seed_", test_seed, "_test.csv"),
                          row.names = 1)

      for (jj in (p+1):(p+max_drop_covariates)) {
        if (paste0('X.', jj) %in% names(test_df)) {
          test_df[[paste0('X.', jj)]] <- NULL
        }
      }

      test_df$test_seed <- test_seed
      full_test_df <- rbind(full_test_df, test_df[1 ,])
    }

    test_seed <- 10000 + seed
    test_df <- read.csv(paste0(base_dir, "data/", datatype, "/n_", n_test, "_seed_", test_seed, "_test.csv"),
                        row.names = 1)

    for (jj in (p+1):(p+max_drop_covariates)) {
      if (paste0('X.', jj) %in% names(test_df)) {
        test_df[[paste0('X.', jj)]] <- NULL
      }
    }

    test_df$test_seed <- test_seed
    full_test_df <- rbind(full_test_df, test_df)

    train_df <- read.csv(paste0(base_dir, "output/", datatype, "/n_", n_patients, "_seed_", seed, "_imputed.csv"))

    Z.train <- (train_df[, paste0("X.", z_loc)] > z_thresh)

    for (jj in (p+1):(p+max_drop_covariates)) {
      if (paste0('X.', jj) %in% names(train_df)) {
        train_df[[paste0('X.', jj)]] <- NULL
      }
    }

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
    if (p == 1) {
      X.train <- matrix(X.train, nrow = length(X.train), ncol = 1)
    }

    W.train <- as.numeric(as.logical(train_df$W))

    D.train <- as.numeric(as.logical(train_df$D))
    Y.train <- train_df[, "Y"]

    x_columns <- paste0("X.", 1:p)

    # P(D=1|X,W) model
    formula <- as.formula(paste("D ~", paste(c(x_columns, "W"), collapse = " + ")))
    ipcw_model <- glm(formula, family = binomial, data=train_df[, c(x_columns, "W", "D")])
    ipcw_scores <- predict(ipcw_model, newdata = train_df[, c(x_columns, "W", "D")], type = "response")
    ipcw_weights <- ifelse(D.train == 1, 1 / ipcw_scores, 0)

    X_non_censored <- X.train[D.train == 1, ]
    Z_non_censored <- Z.train[D.train == 1]
    W_non_censored <- W.train[D.train == 1]
    Y_non_censored <- Y.train[D.train == 1]
    ipcw_non_censored <- ipcw_weights[D.train == 1]


    # Calculate statistics
    min_value <- min(ipcw_non_censored)
    max_value <- max(ipcw_non_censored)
    mean_value <- mean(ipcw_non_censored)
    std_dev <- sd(ipcw_non_censored)

    # Print results
    cat("Minimum:", min_value, "\n")
    cat("Maximum:", max_value, "\n")
    cat("Mean:", mean_value, "\n")
    cat("Standard Deviation:", std_dev, "\n")


    ipcw_instrumental_forest <- instrumental_forest(
      X = X_non_censored,
      W = W_non_censored,
      Y = Y_non_censored,
      Z = Z_non_censored,
      sample.weights = ipcw_non_censored
    )




    X.test <- full_test_df[, paste0("X.", 1:p)]
    if (p == 1) {
      X.test <- matrix(X.test, nrow = length(X.test), ncol = 1)
    }


    c.pred <- predict(ipcw_instrumental_forest, X.test, estimate.variance = TRUE)
    full_test_df[["HTE_IPCW_IV"]] <- c.pred$predictions
    full_test_df[["VAR_HTE_IPCW_IV"]] <- c.pred$variance.estimates



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
      mean_x <- mean(x)  # Calculate the mean of the row if not using MEAN_X column
      sum((x - mean_x)^2) / (length(x) - 1)  # Using the MEAN_X column directly
    })

    full_test_df$TOTAL_VAR_OURS <- full_test_df$MEAN_VAR_HTE + (1 + 1 / n_htes) * full_test_df$IMPUTATION_VAR_HTE


    write.csv(full_test_df, paste0(base_dir, "output/", datatype, "/n_", n_patients, "_seed_", seed, "_test_with_htes_iv.csv"))

    end_time <- Sys.time()
    print(paste("Iteration", seed, "datatype", datatype, "took", as.integer(as.numeric(end_time) - as.numeric(start_time)), "seconds"))
  }
}


for (seed in start_reps:(reps)) {
  for (datatype in datatypes) {
    print(paste0("Starting seed: ", seed, " datatype: ", datatype))
    start_time <- Sys.time()

    full_test_df <- data.frame()

    for (test_seed in 0:n_percs) {
      test_seed <- 5000 + test_seed
      test_df <- read.csv(paste0(base_dir, "data/", datatype, "/n_", n_test, "_seed_", test_seed, "_test.csv"),
                          row.names = 1)

      for (jj in (p+1):(p+max_drop_covariates)) {
        if (paste0('X.', jj) %in% names(test_df)) {
          test_df[[paste0('X.', jj)]] <- NULL
        }
      }

      test_df$test_seed <- test_seed
      full_test_df <- rbind(full_test_df, test_df[1 ,])
    }

    test_seed <- 10000 + seed
    test_df <- read.csv(paste0(base_dir, "data/", datatype, "/n_", n_test, "_seed_", test_seed, "_test.csv"),
                        row.names = 1)

    for (jj in (p+1):(p+max_drop_covariates)) {
      if (paste0('X.', jj) %in% names(test_df)) {
        test_df[[paste0('X.', jj)]] <- NULL
      }
    }

    test_df$test_seed <- test_seed
    full_test_df <- rbind(full_test_df, test_df)

    train_df <- read.csv(paste0(base_dir, "output/", datatype, "/n_", n_patients, "_seed_", seed, "_imputed.csv"))

    for (jj in (p+1):(p+max_drop_covariates)) {
      if (paste0('X.', jj) %in% names(train_df)) {
        train_df[[paste0('X.', jj)]] <- NULL
      }
    }

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
    W.train <- train_df$W

    X.test <- full_test_df[, paste0("X.", 1:p)]

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
      full_test_df[[paste0("HTE_", i-1)]] <- c.pred$predictions
      full_test_df[[paste0("VAR_HTE_", i-1)]] <- c.pred$variance.estimates
    }

    full_test_df[["MEAN_HTE"]] <- rowMeans(full_test_df[, paste0("HTE_", 0:(n_htes-1))])
    full_test_df[["MEAN_VAR_HTE"]] <- rowMeans(full_test_df[, paste0("VAR_HTE_", 0:(n_htes-1))])

    full_test_df[["IMPUTATION_VAR_HTE"]] <- apply(full_test_df[, paste0("HTE_", 0:(n_htes-1))], 1, function(x) {
      mean_x <- mean(x)  # Calculate the mean of the row if not using MEAN_X column
      sum((x - mean_x)^2) / (length(x) - 1)  # Using the MEAN_X column directly
    })

    full_test_df$TOTAL_VAR_OURS <- full_test_df$MEAN_VAR_HTE + (1 + 1 / n_htes) * full_test_df$IMPUTATION_VAR_HTE

    # Find indices where Y column value is greater than or equal to Y.max
    adms_censoring_index <- which(train_df[['Y']] >= Y.max)

    # For those indices, set the Y column value to Y.max
    train_df[adms_censoring_index, 'Y'] <- Y.max

    # For those indices, set the D column value to 1
    train_df[adms_censoring_index, 'D'] <- 1

    Y.train <- train_df$Y
    D.train <- train_df$D

    cs.forest = causal_survival_forest(X.train, Y.train, W.train, D.train, horizon = Y.max,
                                       num.trees = 2000, ci.group.size = 8, min.node.size = 25)
    cs.pred <- predict(cs.forest, X.test, estimate.variance = TRUE)
    full_test_df[["CSF"]] <- cs.pred$predictions
    full_test_df[["TOTAL_VAR_CSF"]] <- cs.pred$variance.estimates

    write.csv(full_test_df, paste0(base_dir, "output/", datatype, "/n_", n_patients, "_seed_", seed, "_test_with_htes_naive.csv"))

    end_time <- Sys.time()
    print(paste("Iteration", seed, "datatype", datatype, "took", as.integer(as.numeric(end_time) - as.numeric(start_time)), "seconds"))
  }
}

