library(grf)

options(warn = 1)

covariates <- c("age", "wtkg", "karnof", "cd40", "cd80", "gender", "homo",
                "race", "symptom", "drugs", "hemo", "str2")

base_output_dir = 'example_data/hiv/'


reps <- 10
time_unit <- 28

n_imputations <- 2500
case <- 60
min_observed_leaf <- 3
min_node_size <- 5
ci_group_size <- 18

arms_list <- list(c(0, 1), c(0, 2), c(0, 3))

for (rep in 1:reps) {
  for (arms in arms_list) {
    seed <- rep + 100*case + arms[2]*10000
    set.seed(seed)
    dfx <- read.csv(paste0(base_output_dir, "hiv_imputed_p", case, "_mol_", min_observed_leaf, "_nimp_", n_imputations, "_arms_", arms[1], arms[2], "_rep_", rep, ".csv"))
    train_df <- na.omit(dfx)

    full_cohort_test_df <- read.csv("example_data/ACTG175.csv")

    n_htes = 200
    
    Y.max <- (840 %/% time_unit) 
    num_trees <- 2000

    duration_col = 'U'
    event_col = 'Delta'
    treatment_col = 'W'

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
      X.train <- train_df[, covariates]
      
      W.train <- train_df$W
      Y.train <- train_df[[paste0("Ti", i-1)]]
      c.forest <- causal_forest(X.train, Y.train, W.train, ci.group.size = ci_group_size, min.node.size = min_node_size, num.trees = num_trees)
      c.pred <- predict(c.forest, X.train, estimate.variance = TRUE)
      train_df[[paste0("HTE_", i-1)]] <- c.pred$predictions
      train_df[[paste0("VAR_HTE_", i-1)]] <- c.pred$variance.estimates

      
      X.test <- full_cohort_test_df[, covariates]
      c.pred <- predict(c.forest, X.test, estimate.variance = TRUE)
      full_cohort_test_df[[paste0("HTE_", i-1)]] <- c.pred$predictions
      full_cohort_test_df[[paste0("VAR_HTE_", i-1)]] <- c.pred$variance.estimates
      
    }
    
    train_df[["MEAN_HTE"]] <- rowMeans(train_df[, paste0("HTE_", 0:(n_htes-1))])
    train_df[["MEAN_VAR_HTE"]] <- rowMeans(train_df[, paste0("VAR_HTE_", 0:(n_htes-1))])
    
    train_df[["IMPUTATION_VAR_HTE"]] <- apply(train_df[, paste0("HTE_", 0:(n_htes-1))], 1, function(x) {
      mean_x <- mean(x)  
      sum((x - mean_x)^2) / (length(x) - 1)  
    })
    
    
    full_cohort_test_df[["MEAN_HTE"]] <- rowMeans(full_cohort_test_df[, paste0("HTE_", 0:(n_htes-1))])
    full_cohort_test_df[["MEAN_VAR_HTE"]] <- rowMeans(full_cohort_test_df[, paste0("VAR_HTE_", 0:(n_htes-1))])
    
    full_cohort_test_df[["IMPUTATION_VAR_HTE"]] <- apply(full_cohort_test_df[, paste0("HTE_", 0:(n_htes-1))], 1, function(x) {
      mean_x <- mean(x)  # Calculate the mean of the row if not using MEAN_X column
      sum((x - mean_x)^2) / (length(x) - 1)  # Using the MEAN_X column directly
    })
    
    train_df$TOTAL_VAR_OURS <- train_df$MEAN_VAR_HTE + (1 + 1 / n_htes) * train_df$IMPUTATION_VAR_HTE
    full_cohort_test_df$TOTAL_VAR_OURS <- full_cohort_test_df$MEAN_VAR_HTE + (1 + 1 / n_htes) * full_cohort_test_df$IMPUTATION_VAR_HTE
    
    
    # Find indices where Y column value is greater than or equal to Y.max
    adms_censoring_index <- which(train_df[['Y']] >= Y.max)
    
    # For those indices, set the Y column value to Y.max
    train_df[adms_censoring_index, 'U'] <- Y.max
    
    # For those indices, set the D column value to 1
    train_df[adms_censoring_index, 'Delta'] <- 1
    
    Y.train <- train_df$U
    D.train <- train_df$Delta
    
    cs.forest = causal_survival_forest(X.train, Y.train, W.train, D.train, horizon = Y.max,
                                       num.trees = num_trees, ci.group.size = ci_group_size, min.node.size = min_node_size) 
    cs.pred <- predict(cs.forest, X.train, estimate.variance = TRUE)
    train_df[["CSF"]] <- cs.pred$predictions
    train_df[["TOTAL_VAR_CSF"]] <- cs.pred$variance.estimates

    
    X.test <- full_cohort_test_df[, covariates]
    cs.pred <- predict(cs.forest, X.test, estimate.variance = TRUE)
    full_cohort_test_df[["CSF"]] <- cs.pred$predictions
    full_cohort_test_df[["TOTAL_VAR_CSF"]] <- cs.pred$variance.estimates
    
    end_time <- Sys.time()
    print(paste("Iteration", seed, "took", as.integer(as.numeric(end_time) - as.numeric(start_time)), "seconds"))
    
    write.csv(train_df, paste0(base_output_dir, "hiv_imputed_with_htes_p", case, "_ci_", ci_group_size, "_node_", min_node_size, "_mol_", min_observed_leaf, "_nimp_", n_imputations, "_arms_", arms[1], arms[2], "_nt_", num_trees, "_rep_", rep, ".csv"))
    write.csv(full_cohort_test_df, paste0(base_output_dir, "hiv_imputed_with_htes_full_cohort_test_p", case, "_ci_", ci_group_size, "_node_", min_node_size, "_mol_", min_observed_leaf, "_nimp_", n_imputations, "_arms_", arms[1], arms[2], "_nt_", num_trees , "_rep_", rep, ".csv"))
  }
}