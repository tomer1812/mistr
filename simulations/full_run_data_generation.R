library(rje)

start_reps = 1
reps = 100
n.mc = 20000
n_percs = 20


generate_cases <- function (n, p, Y.max = NULL, X = NULL, 
                            n.mc = 20000, dgp = c("type1", "type2", "type3", "type4", 
                                                  "type5", "type6", "type7", "type8", "type9", "type10",
                                                  'type140', 'type141','type142','type143','type144', 
                                                  "type200", "type201", "type202", "type203", "type204")) 
{
  dgp <- match.arg(dgp)
  if (dgp == "type1") {
    
    Y.max <- 0.7
    e <- (1 + dbeta(X[, 1], 2, 4))/4
    W <- rbinom(n, 1, e)
    I1 <- X[, 1] < 0.5
    ft <- exp(-1.85 - 0.8 * I1 + 0.7 * sqrt(X[, 2]) + 0.2 * 
                X[, 3] + (0.7 - 0.4 * I1 - 0.4 * sqrt(X[, 2])) * 
                W + rnorm(n))
    
    failure.time <- ft
    
    numerator <- -log(runif(n))
    denominator <- exp(-1.75 - 0.5 * sqrt(X[, 2]) + 0.2 * 
                         X[, 3] + (1.15 + 0.5 * I1 - 0.3 * sqrt(X[, 2])) * 
                         W)
    censor.time <- (numerator/denominator)^(1/2)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    
    eps <- rnorm(n.mc)
    for (i in 1:n) {
      ft0 <- exp(-1.85 - 0.8 * I1[i] + 0.7 * sqrt(X[i, 
                                                    2]) + 0.2 * X[i, 3] + eps)
      ft1 <- exp(-1.85 - 0.8 * I1[i] + 0.7 * sqrt(X[i, 
                                                    2]) + 0.2 * X[i, 3] + 0.7 - 0.4 * I1[i] - 0.4 * 
                   sqrt(X[i, 2]) + eps)
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  else if (dgp == "type2") {
    Y.max <- 0.7
    e <- (1 + dbeta(X[, 1], 2, 4))/4
    W <- rbinom(n, 1, e)
    numerator <- -log(runif(n))
    cox.ft <- (numerator/exp(X[, 1] + (-0.5 + X[, 2]) * W))^2
    failure.time <- cox.ft
    
    censor.time <- 3 * runif(n)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    
    numerator <- -log(runif(n.mc))
    for (i in 1:n) {
      cox.ft0 <- (numerator/exp(X[i, 1] + (-0.5 + X[i, 2]) * 0))^2
      cox.ft1 <- (numerator/exp(X[i, 1] + (-0.5 + X[i, 2]) * 1))^2
      cate[i] <- mean(pmin(cox.ft1, Y.max) - pmin(cox.ft0, Y.max))
    }
  }
  else if (dgp == "type3") {

    Y.max <- 11
    e <- (1 + dbeta(X[, 1], 2, 4))/4

    W <- rbinom(n, 1, e)
    lambda.failure <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3) * W
    
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 12 + log(1 + exp(X[, 3]))
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    cate.prob <- rep(NA, n)
    lambda.failure.0 <- X[, 2]^2 + X[, 3] + 6
    lambda.failure.1 <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3)
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
      cate.prob[i] <- mean(ft1 > y0) - mean(ft0 > y0)
    }
    cate.sign <- sign(sqrt(X[, 1]) - 0.3)
  }
  else if (dgp == "type4") {

    Y.max <- 3
    e <- 1/((1 + exp(-X[, 1])) * (1 + exp(-X[, 2])))
    W <- rbinom(n, 1, e)
    lambda.failure <- X[, 2] + X[, 3] + pmax(0, X[, 1] - 0.3) * W
    
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 1 + log(1 + exp(X[, 3]))
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    
    lambda.failure.0 <- X[, 2] + X[, 3]
    lambda.failure.1 <- X[, 2] + X[, 3] + pmax(0, X[, 1] - 0.3)
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  } 
  else if (dgp == "type5") {
    
    Y.max <- 6
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    I4 <- X[, 4] < 0.5
    
    lambda.failure <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3) * W
    
    failure.time <- rpois(n, lambda = lambda.failure)
    
    censor.time1 <- rep(1000000, n)
    censor.time2 <- 1 + I4 #sample(1:3, n, replace = TRUE)
    
    selection <- runif(n) < 0.6
    censor.time <- ifelse(selection, censor.time1, censor.time2)
    
    
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    
    lambda.failure.0 <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3) * 0
    lambda.failure.1 <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3) * 1
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  } 
  else if (dgp == "type6") {

    Y.max <- 6
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    lambda.failure <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3) * W
    
    failure.time <- rpois(n, lambda = lambda.failure)

    lambda.censor <- 3 + log(1 + exp(2*X[, 2] + X[, 3]))
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    
    lambda.failure.0 <- X[, 2]^2 + X[, 3] + 6
    lambda.failure.1 <- X[, 2]^2 + X[, 3] + 6 + 2 * (sqrt(X[, 1]) - 0.3)
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  else if (dgp == "type7") {
    
    Y.max <- 7
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)

    lambda.failure <- X[, 2]^2 + X[, 3] + 7 + 2 * (sqrt(X[, 1]) - 0.3) * W
    failure.time <- rpois(n, lambda = lambda.failure)

    lambda.censor <- 3 + 4*X[, 6] + 2*X[, 7] 
    censor.time <- rpois(n, lambda = lambda.censor)
    
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    
    lambda.failure.0 <- X[, 2]^2 + X[, 3] + 7 + 2 * (sqrt(X[, 1]) - 0.3) * 0
    lambda.failure.1 <- X[, 2]^2 + X[, 3] + 7 + 2 * (sqrt(X[, 1]) - 0.3) * 1
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  } 
  else if (dgp == "type8") {

    Y.max <- 6
    e <- (1 + dbeta(X[, 1], 2, 4)) / 4
    W <- rbinom(n, 1, e)
    lambda.failure <- X[, 2]^2 + X[, 3] + 7 + 2 * (sqrt(X[, 1]) - 0.3) * W
    
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 3 
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    
    lambda.failure.0 <- X[, 2]^2 + X[, 3] + 7
    lambda.failure.1 <- X[, 2]^2 + X[, 3] + 7 + 2 * (sqrt(X[, 1]) - 0.3)
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  else if (dgp == "type9") {
    
    Y.max <- 0.7
    e <- (1 + dbeta(X[, 1], 2, 4))/4
    W <- rbinom(n, 1, e)
    I1 <- X[, 1] < 0.5
    ft <- exp(0.3 - 0.5 * I1 + 0.5 * sqrt(X[, 2]) + 0.2 * 
                X[, 3] + (1 - 0.8 * I1 - 0.8 * sqrt(X[, 2])) * W + rnorm(n))
    
    failure.time <- ft
    
    numerator <- -log(runif(n))
    denominator <- exp(-0.9 + 2 * sqrt(X[, 2]) + 2 * X[, 3] + (1.15 + 
                                            0.5 * I1 - 0.3 * sqrt(X[, 2])) * W)
    censor.time <- (numerator/denominator)^(1/2)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    
    eps <- rnorm(n.mc)
    for (i in 1:n) {
      ft0 <- exp(0.3 - 0.5 * I1[i] + 0.5 * sqrt(X[i, 2]) + 0.2 * X[i, 3] + eps)

      ft1 <- exp(0.3 - 0.5 * I1[i] + 0.5 * sqrt(X[i, 2]) + 0.2 * 
                   X[i, 3] + (1 - 0.8 * I1[i] - 0.8 * sqrt(X[i, 2])) * 1 + eps)
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  else if (dgp == "type10") {

    Y.max <- 0.7
    e <- (1 + dbeta(X[, 1], 2, 4))/4
    W <- rbinom(n, 1, e)
    numerator <- -log(runif(n))
    cox.ft <- (numerator/exp(X[, 1] + (-0.5 + X[, 2]) * W))^2
    
    failure.time <- cox.ft
    
    censor.time <- 0.1 * runif(n)
    
    censor.time1 <- rep(1000000, n)
    censor.time2 <- 0.05 * runif(n)
    
    selection <- runif(n) < 0.15
    censor.time <- ifelse(selection, censor.time1, censor.time2)
    
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    cate.prob <- rep(NA, n)
    numerator <- -log(runif(n.mc))
    for (i in 1:n) {
      cox.ft0 <- (numerator/exp(X[i, 1] + (-0.5 + X[i, 2]) * 0))^2
      cox.ft1 <- (numerator/exp(X[i, 1] + (-0.5 + X[i, 2]) * 1))^2
      cate[i] <- mean(pmin(cox.ft1, Y.max) - pmin(cox.ft0, Y.max))

    }
  }
  
  else if (dgp == "type140") {
    
    # Mimic IV based covariates
    
    Y.max <- 28
    W <- rbinom(n, 1, 0.5)
    
    # 1. AnionGap	
    # 2. Bicarbonate	
    # 3. CalciumTotal	
    # 4. Chloride	
    # 5. Creatinine	
    # 6. Glucose	
    # 7. Magnesium	
    # 8. Phosphate	
    # 9. Potassium	
    # 10. Sodium	
    # 11. UreaNitrogen	
    # 12. Hematocrit	
    # 13. Hemoglobin	
    # 14. MCH	
    # 15. MCHC	
    # 16. MCV	
    # 17. PlateletCount	
    # 18. RDW	
    # 19. RedBloodCells	
    # 20. WhiteBloodCells	
    # 21. Insurance_Medicare	
    # 22. Insurance_Other	
    # 23. Marital_MARRIED	
    # 24. Marital_SINGLE	
    # 25. Marital_WIDOWED	
    # 26. Ethnicity_BLACK	
    # 27. Ethnicity_HISPANIC	
    # 28. Ethnicity_OTHER	
    # 29. Ethnicity_WHITE	
    # 30. AdmsCount_2	
    # 31. AdmsCount_3up	
    # 32. night_admission	
    # 33. gender	
    # 34. direct_emrgency_flag	
    # 35. last_less_than_diff	
    # 36. standardized_age
    
    clipped_standardized_age <- pmax(pmin(X[, 36], 2), -2)
    
    lambda.failure <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * W)
    
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 21
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)

    
    lambda.failure.0 <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * 0)
    
    lambda.failure.1 <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * 1)
    
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  } 
  else if (dgp == "type141") {
    
    # Mimic IV based covariates

    Y.max <- 28
    W <- rbinom(n, 1, 0.5)
    
    # 1. AnionGap	
    # 2. Bicarbonate	
    # 3. CalciumTotal	
    # 4. Chloride	
    # 5. Creatinine	
    # 6. Glucose	
    # 7. Magnesium	
    # 8. Phosphate	
    # 9. Potassium	
    # 10. Sodium	
    # 11. UreaNitrogen	
    # 12. Hematocrit	
    # 13. Hemoglobin	
    # 14. MCH	
    # 15. MCHC	
    # 16. MCV	
    # 17. PlateletCount	
    # 18. RDW	
    # 19. RedBloodCells	
    # 20. WhiteBloodCells	
    # 21. Insurance_Medicare	
    # 22. Insurance_Other	
    # 23. Marital_MARRIED	
    # 24. Marital_SINGLE	
    # 25. Marital_WIDOWED	
    # 26. Ethnicity_BLACK	
    # 27. Ethnicity_HISPANIC	
    # 28. Ethnicity_OTHER	
    # 29. Ethnicity_WHITE	
    # 30. AdmsCount_2	
    # 31. AdmsCount_3up	
    # 32. night_admission	
    # 33. gender	
    # 34. direct_emrgency_flag	
    # 35. last_less_than_diff	
    # 36. standardized_age
    
    clipped_standardized_age <- pmax(pmin(X[, 36], 2), -2)
    
    lambda.failure <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * W)
    
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 23 
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    
    
    
    lambda.failure.0 <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * 0)
    
    lambda.failure.1 <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * 1)
    
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  } 
  else if (dgp == "type142") {
    
    # Mimic IV based covariates

    Y.max <- 28
    W <- rbinom(n, 1, 0.5)
    
    # 1. AnionGap	
    # 2. Bicarbonate	
    # 3. CalciumTotal	
    # 4. Chloride	
    # 5. Creatinine	
    # 6. Glucose	
    # 7. Magnesium	
    # 8. Phosphate	
    # 9. Potassium	
    # 10. Sodium	
    # 11. UreaNitrogen	
    # 12. Hematocrit	
    # 13. Hemoglobin	
    # 14. MCH	
    # 15. MCHC	
    # 16. MCV	
    # 17. PlateletCount	
    # 18. RDW	
    # 19. RedBloodCells	
    # 20. WhiteBloodCells	
    # 21. Insurance_Medicare	
    # 22. Insurance_Other	
    # 23. Marital_MARRIED	
    # 24. Marital_SINGLE	
    # 25. Marital_WIDOWED	
    # 26. Ethnicity_BLACK	
    # 27. Ethnicity_HISPANIC	
    # 28. Ethnicity_OTHER	
    # 29. Ethnicity_WHITE	
    # 30. AdmsCount_2	
    # 31. AdmsCount_3up	
    # 32. night_admission	
    # 33. gender	
    # 34. direct_emrgency_flag	
    # 35. last_less_than_diff	
    # 36. standardized_age
    
    clipped_standardized_age <- pmax(pmin(X[, 36], 2), -2)
    
    lambda.failure <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * W)
    
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 24.7 
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)

    
    lambda.failure.0 <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * 0)
    
    lambda.failure.1 <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * 1)
    
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  } 
  else if (dgp == "type143") {
    
    # Mimic IV based covariates
    
    Y.max <- 28
    W <- rbinom(n, 1, 0.5)
    
    # 1. AnionGap	
    # 2. Bicarbonate	
    # 3. CalciumTotal	
    # 4. Chloride	
    # 5. Creatinine	
    # 6. Glucose	
    # 7. Magnesium	
    # 8. Phosphate	
    # 9. Potassium	
    # 10. Sodium	
    # 11. UreaNitrogen	
    # 12. Hematocrit	
    # 13. Hemoglobin	
    # 14. MCH	
    # 15. MCHC	
    # 16. MCV	
    # 17. PlateletCount	
    # 18. RDW	
    # 19. RedBloodCells	
    # 20. WhiteBloodCells	
    # 21. Insurance_Medicare	
    # 22. Insurance_Other	
    # 23. Marital_MARRIED	
    # 24. Marital_SINGLE	
    # 25. Marital_WIDOWED	
    # 26. Ethnicity_BLACK	
    # 27. Ethnicity_HISPANIC	
    # 28. Ethnicity_OTHER	
    # 29. Ethnicity_WHITE	
    # 30. AdmsCount_2	
    # 31. AdmsCount_3up	
    # 32. night_admission	
    # 33. gender	
    # 34. direct_emrgency_flag	
    # 35. last_less_than_diff	
    # 36. standardized_age
    
    clipped_standardized_age <- pmax(pmin(X[, 36], 2), -2)
    
    lambda.failure <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * W)
    
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 26.5
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    
    
    lambda.failure.0 <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * 0)
    
    lambda.failure.1 <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * 1)
    
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  } 
  else if (dgp == "type144") {
    
    # Mimic IV based covariates

    Y.max <- 28
    W <- rbinom(n, 1, 0.5)
    
    # 1. AnionGap	
    # 2. Bicarbonate	
    # 3. CalciumTotal	
    # 4. Chloride	
    # 5. Creatinine	
    # 6. Glucose	
    # 7. Magnesium	
    # 8. Phosphate	
    # 9. Potassium	
    # 10. Sodium	
    # 11. UreaNitrogen	
    # 12. Hematocrit	
    # 13. Hemoglobin	
    # 14. MCH	
    # 15. MCHC	
    # 16. MCV	
    # 17. PlateletCount	
    # 18. RDW	
    # 19. RedBloodCells	
    # 20. WhiteBloodCells	
    # 21. Insurance_Medicare	
    # 22. Insurance_Other	
    # 23. Marital_MARRIED	
    # 24. Marital_SINGLE	
    # 25. Marital_WIDOWED	
    # 26. Ethnicity_BLACK	
    # 27. Ethnicity_HISPANIC	
    # 28. Ethnicity_OTHER	
    # 29. Ethnicity_WHITE	
    # 30. AdmsCount_2	
    # 31. AdmsCount_3up	
    # 32. night_admission	
    # 33. gender	
    # 34. direct_emrgency_flag	
    # 35. last_less_than_diff	
    # 36. standardized_age
    
    clipped_standardized_age <- pmax(pmin(X[, 36], 2), -2)
    
    lambda.failure <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * W)
    
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 29
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)

    
    lambda.failure.0 <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * 0)
    
    lambda.failure.1 <- 3*(10 + 0.25 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) + 0.25*clipped_standardized_age +
      0.5*(-0.3 - 0.5 * (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5]) - 0.5*clipped_standardized_age) * 1)
    
    for (i in 1:n) {
      ft0 <- rpois(n.mc, lambda.failure.0[i])
      ft1 <- rpois(n.mc, lambda.failure.1[i])
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  } 

  else if (dgp == "type200") {
    # Instrumental variable

    Y.max <- 8
    Z <- as.numeric(X[, 5] > 0.5)
    
    W.star <- 0.5*X[, 4] + 0.5*Z  + 0.2*rnorm(n)
    W <- as.numeric(W.star > 0.5)

    lambda.failure <- 2*X[, 1] + X[, 2] + 2*X[, 4] + 4 + 2 * (sqrt(X[, 1]) - 0.3) * W

    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 7 
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    
    cate <- rep(NA, n)
    
    rand_u = matrix(runif(n * n.mc), n, n.mc)
    
    for (i in 1:n) {
      lambda.failure.0 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 4
      lambda.failure.1 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 4 + 2 * (sqrt(X[i, 1]) - 0.3)
      ft0 <- rpois(n.mc, lambda.failure.0)
      ft1 <- rpois(n.mc, lambda.failure.1)
      
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  
  else if (dgp == "type201") {
    # 90%+ censoring, discrete
    # Instrumental variable

    Y.max <- 8
    Z <- as.numeric(X[, 5] > 0.5)
    
    W.star <- 0.5*X[, 4] + 0.35*Z  + 0.2*rnorm(n)
    W <- as.numeric(W.star > 0.5)

    lambda.failure <- 2*X[, 1] + X[, 2] + 2*X[, 4] + 4 + 2 * (sqrt(X[, 1]) - 0.3) * W

    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 7 
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    
    rand_u = matrix(runif(n * n.mc), n, n.mc)
    
    for (i in 1:n) {
      lambda.failure.0 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 4
      lambda.failure.1 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 4 + 2 * (sqrt(X[i, 1]) - 0.3)
      ft0 <- rpois(n.mc, lambda.failure.0)
      ft1 <- rpois(n.mc, lambda.failure.1)
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  
  else if (dgp == "type202") {
    # 90%+ censoring, discrete
    # Instrumental variable

    Y.max <- 8
    Z <- as.numeric(X[, 5] > 0.5)
    
    W.star <- 0.5*sqrt(X[, 4]) + 0.35*Z  + 0.2*rnorm(n)
    W <- as.numeric(W.star > 0.5)

    lambda.failure <- 2*X[, 1] + X[, 2] + 3*sqrt(X[, 4]) + 3 + 2 * (sqrt(X[, 1]) - 0.3) * W
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 7 
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    
    rand_u = matrix(runif(n * n.mc), n, n.mc)
    
    for (i in 1:n) {
      lambda.failure.0 <- 2 * X[i, 1] + X[i, 2] + 3*sqrt(rand_u[i, ]) + 3
      lambda.failure.1 <- 2 * X[i, 1] + X[i, 2] + 3*sqrt(rand_u[i, ]) + 3 + 2 * (sqrt(X[i, 1]) - 0.3)
      ft0 <- rpois(n.mc, lambda.failure.0)
      ft1 <- rpois(n.mc, lambda.failure.1)
      
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  
  else if (dgp == "type203") {
    # 90%+ censoring, discrete
    # Instrumental variable
    

    Y.max <- 8
    Z <- as.numeric(X[, 5] > 0.5)
    
    W.star <- 0.5*X[, 4] + 0.5*Z  + 0.2*rnorm(n)
    W <- as.numeric(W.star > 0.5)

    lambda.failure <- 2*X[, 1] + X[, 2] + 2*X[, 4] + 5 + 2 * (sqrt(X[, 1]) - 0.3) * W
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 6 
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)
    
    rand_u = matrix(runif(n * n.mc), n, n.mc)
    
    for (i in 1:n) {
      lambda.failure.0 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 5
      lambda.failure.1 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 5 + 2 * (sqrt(X[i, 1]) - 0.3)
      ft0 <- rpois(n.mc, lambda.failure.0)
      ft1 <- rpois(n.mc, lambda.failure.1)
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  
  else if (dgp == "type204") {
    # 90%+ censoring, discrete
    # Instrumental variable

    Y.max <- 8
    Z <- as.numeric(X[, 5] > 0.5)
    
    W.star <- 0.5*X[, 4] + 0.5*Z  + 0.2*rnorm(n)
    W <- as.numeric(W.star > 0.5)

    lambda.failure <- 2*X[, 1] + X[, 2] + 2*X[, 4] + 6 + 2 * (sqrt(X[, 1]) - 0.3) * W
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 4
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)

    rand_u = matrix(runif(n * n.mc), n, n.mc)
    
    for (i in 1:n) {
      lambda.failure.0 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 6
      lambda.failure.1 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 6 + 2 * (sqrt(X[i, 1]) - 0.3)
      ft0 <- rpois(n.mc, lambda.failure.0)
      ft1 <- rpois(n.mc, lambda.failure.1)
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  
  else if (dgp == "type205") {
    # 90%+ censoring, discrete
    # Instrumental variable
    # Weak instrument analysis

    Y.max <- 8
    Z <- as.numeric(X[, 5] > 0.5)
    
    W.star <- 0.5*X[, 4] + 0.4*Z  + 0.2*rnorm(n)
    W <- as.numeric(W.star > 0.5)

    lambda.failure <- 2*X[, 1] + X[, 2] + 2*X[, 4] + 4 + 2 * (sqrt(X[, 1]) - 0.3) * W
    
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 7 
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    
    cate <- rep(NA, n)
    rand_u = matrix(runif(n * n.mc), n, n.mc)
    
    for (i in 1:n) {
      lambda.failure.0 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 4
      lambda.failure.1 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 4 + 2 * (sqrt(X[i, 1]) - 0.3)
      ft0 <- rpois(n.mc, lambda.failure.0)
      ft1 <- rpois(n.mc, lambda.failure.1)
      
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  
  else if (dgp == "type206") {
    # 90%+ censoring, discrete
    # Instrumental variable
    # Weak instrument analysis

    Y.max <- 8
    Z <- as.numeric(X[, 5] > 0.5)
    
    W.star <- 0.5*X[, 4] + 0.3*Z  + 0.2*rnorm(n)
    W <- as.numeric(W.star > 0.5)

    lambda.failure <- 2*X[, 1] + X[, 2] + 2*X[, 4] + 4 + 2 * (sqrt(X[, 1]) - 0.3) * W
    
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 7 
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    
    cate <- rep(NA, n)
    rand_u = matrix(runif(n * n.mc), n, n.mc)
    
    for (i in 1:n) {
      lambda.failure.0 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 4
      lambda.failure.1 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 4 + 2 * (sqrt(X[i, 1]) - 0.3)
      ft0 <- rpois(n.mc, lambda.failure.0)
      ft1 <- rpois(n.mc, lambda.failure.1)
      
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  
  else if (dgp == "type207") {
    # 90%+ censoring, discrete
    # Instrumental variable
    # Weak instrument analysis
  
    Y.max <- 8
    Z <- as.numeric(X[, 5] > 0.5)
    
    W.star <- 0.5*X[, 4] + 0.2*Z  + 0.2*rnorm(n)
    W <- as.numeric(W.star > 0.5)

    lambda.failure <- 2*X[, 1] + X[, 2] + 2*X[, 4] + 4 + 2 * (sqrt(X[, 1]) - 0.3) * W
    
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 7 
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    
    cate <- rep(NA, n)
    rand_u = matrix(runif(n * n.mc), n, n.mc)
    
    for (i in 1:n) {
      lambda.failure.0 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 4
      lambda.failure.1 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 4 + 2 * (sqrt(X[i, 1]) - 0.3)
      ft0 <- rpois(n.mc, lambda.failure.0)
      ft1 <- rpois(n.mc, lambda.failure.1)
      
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  
  else if (dgp == "type208") {
    # 90%+ censoring, discrete
    # Instrumental variable
    # Weak instrument analysis

    Y.max <- 8
    Z <- as.numeric(X[, 5] > 0.5)
    
    W.star <- 0.5*X[, 4] + 0.1*Z  + 0.2*rnorm(n)
    W <- as.numeric(W.star > 0.5)

    lambda.failure <- 2*X[, 1] + X[, 2] + 2*X[, 4] + 4 + 2 * (sqrt(X[, 1]) - 0.3) * W
    
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 7 
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    
    cate <- rep(NA, n)
    rand_u = matrix(runif(n * n.mc), n, n.mc)
    
    for (i in 1:n) {
      lambda.failure.0 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 4
      lambda.failure.1 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 4 + 2 * (sqrt(X[i, 1]) - 0.3)
      ft0 <- rpois(n.mc, lambda.failure.0)
      ft1 <- rpois(n.mc, lambda.failure.1)
      
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  
  else if (dgp == "type210") {
    # 90%+ censoring, discrete
    # Instrumental variable
    # Weak Instrument analysis

    Y.max <- 8
    Z <- as.numeric(X[, 5] > 0.5)
    
    W.star <- 0.5*X[, 4] + 0.4*Z  + 0.2*rnorm(n)
    W <- as.numeric(W.star > 0.5)

    lambda.failure <- 2*X[, 1] + X[, 2] + 2*X[, 4] + 6 + 2 * (sqrt(X[, 1]) - 0.3) * W
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 4
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)

    rand_u = matrix(runif(n * n.mc), n, n.mc)
    
    for (i in 1:n) {
      lambda.failure.0 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 6
      lambda.failure.1 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 6 + 2 * (sqrt(X[i, 1]) - 0.3)
      ft0 <- rpois(n.mc, lambda.failure.0)
      ft1 <- rpois(n.mc, lambda.failure.1)
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  else if (dgp == "type211") {
    # 90%+ censoring, discrete
    # Instrumental variable
    # Weak Instrument analysis

    Y.max <- 8
    Z <- as.numeric(X[, 5] > 0.5)
    
    W.star <- 0.5*X[, 4] + 0.3*Z  + 0.2*rnorm(n)
    W <- as.numeric(W.star > 0.5)

    lambda.failure <- 2*X[, 1] + X[, 2] + 2*X[, 4] + 6 + 2 * (sqrt(X[, 1]) - 0.3) * W
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 4
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)

    rand_u = matrix(runif(n * n.mc), n, n.mc)
    
    for (i in 1:n) {
      lambda.failure.0 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 6
      lambda.failure.1 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 6 + 2 * (sqrt(X[i, 1]) - 0.3)
      ft0 <- rpois(n.mc, lambda.failure.0)
      ft1 <- rpois(n.mc, lambda.failure.1)
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
  else if (dgp == "type212") {
    # 90%+ censoring, discrete
    # Instrumental variable
    # Weak Instrument analysis

    Y.max <- 8
    Z <- as.numeric(X[, 5] > 0.5)
    
    W.star <- 0.5*X[, 4] + 0.2*Z  + 0.2*rnorm(n)
    W <- as.numeric(W.star > 0.5)

    lambda.failure <- 2*X[, 1] + X[, 2] + 2*X[, 4] + 6 + 2 * (sqrt(X[, 1]) - 0.3) * W
    failure.time <- rpois(n, lambda = lambda.failure)
    
    lambda.censor <- 4
    censor.time <- rpois(n, lambda = lambda.censor)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    cate <- rep(NA, n)

    rand_u = matrix(runif(n * n.mc), n, n.mc)
    
    for (i in 1:n) {
      lambda.failure.0 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 6
      lambda.failure.1 <- 2 * X[i, 1] + X[i, 2] + 2 * rand_u[i, ] + 6 + 2 * (sqrt(X[i, 1]) - 0.3)
      ft0 <- rpois(n.mc, lambda.failure.0)
      ft1 <- rpois(n.mc, lambda.failure.1)
      cate[i] <- mean(pmin(ft1, Y.max) - pmin(ft0, Y.max))
    }
  }
 
  list(X = X, Y = Y, W = W, D = D, cate = cate, dgp = dgp, Y.max = Y.max)
}




base_dir = "examples_data/"



for (datatype in c(1)) {

  start_time <- Sys.time()

  n = 5000
  n_test = 5000
  p = 5
  dgp = paste0("type", datatype)
  print(dgp)
  output_dir = paste0(base_dir, "data/", dgp)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  l <- vector("list", reps+1)
  for (seed in 0:reps) {
    rep_start_time <- Sys.time()
    set.seed(seed)

    if (datatype %in% c(7)) {
      X.train <- matrix(runif(n * (p+2)), n, p+2)
    } else {
      X.train <- matrix(runif(n * p), n, p)
    }
    data = generate_cases(n, p, dgp = dgp, n.mc = n.mc, X = X.train)
    Y.max = data$Y.max

    write.csv(data, paste0(output_dir, "/n_", n, "_seed_", seed, ".csv"))

    l[[seed+1]] = (100*sum(data$D == 0) / n)


    set.seed(10000+seed)

    if (datatype %in% c(7)) {
      X.test <- matrix(runif(n_test * (p+2)), n_test, p+2)
    } else {
      X.test <- matrix(runif(n_test * p), n_test, p)
    }
    data = generate_cases(n_test, p, dgp = dgp, n.mc = n.mc, X = X.test)

    write.csv(data, paste0(output_dir, "/n_", n_test, "_seed_", 10000+seed, "_test.csv"))
    rep_end_time <- Sys.time()
    print(paste("Type", datatype, "rep", seed, "took", as.integer(as.numeric(rep_end_time) - as.numeric(rep_start_time)), "seconds"))
  }

  print(mean(unlist(l)))

  for (seed in 0:n_percs) {

    set.seed(5000+seed)

    if (datatype %in% c(7)) {
      X.test <- matrix(seed / n_percs, 3, p+2)
    } else {
      X.test <- matrix(seed / n_percs, 3, p)
    }

    data = generate_cases(3, p, dgp = dgp, n.mc = n.mc, X = X.test)
    print(paste(dgp, 5000+seed))
    write.csv(data, paste0(output_dir, "/n_", n_test, "_seed_", 5000+seed, "_test.csv"))
  }

  end_time <- Sys.time()
  print(paste("Type", datatype, "took", as.integer(as.numeric(end_time) - as.numeric(start_time)), "seconds"))
}


# MIMIC data should be already pre-processed and in data dir 

preprocessed_filename = "fitters_table.csv"

for (datatype in c(140,141,142,143,144)) {

  start_time <- Sys.time()

  # In MIMIC data n = n_test
  n = 5000
  n_test = 5000

  n_folds = 5
  p = 36
  dgp = paste0("type", datatype)
  print(dgp)

  output_dir = paste0(base_dir, "data/", dgp)
  data_dir = paste0(base_dir, "data/")

  fitters_df <- read.csv(paste0(data_dir, preprocessed_filename), row.names = 1)

  seed = 100
  set.seed(seed)
  fitters_df <- fitters_df[sample(nrow(fitters_df)), ]

  fitters_df$J <- NULL
  fitters_df$X <- NULL
  fitters_df$pid <- NULL
  full_data = generate_cases(nrow(fitters_df), p, dgp = dgp, n.mc = n.mc, X = fitters_df)

  fold_size <- floor(nrow(fitters_df) / n_folds)
  for (i in 1:(n_folds)) {
    if (i == n_folds) {
      full_data$fold[((i-1) * fold_size + 1):nrow(fitters_df)] <- i
    } else {
      full_data$fold[((i-1) * fold_size + 1):(i * fold_size)] <- i
    }
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  write.csv(full_data, paste0(output_dir, "/full_data", "_seed_", seed, ".csv"))

  for (seed in 1:n_folds) {
    rep_start_time <- Sys.time()

    fold_indices <- which(full_data$fold == seed)
    data <- lapply(full_data, function(x) {
      if (is.matrix(x) || is.data.frame(x)) {
        return(x[fold_indices, ])
      } else if (length(x) > 1) {
        return(x[fold_indices])
      } else {
        return(x)
      }
    })

    data$fold <- NULL

    Y.max = data$Y.max
    write.csv(data, paste0(output_dir, "/n_", n, "_seed_", seed, ".csv"))


    set.seed(10000+seed)
    if (seed == n_folds) {
      fold_indices <- which(full_data$fold == 1)
    } else {
      fold_indices <- which(full_data$fold == (seed+1))
    }

    data <- lapply(full_data, function(x) {
      if (is.matrix(x) || is.data.frame(x)) {
        return(x[fold_indices, ])
      } else if (length(x) > 1) {
        return(x[fold_indices])
      } else {
        return(x)
      }
    })

    data$fold <- NULL

    write.csv(data, paste0(output_dir, "/n_", n_test, "_seed_", 10000+seed, "_test.csv"))
    rep_end_time <- Sys.time()
    print(paste("Type", datatype, "rep", seed, "took", as.integer(as.numeric(rep_end_time) - as.numeric(rep_start_time)), "seconds"))
  }

  end_time <- Sys.time()
  print(paste("Type", datatype, "took", as.integer(as.numeric(end_time) - as.numeric(start_time)), "seconds"))
}





# Instrumental Variable


for (datatype in c(200)) {

  start_time <- Sys.time()

  n = 5000
  n_test = 5000
  p = 3
  dgp = paste0("type", datatype)
  print(dgp)

  output_dir = paste0(base_dir, "data/", dgp)

  if (!dir.exists(output_dir)) {
      dir.create(output_dir)
  }

  for (seed in start_reps:reps) {
    rep_start_time <- Sys.time()
    set.seed(seed)

    if (datatype %in% c(200, 201, 202, 203, 204, 205, 206, 207, 208, 210, 211, 212)) {
      X.train <- matrix(runif(n * (p+2)), n, p+2)
    } else {
      X.train <- matrix(runif(n * p), n, p)
    }
    
    data = generate_cases(n, p, dgp = dgp, n.mc = n.mc, X = X.train)
    Y.max = data$Y.max

    write.csv(data, paste0(output_dir, "/n_", n, "_seed_", seed, ".csv"))


    set.seed(10000+seed)

    if (datatype %in% c(200, 201, 202, 203, 204, 205, 206, 207, 208, 210, 211, 212)) {
      X.test <- matrix(runif(n_test * (p+2)), n_test, p+2)
    } else {
      X.test <- matrix(runif(n_test * p), n_test, p)
    }
    
    data = generate_cases(n_test, p, dgp = dgp, n.mc = n.mc, X = X.test)
    
    write.csv(data, paste0(output_dir, "/n_", n_test, "_seed_", 10000+seed, "_test.csv"))
    rep_end_time <- Sys.time()
    print(paste("Type", datatype, "rep", seed, "took", as.integer(as.numeric(rep_end_time) - as.numeric(rep_start_time)), "seconds"))
  }

  for (seed in 0:n_percs) {

    set.seed(5000+seed)

    if (datatype %in% c(200, 201, 202, 203, 204, 205, 206, 207, 208, 210, 211, 212)) {
      X.test <- matrix(seed / n_percs, 3, p+2)
    } else {
      X.test <- matrix(seed / n_percs, 3, p)
    }

    data = generate_cases(3, p, dgp = dgp, n.mc = n.mc, X = X.test)
    print(paste(dgp, 5000+seed))
    write.csv(data, paste0(output_dir, "/n_", n_test, "_seed_", 5000+seed, "_test.csv"))
  }

  end_time <- Sys.time()
  print(paste("Type", datatype, "took", as.integer(as.numeric(end_time) - as.numeric(start_time)), "seconds"))
}