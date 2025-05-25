library(rje)

generate_cases <- function (n, p, X = NULL, n.mc = 20000, dgp = c("type8", "type200"))
{
  dgp <- match.arg(dgp)

  if (dgp == "type8") {

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

  else if (dgp == "type200") {

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
  
  list(X = X, Y = Y, W = W, D = D, cate = cate, dgp = dgp, Y.max = Y.max)
}


n = 5000
n_test = 5000
p = 5
n.mc = 20000
datatype = 8
dgp = paste0("type", datatype)
print(dgp)
base_dir = "example_data/data/"
output_dir = paste0(base_dir, dgp)

if (!dir.exists(base_dir)) {
  dir.create(base_dir)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

seed = 1
set.seed(seed)

X.train <- matrix(runif(n * p), n, p)
data = generate_cases(n, p, dgp = dgp, n.mc = n.mc, X = X.train)
write.csv(data, paste0(output_dir, "/n_", n, "_seed_", seed, ".csv"))

set.seed(10000+seed)
X.test <- matrix(runif(n_test * p), n_test, p)
data = generate_cases(n_test, p, dgp = dgp, n.mc = n.mc, X = X.test)
write.csv(data, paste0(output_dir, "/n_", n_test, "_seed_", 10000+seed, "_test.csv"))


# Instrumental Variable

datatype = 200

n = 5000
n_test = 5000
p = 3
dgp = paste0("type", datatype)
output_dir = paste0(base_dir, dgp)
print(dgp)

if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

set.seed(seed)
X.train <- matrix(runif(n * (p+2)), n, p+2)
data = generate_cases(n, p, dgp = dgp, n.mc = n.mc, X = X.train)
write.csv(data, paste0(output_dir, "/n_", n, "_seed_", seed, ".csv"))

set.seed(10000+seed)
X.test <- matrix(runif(n_test * (p+2)), n_test, p+2)
data = generate_cases(n_test, p, dgp = dgp, n.mc = n.mc, X = X.test)
write.csv(data, paste0(output_dir, "/n_", n_test, "_seed_", 10000+seed, "_test.csv"))
