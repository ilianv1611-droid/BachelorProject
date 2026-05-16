################################################################################
# Title: Test procedure on ACTG175
# Author: Ilian Verlee
# Date: 28-02-2026
################################################################################

library(car)
library(glmnet)
library(BART)
library(parallel)
library(pbapply)

data("ACTG175")
data <- ACTG175


baseline_vars <- c( "age", "wtkg", "hemo", "homo", "karnof", "drugs"
                   , "oprior", "z30", "preanti", "race","gender", "strat",
                   "symptom", "treat", "cd40", "cd80")

data$hemo     <- as.factor(data$hemo)
data$homo     <- as.factor(data$homo)
data$drugs    <- as.factor(data$drugs)
data$oprior   <- as.factor(data$oprior)
data$z30      <- as.factor(data$z30)
data$race     <- as.factor(data$race)
data$gender   <- as.factor(data$gender)
data$strat    <- as.factor(data$strat)
data$symptom  <- as.factor(data$symptom)
data$treat    <- as.numeric(data$treat)

#Defining functions:


#TODO fix deze fuck ah functies, want model_matrix maakt wacky shit


#ATE calculator for a linear model using BIC, CV and EBIC


ATECalculator_lm <- function(data, x, y, outcome,
                             type = c("BIC", "CV", "EBIC")){
  type <- match.arg(type, c("BIC", "CV", "EBIC"))
  
  n <- length(y)
  p <- ncol(x)
  
  p.fac <- rep(1, ncol(x))
  p.fac[which(colnames(x) == "treat")] <- 0
  
  lambda <- NULL
  
  if (type == "CV") {
    
    cv.out <- cv.glmnet(x, y, alpha = 1, nfolds = 5)
    lambda <- cv.out$lambda.min
    fit <- cv.out$glmnet.fit 
    
  } else if (type == "BIC") {
    
    fit <- glmnet(x, y, alpha = 1)
    
    pred <- predict(fit, newx = x)
    
    rss <- colSums((matrix(y, n, length(fit$lambda)) - pred)^2)
    
    df <- fit$df
    
    BIC_vals <- n * log(rss / n) + log(n) * df
    
    lambda <- fit$lambda[which.min(BIC_vals)]
  }
  else if (type == "EBIC") {
    
    fit <- glmnet(x, y, alpha = 1)
    
    pred <- predict(fit, newx = x)
    
    rss <- colSums((matrix(y, n, length(fit$lambda)) - pred)^2)
    
    df <- fit$df
    
    log_binom <- lgamma(p + 1) - lgamma(df + 1) - lgamma(p - df + 1)
    
    EBIC_vals <- n * log(rss / n) + log(n) * df + 2 * log_binom
    
    lambda <- fit$lambda[which.min(EBIC_vals)]
  }
  
  
  fit <- glmnet(x, y, alpha = 1, lambda = lambda, penalty.factor = p.fac)
  
  
  coef_lasso <- coef(fit)
  
  ATE_lasso <- coef_lasso["treat", ]
  
  vars <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  vars <- vars[vars != "(Intercept)"]
  
  vars_new <- unique(c(vars, "treat"))
  
  formula_selected <- reformulate(vars_new, response = outcome)
  df <- data.frame(y,x)
  colnames(df)[1] <- outcome
  
  m <- lm(formula_selected, df)
  
  ATE_post <- m$coefficients["treat"]
  
  return(list(ATE_post = ATE_post,
              ATE_lasso = ATE_lasso))
}












#LASSO and post-LASSO p-values using permutation, semi-permutation or
#randomization tests.


dual_permutation_test <- function(data, permut, outcome,
                                  testType = c("perm", "sem", "rand"),
                                  type = c("BIC", "CV", "EBIC")) {
  testType <- match.arg(testType)
  type <- match.arg(type)
  n <- nrow(data)
  
  form <- reformulate(".", response = outcome)
  x <- model.matrix(form, data)[, -1]
  y <- data[[outcome]]
  
  or <- ATECalculator_lm(data, x, y, outcome, type)
  post_lasso_or <- or$ATE_post
  lasso_or <- or$ATE_lasso
  
  k_post_lasso <- numeric(permut)
  k_lasso <- numeric(permut)
  
  A_col <- which(colnames(x) == "treat")
  
  for (i in 1:permut) {
    
    if (i %% 1000 == 0) print(i)
    
    data_perm <- data
    x_perm <- x
    
    if (testType == "perm") {
      permuted <- sample(data$treat)
      data_perm$treat <- permuted
      x_perm[,A_col] <- permuted
      
    } else if (testType == "sem") {
      par <- mean(data$treat)
      permuted <- c(rbinom(n, 1, par))
      data_perm$treat <- permuted
      x_perm[,A_col] <- permuted
      
    } else if (testType == "rand") {
      permuted <- c(rbinom(n, 1, 0.5))
      data_perm$treat <- permuted
      x_perm[,A_col] <- permuted
    }
    perm <- ATECalculator_lm(data_perm, x_perm, y, outcome, type)
    
    k_post_lasso[i] <- perm$ATE_post
    k_lasso[i] <- perm$ATE_lasso
    
    if (is.na(k_post_lasso[i])) k_post_lasso[i] <- 0
    if (is.na(k_lasso[i])) k_lasso[i] <- 0
  }
  
  p_post_lasso <- (1+sum(abs(k_post_lasso) >= abs(post_lasso_or)))/(permut+1)
  p_lasso <- (1+sum(abs(k_lasso) >= abs(lasso_or)))/(permut+1)
  
  list(
    p_post_lasso = p_post_lasso,
    p_lasso = p_lasso
  )
}
















ATE_Calculator_logit <- function(data, x , y , outcome,
                                 type = c("CV", "BIC", "EBIC" , "AIC")) {
  
  type <- match.arg(type)
  
  
  n <- length(y)
  p <- ncol(x)
  
  
  if (type == "CV") {
    
    cv.out <- cv.glmnet(x, y, alpha = 1, family = "binomial", nfolds = 5)
    lambda <- cv.out$lambda.min
    fit <- cv.out$glmnet.fit
    
  } else {
    
    fit <- glmnet(x, y, alpha = 1, family = "binomial")
    
    loglik <- fit$nulldev * (1 - fit$dev.ratio)
    
    beta <- as.matrix(coef(fit))
    df <- apply(beta[-1, ], 2, function(b) sum(b != 0))
    
    if (type == "BIC") {
      crit <- -loglik + log(n) * df
    }
    
    
    if (type == "AIC") {
      crit <- -loglik + 2 * df 
    }
    
    if (type == "EBIC") {
      log_binom <- lgamma(p + 1) - lgamma(df + 1) - lgamma(p - df + 1)
      crit <- -loglik + log(n) * df + 2 * log_binom
    }
    
    lambda <- fit$lambda[which.min(crit)]
  }
  
  p.fac <- rep(1, ncol(x))
  p.fac[which(colnames(x) == "treat")] <- 0
  
  
  fit <- glmnet(x, y, alpha = 1, family = "binomial", lambda = lambda,
                penalty.factor = p.fac)
  
  
  data1 <- data
  data1$treat <- 1
  
  data0 <- data
  data0$treat <- 0
  
  form <- reformulate(".", response = outcome)
  
  x1 <- model.matrix(form, data1)[, -1]
  x0 <- model.matrix(form, data0)[, -1]
  
  p1 <- predict(fit, newx = x1, s = lambda, type = "response")
  p0 <- predict(fit, newx = x0, s = lambda, type = "response")
  
  ATE_lasso <- mean(p1 - p0)
  
  coef_lasso <- coef(fit)
  
  vars <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  vars <- vars[vars != "(Intercept)"]
  
  vars_new <- unique(c(vars, "treat"))
  
  formula_selected <- reformulate(vars_new, response = outcome)
  df <- data.frame(y,x)
  colnames(df)[1] <- outcome
  
  
  fit <- glm(formula_selected,  family = "binomial", data = df)
  
  data1 <- df
  data1$treat <- 1
  
  data0 <- df
  data0$treat <- 0
  
  p1 <- predict(fit, newdata = data1,  type = "response")
  p0 <- predict(fit, newdata = data0, type = "response")
  
  ATE_post <- mean(p1 - p0)
  
  
  return(list(ATE_post = ATE_post,
              ATE_lasso = ATE_lasso))
}









dual_permutation_test_logit <- function(data, permut, outcome,
                                        testType = c("perm", "sem", "rand"),
                                        type = c("BIC", "CV", "EBIC", "AIC")) {
  testType <- match.arg(testType)
  type <- match.arg(type)
  n <- nrow(data)
  
  form <- reformulate(".", response = outcome)
  x <- model.matrix(form, data)[, -1]
  y <- data[[outcome]]
  
  or <- ATE_Calculator_logit(data, x, y, outcome, type)
  post_lasso_or <- or$ATE_post
  lasso_or <- or$ATE_lasso
  
  k_post_lasso <- numeric(permut)
  k_lasso <- numeric(permut)
  
  A_col <- which(colnames(x) == "treat")
  
  for (i in 1:permut) {
    
    if (i %% 1000 == 0) print(i)
    
    data_perm <- data
    x_perm <- x
    
    if (testType == "perm") {
      permuted <- sample(data$treat)
      data_perm$treat <- permuted
      x_perm[,A_col] <- permuted
      
    } else if (testType == "sem") {
      par <- mean(data$treat)
      permuted <- c(rbinom(n, 1, par))
      data_perm$treat <- permuted
      x_perm[,A_col] <- permuted
      
    } else if (testType == "rand") {
      permuted <- c(rbinom(n, 1, 0.5))
      data_perm$treat <- permuted
      x_perm[,A_col] <- permuted
    }
    perm <- ATE_Calculator_logit(data_perm, x_perm, y, outcome, type)
    
    k_post_lasso[i] <- perm$ATE_post
    k_lasso[i] <- perm$ATE_lasso
    
    if (is.na(k_post_lasso[i])) k_post_lasso[i] <- 0
    if (is.na(k_lasso[i])) k_lasso[i] <- 0
  }
  
  p_post_lasso <- (1+sum(abs(k_post_lasso) >= abs(post_lasso_or)))/(permut+1)
  p_lasso <- (1+sum(abs(k_lasso) >= abs(lasso_or)))/(permut+1)
  
  list(
    p_post_lasso = p_post_lasso,
    p_lasso = p_lasso
  )
}

two_prop_z <- function(data, outcome){
  n1 <- sum(data$treat == 1)
  n0 <- sum(data$treat == 0)
  
  y1 <- sum(data[[outcome]][data$treat == 1] == 1)
  y0 <- sum(data[[outcome]][data$treat == 0] == 1)
  return(prop.test(x = c(y1, y0), n = c(n1, n0), correct = FALSE)$p.value)
}






#---------------------Test cd420 lineair model----------------------------------


# On full data


data_cd4 <- data[,c("cd420", baseline_vars)]
permutations <- 10000
type<- "BIC"
p <- dual_permutation_test(data_cd4, permutations, "cd420", testType = "rand", type)
summary(lm(cd420 ~ treat, data_cd4))
model_step<-step(lm(cd420 ~.,data_cd4))
summary(model_step)
p


# Selected model on 99.9% of the data



  # Help function to ensure at least 1 of each factor is inside


ensure_factor_levels <- function(data, prop = 0.01) {
  
  n <- nrow(data)
  target_size <- floor(prop * n)
  
  mandatory_idx <- c()
  
  
  factor_cols <- names(data)[sapply(data, is.factor)]
  
  for(col in factor_cols) {
    
    levs <- levels(data[[col]])
    
    for(lv in levs) {
      
      idx <- which(data[[col]] == lv)
      
      mandatory_idx <- c(mandatory_idx, sample(idx, 1))
    }
  }
  
  mandatory_idx <- unique(mandatory_idx)
  
  
  remaining <- setdiff(seq_len(n), mandatory_idx)
  
  extra_needed <- max(0, target_size - length(mandatory_idx))
  
  extra_idx <- sample(remaining, extra_needed)
  
  test_idx <- c(mandatory_idx, extra_idx)
  
  train_idx <- setdiff(seq_len(n), test_idx)
  
  list(
    train_data = data[train_idx, ],
    test_data  = data[test_idx, ]
  )
}









  # ATE calculator for selected model


ATECalculator_lm_true_model <- function(data){
  m   <- lm(cd420 ~ strat + treat + 
              cd40, data= data)
  ATE <- m$coefficients["treat"]
  
  return(ATE)
}









  # Permutation test for this model


permutation_test_true_model <- function(data, permut,
                                        testType = c("perm", "sem", "rand")) {
  testType <- match.arg(testType)
  n <- nrow(data)
  
  
  or <- ATECalculator_lm_true_model(data)
  
  k <- numeric(permut)
  
  for (i in 1:permut) {
    
    if (i %% 1000 == 0) print(i)
    
    data_perm <- data
    
    if (testType == "perm") {
      permuted <- c(sample(data[,"treat"]))
      data_perm[,"treat"] <- permuted
      
    } else if (testType == "sem") {
      par <- mean(data[,"treat"])
      permuted <- c(rbinom(n, 1, par))
      data_perm[,"treat"] <- permuted
      
    } else if (testType == "rand") {
      permuted <- c(rbinom(n, 1, 0.5))
      data_perm[,"treat"]<- permuted
    }
    k[i] <- ATECalculator_lm_true_model(data_perm)
    
    
    if (is.na(k[i])) k[i] <- 0
  }
  
  p <- mean(abs(k) >= abs(or))
  
  return(p)
}









# Parallell for prop = 0.02


tester.seed <- function(data_cd4, prop = 0.02, permut = 10000){
  split <- ensure_factor_levels(data_cd4, prop)
  
  train_data <- split$train_data
  test_data  <- split$test_data
  
  
  dualATE <- dual_permutation_test(test_data, permutations,"cd420", testType = "rand", "BIC")
  trueATE <- permutation_test_true_model(test_data, 10000, "perm")
  
  p_LASSO      <-  dualATE$p_lasso
  p_post_LASSO <- dualATE$p_post_lasso
  
  m       <-lm(cd420 ~ treat , data= test_data)
  p_m     <-summary(m)$coefficients["treat", "Pr(>|t|)"]
  
  l      <-lm(cd420 ~ strat + treat + 
                 cd40 + karnof , data= test_data)
  p_l     <-summary(m)$coefficients["treat", "Pr(>|t|)"]
  
  
  p <- data.frame(
    p_LASSO      = p_LASSO,
    p_post_LASSO = p_post_LASSO,
    p_perm       = trueATE,
    p_rlm        = p_m,
    p_tlm        = p_l
  )
}



n.seed <- 500
prop   <- 0.02

params <- expand_grid(
  seed = 1:n.seed
)
param_list <- split(params, seq_len(nrow(params)))

n_cores <- max(1, detectCores() - 2)
cl <- makeCluster(n_cores)


clusterEvalQ(cl, {
  library(glmnet)
  library(tidyverse)
})


clusterExport(cl, varlist = c(
  "data_cd4",
  "permutations",
  "ensure_factor_levels",
  "dual_permutation_test",
  "permutation_test_true_model",
  "tester.seed",
  "ATECalculator_lm",
  "ATECalculator_lm_true_model"
))
# ----- Run simulation -----

resultsLMrand <- pblapply(param_list, cl = cl, FUN = function(param) {
  tester.seed(data_cd4 = data_cd4, permut =  permutations)
})

stopCluster(cl)




results.sim <- data.frame(do.call(rbind, resultsLMrand))
colnames(results.sim) <- c("LASSO",
                           "Post_LASSO",
                           "Perm_true_model",
                           "restricted_model",
                           "true_model"
)

p_mean <- c(colMeans(results.sim[,], na.rm = TRUE))
power <- colMeans(results.sim[, ]
                  <= alpha)

resultsRandNorm <- if(bet1 ==0 ) data.frame( p_mean, type1 = power) else data.frame(p_mean, power)

rownames(results.sim) <- NULL

print(resultsRandNorm)









#---------------------Test cens logistich model---------------------------------


data_cens <- data[,c("cens", baseline_vars)]
outcome<- "cens"
permutations <- 10000
type<- "BIC"
dual_permutation_test_logit(data_cens, permutations, outcome, testType = "rand", type)
two_prop_z(data_cens, outcome)









# Selected model on 99.9% of the data



# Help function to ensure at least 1 of each factor is inside


ensure_factor_levels <- function(data, prop = 0.01) {
  
  n <- nrow(data)
  target_size <- floor(prop * n)
  
  mandatory_idx <- c()
  
  
  factor_cols <- names(data)[sapply(data, is.factor)]
  
  for(col in factor_cols) {
    
    levs <- levels(data[[col]])
    
    for(lv in levs) {
      
      idx <- which(data[[col]] == lv)
      
      mandatory_idx <- c(mandatory_idx, sample(idx, 1))
    }
  }
  
  mandatory_idx <- unique(mandatory_idx)
  
  
  remaining <- setdiff(seq_len(n), mandatory_idx)
  
  extra_needed <- max(0, target_size - length(mandatory_idx))
  
  extra_idx <- sample(remaining, extra_needed)
  
  test_idx <- c(mandatory_idx, extra_idx)
  
  train_idx <- setdiff(seq_len(n), test_idx)
  
  list(
    train_data = data[train_idx, ],
    test_data  = data[test_idx, ]
  )
}









# ATE calculator for selected model


ATE_Calculator_logit_true_model <- function(data){
  
  fit <- glm(cens ~ treat + cd40 + cd80 + symptom,
             family = binomial,
             data)
  
  data1 <- data; data1[,"treat"] <- 1
  data0 <- data; data0[,"treat"] <- 0
  
  p1 <- predict(fit, newdata = data1, type = "response")
  p0 <- predict(fit, newdata = data0, type = "response")
  return(mean(p1-p0))
}









# Permutation test for this model


permutation_test_logit_true_model<- function(data, permut,
                                             testType = c("perm", "sem", "rand")){
  
  n <- nrow(data)
  
  
  k<- numeric(permut)
  
  ATE_or <- ATE_Calculator_logit_true_model(data)
  for (i in 1:permut) {
    
    if (i %% 1000 == 0) print(i)
    
    data_perm <- data
    
    if (testType == "perm") {
      permuted <- c(sample(data[,"treat"]))
      data_perm[,"treat"] <- permuted
      
    } else if (testType == "sem") {
      par <- mean(data[,"treat"])
      permuted <- c(rbinom(n, 1, par))
      data_perm[,"treat"] <- permuted
      
    } else if (testType == "rand") {
      permuted <- c(rbinom(n, 1, 0.5))
      data_perm[,"treat"] <- permuted
    }
    k[i] <- ATE_Calculator_logit_true_model(data_perm)
    
    if (is.na(k[i])) k[i] <- 0
  }
  return(mean(abs(k) >= abs(ATE_or)))
}


split <- ensure_factor_levels(data_cens, prop)

train_data <- split$train_data
test_data  <- split$test_data





m <- glm(cens ~ ., family = binomial, data = train_data)

step_model <- step(m)

top4 <- rownames(
  summary(step_model)$coefficients[-1, ]
)[
  order(summary(step_model)$coefficients[-1, "Pr(>|z|)"])[1:4]
]
top4














# Parallell for prop = 0.02


tester.seed <- function(data, prop = 0.02, permut = 10000){
  split <- ensure_factor_levels(data_cens, prop)
  
  train_data <- split$train_data
  test_data  <- split$test_data
  
  
  #dualATE <- dual_permutation_test_logit(test_data, permutations,"cens", testType = "perm", "BIC")
  trueATE <- permutation_test_logit_true_model(test_data, 10000, "perm")
  
  two_p   <- two_prop_z(test_data, outcome)
  
  #p_LASSO      <-  dualATE$p_lasso
  #p_post_LASSO <- dualATE$p_post_lasso
  
  
  p <- data.frame(
    #p_LASSO      = p_LASSO,
    #p_post_LASSO = p_post_LASSO,
    p_perm       = trueATE,
    two_p        = two_p
  )
}



n.seed <- 500
prop   <- 0.02

params <- expand_grid(
  seed = 1:n.seed
)
param_list <- split(params, seq_len(nrow(params)))

n_cores <- max(1, detectCores() - 2)
cl <- makeCluster(n_cores)


clusterEvalQ(cl, {
  library(glmnet)
  library(tidyverse)
})


clusterExport(cl, varlist = c(
  "data_cens",
  "permutations",
  "ensure_factor_levels",
  "dual_permutation_test_logit",
  "permutation_test_logit_true_model",
  "tester.seed",
  "ATE_Calculator_logit_true_model",
  "ATE_Calculator_logit",
  "two_prop_z",
  "outcome"
))
# ----- Run simulation -----

resultsLMrand <- pblapply(param_list, cl = cl, FUN = function(param) {
  tester.seed(data = data_cens, permut =  permutations)
})

stopCluster(cl)




results.sim <- data.frame(do.call(rbind, resultsLMrand))
colnames(results.sim) <- c("Perm_true_model",
                           "two_prop"
)

p_mean <- c(colMeans(results.sim[,], na.rm = TRUE))
power <- colMeans(results.sim[, ]
                  <= alpha)

resultsRandNorm <- if(bet1 ==0 ) data.frame( p_mean, type1 = power) else data.frame(p_mean, power)

rownames(results.sim) <- NULL

print(resultsRandNorm)
