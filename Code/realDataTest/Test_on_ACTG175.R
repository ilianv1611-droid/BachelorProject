################################################################################
# Title: Test procedure on ACTG175
# Author: Ilian Verlee
# Date: 28-02-2026
################################################################################


library(glmnet)


data("ACTG175")
data <- ACTG175


baseline_vars <- c( "age", "wtkg", "hemo", "homo", "karnof", "drugs"
                   , "oprior", "z30", "preanti", "race","gender", "strat",
                   "symptom", "treat", "cd40", "cd80")

data$hemo     <- as.factor(data$hemo)
data$homo     <- as.factor(data$homo)
data$drugs    <- as.factor(data$drugs)
data$oprior   <- as.factor(data$oprior)
data$z30   <- as.factor(data$z30)
data$race     <- as.factor(data$race)
data$gender   <- as.factor(data$gender)
data$strat    <- as.factor(data$strat)
data$symptom  <- as.factor(data$symptom)
data$treat    <- as.factor(data$treat)


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
  
  ATE_lasso <- coef_lasso["treat1", ]
  
  vars <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  vars <- vars[vars != "(Intercept)"]
  
  vars_new <- unique(c(vars, "treat1"))
  
  formula_selected <- reformulate(vars_new, response = outcome)
  df <- data.frame(y,x)
  colnames(df)[1] <- outcome
  
  m <- lm(formula_selected, df)
  
  ATE_post <- m$coefficients["treat1"]
  
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
  
  data_perm <- data
  x_perm <- x
  
  for (i in 1:permut) {
    
    if (i %% 1000 == 0) print(i)
    
    data_perm <- data
    
    if (testType == "perm") {
      permuted <- c(sample(data$treat))
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
  
  p_post_lasso <- mean(abs(k_post_lasso) >= abs(post_lasso_or))
  p_lasso <- mean(abs(k_lasso) >= abs(lasso_or))
  
  list(
    p_post_lasso = p_post_lasso,
    p_lasso = p_lasso
  )
}
















ATE_Calculator_logit <- function(data, x , y , type = c("CV", "BIC", "EBIC" , "AIC")) {
  
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
  data1$A <- 1
  
  data0 <- data
  data0$A <- 0
  
  x1 <- model.matrix(Y ~ ., data1)[, -1]
  x0 <- model.matrix(Y ~ ., data0)[, -1]
  
  p1 <- predict(fit, newx = x1, s = lambda, type = "response")
  p0 <- predict(fit, newx = x0, s = lambda, type = "response")
  
  ATE_lasso <- mean(p1 - p0)
  
  coef_lasso <- coef(fit)
  
  vars <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  vars <- vars[vars != "(Intercept)"]
  
  vars_new <- unique(c(vars, "A"))
  
  formula_selected <- reformulate(vars_new, response = "Y")
  
  fit <- glm(formula_selected,  family = "binomial", data = data)
  
  p1 <- predict(fit, newdata = data1,  type = "response")
  p0 <- predict(fit, newdata = data0, type = "response")
  
  ATE_post <- mean(p1 - p0)
  
  
  return(list(ATE_post = ATE_post,
              ATE_lasso = ATE_lasso))
}









dual_permutation_test_logit <- function(data, permut,
                                        testType = c("perm", "sem", "rand"),
                                        type = c("BIC", "CV", "EBIC", "AIC")) {
  testType <- match.arg(testType)
  type <- match.arg(type)
  n <- nrow(data)
  
  x <- as.matrix(model.matrix(Y ~ ., data)[, -1])
  y <- data$Y
  
  or <- ATE_Calculator_logit(data, x, y, type)
  post_lasso_or <- or$ATE_post
  lasso_or <- or$ATE_lasso
  
  k_post_lasso <- numeric(permut)
  k_lasso <- numeric(permut)
  
  A_col <- which(colnames(x) == "A")
  
  
  
  for (i in 1:permut) {
    
    if (i %% 1000 == 0) print(i)
    
    data_perm <- data
    x_perm <- x
    
    if (testType == "perm") {
      permuted <- c(sample(data$A))
      data_perm$A <- permuted
      x_perm[,A_col] <- permuted
      
    } else if (testType == "sem") {
      par <- mean(data$A)
      permuted <- c(rbinom(n, 1, par))
      data_perm$A <- permuted
      x_perm[,A_col] <- permuted
      
    } else if (testType == "rand") {
      permuted <- c(rbinom(n, 1, 0.5))
      data_perm$A <- permuted
      x_perm[,A_col] <- permuted
    }
    perm <- ATE_Calculator_logit(data_perm, x_perm, y, type)
    
    k_post_lasso[i] <- perm$ATE_post
    k_lasso[i] <- perm$ATE_lasso
    
    if (is.na(k_post_lasso[i])) k_post_lasso[i] <- 0
    if (is.na(k_lasso[i])) k_lasso[i] <- 0
  }
  
  p_post_lasso <- mean(abs(k_post_lasso) >= abs(post_lasso_or))
  p_lasso <- mean(abs(k_lasso) >= abs(lasso_or))
  
  list(
    p_post_lasso = p_post_lasso,
    p_lasso = p_lasso
  )
}







#---------------------Test cd420 lineair model----------------------------------
data_cd4 <- data[,c("cd420", baseline_vars)]
data <- data_cd4
permutations <- 10000
dual_permutation_test(data_cd4, permutations, "cd420", testType = "perm", 
                      type=  "BIC")
summary(lm(cd420 ~ treat, data_cd4))
















#---------------------Test cens logistich model---------------------------------
data_cd4 <- data[,c("cens", baseline_vars)]
