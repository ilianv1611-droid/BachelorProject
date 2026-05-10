################################################################################
# Title: simulation study of randomization and permutation tests in G-computation
# Author: Ilian Verlee
# Date: 28-02-2026
################################################################################

library(parallel)
library(pbapply)
library(glmnet)
library(tidyverse)


set.seed(100)



n <- 30
permutations <- 10000


#Simulation for a LM with 1 covariate which is normally distributed + treatment
simulateDataLM1 <- function(n,a, b, c){
  L <- rnorm(n, 0, 1)
  A <- rbinom(n, 1, 0.5)
  Y <- a+b*A+c*L + rnorm(n,0,0.5)
  return(data.frame(Y = Y, A = A, L = L))
}


#Simulation for a LM with 1 covariate which is exp distributed + treatment

simulateDataLMexp <- function(n, a, b, c){
  L <- rexp(n, 1)
  A <- rbinom(n, 1, 0.5)
  Y <- a+b*A+c*L + rnorm(n,0,0.5)
  return(data.frame(Y = Y, A = A, L = L))
}


#Simulation for a LM with k covariates which are normally distributed + treatment / No Noise

simulateDataLMRand <- function(n,inter, a, k, meanX, varX){
  A <- rbinom(n, 1, 0.5)
  X <- matrix(rnorm(n * k, meanX, varX), nrow = n, ncol = k)
  colnames(X) <- paste0("x", 1:k)
  b <- c(rep(1, k))
  Y <- inter + a * A + (X %*% b) + rnorm(n, 0, 0.5)
  return(data.frame(Y = Y, A, X))
}


#Simulation for a LM with k covariates which are normally distributed + treatment / With Noise

simulateDataLMRand_Noise <- function(n,inter, a, k, meanX, varX){
  A <- rbinom(n, 1, 0.5)
  X <- matrix(rnorm(n * k, meanX, varX), nrow = n, ncol = k)
  U <- matrix(rnorm(n *2 *  k, meanX, varX), nrow = n, ncol = 2*k)
  colnames(U) <- paste0("u", 1:(2*k))
  colnames(X) <- paste0("x", 1:k)
  b <- c(rep(1, k))
  Y <- inter + a * A + (X %*% b) + rnorm(n, 0, 0.5)
  return(data.frame(Y = Y, A, X, U))
  
}









# ATE Calculator using LASSO and post LASSO

ATECalculator_lm <- function(data, x, y, type = c("BIC", "CV", "EBIC")){
  type <- match.arg(type, c("BIC", "CV", "EBIC"))
  
  n <- length(y)
  p <- ncol(x)
  
  p.fac <- rep(1, ncol(x))
  p.fac[which(colnames(x) == "A")] <- 0
  
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
  
  ATE_lasso <- coef(fit)["A", ]
  
  
  coef_lasso <- coef(fit)
  
  vars <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  vars <- vars[vars != "(Intercept)"]
  
  vars_new <- unique(c(vars, "A"))
  
  formula_selected <- reformulate(vars_new, response = "Y")
  
  m <- lm(formula_selected, data)
  
  ATE_post <- m$coefficients["A"]
  
  return(list(ATE_post = ATE_post,
              ATE_lasso = ATE_lasso))
}


#Formula selector for a lm

formula_selector_lm <- function(data, x, y, type = c("BIC", "CV", "EBIC")){
  type <- match.arg(type, c("BIC", "CV", "EBIC"))
  
  n <- length(y)
  p <- ncol(x)
  
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
  
  
  fit <- glmnet(x, y, alpha = 1, lambda = lambda)
  
  
  coef_lasso <- coef(fit)
  
  vars <- rownames(coef_lasso)[coef_lasso[,1] != 0]
  vars <- vars[vars != "(Intercept)"]
  
  vars_new <- unique(c(vars, "A"))
  
  formula_selected <- reformulate(vars_new, response = "Y")
  
  return(formula_selected)
}




post_lasso_coef_p <- function(data, type = c("BIC", "CV", "EBIC")){
  type <- match.arg(type, c("BIC", "CV", "EBIC"))
  x <- model.matrix(Y ~ ., data)[, -1]
  y <- data$Y
  formula <-formula_selector_lm(data, x, y, type)
  m <- lm(formula, data)
  p <- summary(m)$coefficients["A", "Pr(>|t|)"]
  return(p)
}


#Testing for 1 covariate + treatment

testing_LM1 <- function(data, permut, testType = c("perm", "sem", "rand")){
  testType <- match.arg(testType)
  m <- lm(Y ~ A + L, data)
  estimate <- m$coefficients["A"]
  k <- rep(0,permut)
  if( testType == "perm"){
    for(i in 1:permut){
      m <- lm(Y ~ A + L, data.frame (Y =data$Y , A= sample(data$A), L =data$L))
      coef <- m$coefficients["A"]
      k[i]<- if_else(is.na(coef), 0, coef)
    }
  }
  else if(testType == "sem"){
    par <- mean(data$A)
    for(i in 1:permut){
      m <- lm(Y ~ A + L, data.frame (Y =data$Y , A= rbinom(n, 1, par), L =data$L))
      coef <- m$coefficients["A"]
      k[i]<- if_else(is.na(coef), 0, coef)
    }
  }
  else if(testType == "rand"){
    for(i in 1:permut){
      m <- lm(Y ~ A + L, data.frame (Y =data$Y , A= rbinom(n, 1, 0.5), L =data$L))
      coef <- m$coefficients["A"]
      k[i]<- if_else(is.na(coef), 0, coef)
    }
  }
  p <- mean(abs(k) >= abs(estimate))
  return(p)
}



#Test lasso and post_lasso where a new formula is selected for each permutation

dual_permutation_test <- function(data, permut,
                                  testType = c("perm", "sem", "rand"),
                                  type = c("BIC", "CV", "EBIC")) {
  
  testType <- match.arg(testType)
  type <- match.arg(type)
  n <- nrow(data)
  
  x <- model.matrix(Y ~ ., data)[, -1]
  y <- data$Y
  
  or <- ATECalculator_lm(data, x, y, type)
  post_lasso_or <- or$ATE_post
  lasso_or <- or$ATE_lasso
  
  k_post_lasso <- numeric(permut)
  k_lasso <- numeric(permut)
  
  A_col <- which(colnames(x) == "A")
  
  data_perm <- data
  x_perm <- x
  
  for (i in 1:permut) {
    
    if (i %% 1000 == 0) print(i)
    
    data_perm <- data
    
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
    perm <- ATECalculator_lm(data_perm, x_perm, y, type)
    
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



#Test for which only 1 formula is selected in the beginning using LASSO

permutation_test_1_formula <- function(data, permut,
                                      testType = c("perm", "sem", "rand"),
                                      type = c("BIC", "CV", "EBIC")) {
  
  testType <- match.arg(testType)
  type <- match.arg(type)
  n <- nrow(data)
  
  x <- model.matrix(Y ~ ., data)[, -1]
  y <- data$Y
  
  formula_selected <- formula_selector_lm(data, x, y, type)
  or <- lm(formula_selected, data)$coefficients["A"]
  
  
  k <- numeric(permut)
  
  A_col <- which(colnames(x) == "A")
  
  data_perm <- data
  
  for (i in 1:permut) {
    
    if (i %% 1000 == 0) print(i)
    
    data_perm <- data
    
    if (testType == "perm") {
      permuted <- c(sample(data$A))
      data_perm$A <- permuted
      
    } else if (testType == "sem") {
      par <- mean(data$A)
      permuted <- c(rbinom(n, 1, par))
      data_perm$A <- permuted
      
    } else if (testType == "rand") {
      permuted <- c(rbinom(n, 1, 0.5))
      data_perm$A <- permuted
    }
    k[i] <- lm(formula_selected, data_perm)$coefficients["A"]
    
    if (is.na(k[i])) k[i] <- 0
  }
  
  p <- mean(abs(k) >= abs(or))
  
  return(p)
}


histMaker_LM1 <- function(data, permut,param, testType = c("perm", "sem", "rand")){
  testType <- match.arg(testType, c("perm", "sem", "rand"))
  m <- lm(Y ~ A + L, data)
  estimate <- m$coefficients["A"]
  k <- rep(0,permut)
  
  if( testType == "perm"){
    for(i in 1:permut){
      m <- lm(Y ~ A + L, data.frame (Y =data$Y , A= sample(data$A), L =data$L))
      k[i]<- m$coefficients["A"]
    }
  }
  else if(testType == "sem"){
    par <- mean(data$A)
    for(i in 1:permut){
      m <- lm(Y ~ A + L, data.frame (Y =data$Y , A= rbinom(n, 1, par), L =data$L))
      k[i]<- m$coefficients["A"]
    }
  }
  else if(testType == "rand"){
    for(i in 1:permut){
      m <- lm(Y ~ A + L, data.frame (Y =data$Y , A= rbinom(n, 1, 0.5), L =data$L))
      k[i]<- m$coefficients["A"]
    }
  }
  hist(k, ylab="aantal permutaties",  xlab="geschatte ATE", breaks= seq(min(c(min(k, na.rm = TRUE)-1, 0)),max(c(max(k, na.rm = TRUE)+1, param)),0.05), main="Verdeling ATE")
  abline(v = estimate, col ='red')
}



histMaker_Dual <- function(data, permut, type = c("AIC", "BIC", "CV"), 
                           testType = c("perm", "sem", "rand"), 
                           true_param = NULL) {
  
  testType <- match.arg(testType)
  type <- match.arg(type)
  n <- nrow(data)
  
  x <- model.matrix(Y ~ ., data)[, -1]
  y <- data$Y
  A_col <- which(colnames(x) == "A")
  
  
  or <- ATECalculator_lm(data, x, y, type)
  obs_post <- or$ATE_post
  obs_lasso <- or$ATE_lasso
  
  k_post <- numeric(permut)
  k_lasso <- numeric(permut)
  
  for (i in 1:permut) {
    data_p <- data
    x_p <- x
    
    if (testType == "perm") {
      p_vals <- sample(data$A)
    } else if (testType == "sem") {
      p_vals <- rbinom(n, 1, mean(data$A))
    } else {
      p_vals <- rbinom(n, 1, 0.5)
    }
    
    data_p$A <- p_vals
    x_p[, A_col] <- p_vals
    
    res <- ATECalculator_lm(data_p, x_p, y, type)
    k_post[i] <- res$ATE_post
    k_lasso[i] <- res$ATE_lasso
  }
  
  par(mfrow = c(1, 2))
  
  
  hist(k_lasso, 
       main = paste("Lasso Null Dist (", type, ")"),
       xlab = "Permuted ATE", 
       col = "lightgrey", 
       border = "white",
       xaxt = "n",     
       breaks= seq(min(c(min(k_post, na.rm = TRUE)-1, 0))
                   ,max(c(max(k_post, na.rm = TRUE)+1, true_param)),0.05))
  
  axis(1, pos = 0)
  
  abline(v = obs_lasso, col = "red", lwd = 2, lty = 2)
  abline(v = -obs_lasso, col = "red", lwd = 2, lty = 2)
  
  
  hist(k_post, main = paste("Post-Lasso Null Dist (", type, ")"),
       xlab = "Permuted ATE", col = "lightblue", border = "white",
       yaxs = "i", xaxs ="i",
       breaks= seq(min(c(min(k_post, na.rm = TRUE)-1, 0)),max(c(max(k_post, na.rm = TRUE)+1, true_param)),0.05))
  abline(v = obs_post, col = "red", lwd = 2, lty = 2)
  abline(v = -obs_post, col = "red", lwd = 2, lty = 2)
  if(!is.null(true_param)) abline(v = true_param, col = "blue", lwd = 2) 

  
  par(mfrow = c(1, 1))
  
  return(list(p_lasso = mean(abs(k_lasso) >= abs(obs_lasso)),
              p_post = mean(abs(k_post) >= abs(obs_post))))
}

# -------------------------------- Permutation tests ------------------------------------



#type 1 error is meestal laag bij een 1 covariaat

#Type 2 error daaraantegen komt wel voor bij kleine treatment effecten rond de 0.1-0.5



bet0<- 0.1
bet1<- 1
mu<- 0.5


#Basic lineair model
LM <- simulateDataLM1(n = n, bet0, bet1, mu )

testing_LM1(LM, permutations, "perm")



dual_permutation_test(LM, permutations, "perm", "CV")
dual_permutation_test(LM, permutations, "perm", "BIC")
dual_permutation_test(LM, permutations, "perm", "EBIC")

permutation_test_1_formula(LM, permutations, "perm", "BIC")
post_lasso_coef_p(LM, "BIC")

histMaker_LM1(LM, permutations, bet1, "perm")

histMaker_Dual(LM, permutations, true_param = bet1, testType = "sem", type = "BIC")

#Lin model with exp distr. covariate
LMexp <- simulateDataLMexp(n = n, bet0, bet1, mu)


testing_LM1(LMexp, permutations, "perm")
testing_LM(LMexp, permutations, "perm")

histMaker_LM1(LMexp, permutations, bet1, "perm")



#Rand lin model with normal distr. covariate

LMrand <- simulateDataLMRand(n, inter = bet0, a = bet1, k = 3 , meanX = 10, varX = 2)
model <- lm(Y ~ ., data = LMrand)
summary(model)$coefficients["A", "Pr(>|t|)"]

dual_permutation_test(LMrand, permutations, "perm")

histMaker_LM(LMrand,permutations, bet1, "perm")



# --------------------------------- Semi-Permutation tests -------------------------------------


bet0<- 0.1
bet1<- 200
mu<- 200
LM <- simulateDataLMRand(n = n, bet0, bet1, k=3, meanX = 10, varX = 2 )


testing_LM1(LM, permutations, "sem")
testing_LM(LM, permutations, "sem")


histMaker_LM1(LM , permutations, bet1, "sem")
histMaker_LM(LM, permutations, bet1, "sem")


# -------------------------- Randomization tests ----------------------------------






bet0<- 0.1
bet1<- 0.1
mu<- -100
LM <- simulateDataLM1(n = n, bet0, bet1, mu )


testing_LM1(LM, permutations, "rand")
histMaker_LM1(LM, permutations, bet1, "rand")





# ----------------------------------- LM parallell of power for Rand. Norm without noise --------------------------------------------


# CV Time: permutations* n.seed * 1.5/10.000 min = time
# BIC Time: permutations * n.seed * 0.3 min / 10000 min = time
# EBIC Time: permutations * n.seed * 0.3 min / 10000 min = time


#Run at bet1 = 0 for type 1 error rate


# ------ Wrapper for parallel ------
bet0<-0.1
bet1<- 0.5
alpha <- 0.05
permutations <- 10000
type <- "BIC"


tester.seedlm <- function(n, permut, type){
  LMrand <- simulateDataLMRand(n, bet0, bet1, k=3, meanX = 10, varX = 2)
  model <- lm(Y ~ ., data = LMrand)
  
  perm <- dual_permutation_test(LMrand, permutations, "perm", type)
  sem <- dual_permutation_test(LMrand, permutations, "sem", type)
  rand <- dual_permutation_test(LMrand,permutations, "rand", type)
  
  perm_1_formula <- permutation_test_1_formula(LMrand, permutations, "perm", type)
  sem_1_formula <- permutation_test_1_formula(LMrand, permutations, "sem", type)
  rand_1_formula <- permutation_test_1_formula(LMrand, permutations, "rand", type)
  
  
  p<- t(matrix(c(perm$p_lasso, sem$p_lasso, rand$p_lasso, perm$p_post_lasso,
               sem$p_post_lasso, rand$p_post_lasso,
               summary(model)$coefficients["A", "Pr(>|t|)"], t.test(Y~A, data = LMrand)$p.value,
               perm_1_formula, sem_1_formula, rand_1_formula, post_lasso_coef_p(LMrand,type))))
  return(data.frame(p))
}

# ------ Parallel setup ------
n.seed <- 10

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


clusterExport(
  cl,
    c("bet0","bet1", "alpha", "tester.seedlm", "simulateDataLMRand", "dual_permutation_test", "permutations",
      "n", "ATECalculator_lm", "type", "formula_selector_lm", "permutation_test_1_formula", "post_lasso_coef_p")
)

# ----- Run simulation -----

resultsLMrand <- pblapply(param_list, cl = cl, FUN = function(param) {
  set.seed(param$seed)
  tester.seedlm(n = n, permut =  permutations, type = type)
})

stopCluster(cl)

results.sim <- data.frame(do.call(rbind, resultsLMrand))
colnames(results.sim) <- c("Perm_lasso", "Sem_Perm_lasso", "Rand_lasso" ,
                           "Perm_post_lasso", "Sem_Perm_post_lasso", "Rand_post_lasso", "cov_coef", "t_test",
                           "Perm_1_form", "Sem_1_form", "Rand_1_form", "post_lasso_p")

p_mean <- c(colMeans(results.sim[,], na.rm = TRUE))
power <- colMeans(results.sim[, c("Perm_lasso", "Sem_Perm_lasso", "Rand_lasso" ,
                                  "Perm_post_lasso", "Sem_Perm_post_lasso", "Rand_post_lasso", "cov_coef", "t_test",
                                  "Perm_1_form", "Sem_1_form", "Rand_1_form", "post_lasso_p")]
                  <= alpha)

resultsRandNorm <- if(bet1 ==0 ) data.frame( p_mean, type1 = power) else data.frame(p_mean, power)

results.sim <- data.frame(results.sim[, 0:6])
results.sim$seed <- 1:n.seed
rownames(results.sim) <- NULL

print(resultsRandNorm)











# ----------------------------------- LM parallell of power for Rand. Norm with noise --------------------------------------------


# permutations* n.seed * 1.5/10.000 min = time

#Run at bet1 = 0 for type 1 error rate


# ------ Wrapper for parallel ------
bet0<-0.1
bet1<- 0.5
alpha <- 0.05
permutations <- 10000
type <- "BIC"


tester.seedlm <- function(n, permut, type){
  LMrand <- simulateDataLMRand_Noise(n, bet0, bet1, k=3, meanX = 10, varX = 2)
  model <- lm(Y ~ ., data = LMrand)
  
  perm <- dual_permutation_test(LMrand, permutations, "perm", type)
  sem <- dual_permutation_test(LMrand, permutations, "sem", type)
  rand <- dual_permutation_test(LMrand,permutations, "rand", type)
  
  perm_1_formula <- permutation_test_1_formula(LMrand, permutations, "perm", type)
  sem_1_formula <- permutation_test_1_formula(LMrand, permutations, "sem", type)
  rand_1_formula <- permutation_test_1_formula(LMrand, permutations, "rand", type)
  
  
  p<- t(matrix(c(perm$p_lasso, sem$p_lasso, rand$p_lasso, perm$p_post_lasso,
                 sem$p_post_lasso, rand$p_post_lasso,
                 summary(model)$coefficients["A", "Pr(>|t|)"], t.test(Y~A, data = LMrand)$p.value,
                 perm_1_formula, sem_1_formula, rand_1_formula, post_lasso_coef_p(LMrand,type))))
  return(data.frame(p))
}

# ------ Parallel setup ------
n.seed <- 1000

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


clusterExport(
  cl,
  c("bet0","bet1", "alpha", "tester.seedlm", "simulateDataLMRand_Noise", "dual_permutation_test", "permutations",
    "n", "ATECalculator_lm", "type", "formula_selector_lm", "permutation_test_1_formula", "post_lasso_coef_p")
)

# ----- Run simulation -----

resultsLMrand <- pblapply(param_list, cl = cl, FUN = function(param) {
  set.seed(param$seed)
  tester.seedlm(n = n, permut =  permutations, type = type)
})

stopCluster(cl)

results.sim <- data.frame(do.call(rbind, resultsLMrand))
colnames(results.sim) <- c("Perm_lasso", "Sem_Perm_lasso", "Rand_lasso" ,
                           "Perm_post_lasso", "Sem_Perm_post_lasso", "Rand_post_lasso", "cov_coef", "t_test",
                           "Perm_1_form", "Sem_1_form", "Rand_1_form", "post_lasso_coef_p")

p_mean <- c(colMeans(results.sim[,], na.rm = TRUE))
power <- colMeans(results.sim[, c("Perm_lasso", "Sem_Perm_lasso", "Rand_lasso" ,
                                  "Perm_post_lasso", "Sem_Perm_post_lasso", "Rand_post_lasso", "cov_coef", "t_test",
                                  "Perm_1_form", "Sem_1_form", "Rand_1_form", "post_lasso_coef_p")]
                  <= alpha)

resultsRandNorm <- if(bet1 ==0 ) data.frame( p_mean, type1 = power) else data.frame(p_mean, power)

results.sim <- data.frame(results.sim[, 0:6])
results.sim$seed <- 1:n.seed
rownames(results.sim) <- NULL

print(resultsRandNorm)









# ---------------------------- Exp LM parallell ------------------------------------------------------

#Run at bet1 = 0 for type 1 error rate
bet0<- 0.1
bet1<- 0.3
mu <- 10
alpha <- 0.05
permutations <- 50000

#For bet1 = 0.3 and 50 000 permutations the power of Perm is the best with 0.23 (Sem/t_ = 0.20, Rand = 0.21) 
#and the type1_error is the lowest for Sem with 0.17 (Perm = 0.18, Rand/t_ = 0.19)

tester.seedlm <- function(n, permut){
  LMexp <- simulateDataLMexp(n, bet0, bet1, mu)
  model <- lm(Y ~ ., data = LMexp)
  p<- t(matrix(c(testing_LM_lasso(LMexp,permutations, "perm"), testing_LM_lasso(LMexp,permutations, "sem"),
                 testing_LM_lasso(LMexp,permutations, "rand"), testing_LM_post_lasso(LMexp,permutations, "perm"),
                 testing_LM_post_lasso(LMexp,permutations, "sem"), testing_LM_post_lasso(LMexp,permutations, "rand"),
                 summary(model)$coefficients["A", "Pr(>|t|)"])))
  return(data.frame(p))
}

# ------ Parallel setup ------
n.seed <- 100

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


clusterExport(
  cl,
  c("bet0","bet1", "mu", "alpha", "tester.seedlm", "simulateDataLMRand", "testing_LM_lasso", "testing_LM_post_lasso"
    , "permutations", "n",  "covariateSelector", "ATE_Calculator_lm")
)

# ----- Run simulation -----

resultsLMexp <- pblapply(param_list, cl = cl, FUN = function(param) {
  set.seed(param$seed)
  tester.seedlm(n = n, permut =  permutations)
})

stopCluster(cl)

results.sim <- data.frame(do.call(rbind, resultsLMexp))
colnames(results.sim) <- c("Perm", "Sem_Perm", "Rand" ,
                           "Perm2", "Sem_Perm2", "Rand2", "cov_coef")

p_mean <- c(colMeans(results.sim[,], na.rm = TRUE))
power <- colMeans(results.sim[, c("Perm", "Sem_Perm", "Rand", "Perm2", "Sem_Perm2", "Rand2" , "cov_coef")] <= alpha)
resultsRandexp <- data.frame( p_mean, power)

results.sim <- data.frame(results.sim[, 0:4])
results.sim$seed <- 1:n.seed
rownames(results.sim) <- NULL

head(resultsRand)

















# ------------------------------------------------ Logistic Regression -------------------------------------------------


simulateDataLog_Normal <- function(n,inter, a, k, meanX, varX){
  A <- rbinom(n, 1, 0.5)
  X <- matrix(rnorm(n * k, meanX, varX), nrow = n, ncol = k)
  colnames(X) <- paste0("x", 1:k)
  b <- c(rep(1, k))
  linearY <- inter + a * A + (X %*% b) + rnorm(n, 0, 0.5)
  Y <- rbinom(n, 1, plogis(linearY))
  return(data.frame(Y = Y, A, X))
}

simulateDataLog_Normal_Noise <- function(n,inter, a, k, meanX, varX){
  A <- rbinom(n, 1, 0.5)
  X <- matrix(rnorm(n * k, meanX, varX), nrow = n, ncol = k)
  U <- matrix(rnorm(n *2 *  k, meanX, varX), nrow = n, ncol = 2*k)
  colnames(U) <- paste0("u", 1:(2*k))
  colnames(X) <- paste0("x", 1:k)
  b <- c(rep(1, k))
  linearY <- inter + a * A + (X %*% b) + rnorm(n, 0, 0.5)
  Y <- rbinom(n, 1, plogis(linearY))
  return(data.frame(Y = Y, A, X, U))
}

#simulateDataLog_Exp <- function(n,inter, a, k, meanX, varX){
#  A <- rbinom(n, 1, 0.5)
#  X <- matrix(rexp(n * k), nrow = n, ncol = k)
#  U <- matrix(rexp(n *2 *  k), nrow = n, ncol = 2*k)
#  colnames(U) <- paste0("u", 1:(2*k))
# colnames(X) <- paste0("x", 1:k)
# b <- c(rep(1, k))
# linearY <- inter + a * A + (X %*% b) + rnorm(n, 0, 0.5)
# Y <- rbinom(n, 1, plogis(linearY))
# return(data.frame(Y = Y, A, X, U))
#}

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
  p.fac[which(colnames(x) == "A")] <- 0
  
  
  
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


two_prop_z <- function(data){
  n1 <- sum(data$A == 1)
  n0 <- sum(data$A == 0)
  
  y1 <- sum(data$Y[data$A == 1] == 1)
  y0 <- sum(data$Y[data$A == 0] == 1)
  return(prop.test(x = c(y1, y0), n = c(n1, n0), correct = FALSE)$p.value)
}


# Simulation:



b0<-0
b1 <- 0
data <- simulateDataLog_Normal(n, b0, b1, 3, 0, 2)

dual_permutation_test_logit(data, permutations, "perm", "BIC")

two_prop_z(data)

testing_Logit(data_logit1, permutations, "sem")



testing_Logit(data_logit1, permutations, "rand")








# ----------------------------------- Logit parallell of power for Rand. Norm with noise --------------------------------------------


# permutations* n.seed * 2/10.000 min = time

#Run at bet1 = 0 for type 1 error rate


# ------ Wrapper for parallel ------
bet0<-0.1
bet1<- 5
alpha <- 0.05
permutations <- 10000
type <- "BIC"


tester.seed_logit <- function(n, permut, type){
  LogRand <- simulateDataLog_Normal(n, bet0, bet1, k=3, meanX = 0, varX = 2)
  
  
  perm <- dual_permutation_test_logit(LogRand, permutations, "perm", type)
  sem <- dual_permutation_test_logit(LogRand, permutations, "sem", type)
  rand <- dual_permutation_test_logit(LogRand,permutations, "rand", type)
  
  
  p<- t(matrix(c(perm$p_lasso, sem$p_lasso, rand$p_lasso, perm$p_post_lasso,
                 sem$p_post_lasso, rand$p_post_lasso, two_prop_z(LogRand))))
  return(data.frame(p))
}

# ------ Parallel setup ------
n.seed <- 1000

params <- expand_grid(
  seed = 1:n.seed
)
param_list <- split(params, seq_len(nrow(params)))

n_cores <- max(1, detectCores() - 2)
cl <- makeCluster(n_cores)


clusterSetRNGStream(cl)

clusterEvalQ(cl, {
  library(glmnet)
  library(tidyverse)
})


clusterExport(
  cl,
  c("bet0","bet1", "alpha", "tester.seed_logit", "simulateDataLog_Normal", "dual_permutation_test_logit"
    , "permutations", "n", "ATE_Calculator_logit", "type", "two_prop_z")
)

# ----- Run simulation -----

results_Logit_rand <- pblapply(param_list, cl = cl, FUN = function(param) {
  tester.seed_logit(n = n, permut =  permutations, type = type)
})

stopCluster(cl)

results.sim <- data.frame(do.call(rbind, results_Logit_rand))
colnames(results.sim) <- c("Perm_lasso", "Sem_Perm_lasso", "Rand_lasso" ,
                           "Perm_post_lasso", "Sem_Perm_post_lasso", "Rand_post_lasso", "two_prop_z")

p_mean <- c(colMeans(results.sim[,], na.rm = TRUE))
power <- colMeans(results.sim[, c("Perm_lasso", "Sem_Perm_lasso", "Rand_lasso" ,
                                  "Perm_post_lasso", "Sem_Perm_post_lasso", "Rand_post_lasso", "two_prop_z")]
                  <= alpha)

resultsRandNorm <- if(bet1 ==0 ) data.frame( p_mean, type1 = power) else data.frame(p_mean, power)

results.sim <- data.frame(results.sim[, 0:6])
results.sim$seed <- 1:n.seed
rownames(results.sim) <- NULL

print(resultsRandNorm)


