'''Algorithm: Two identical regresssion problems.
1. Learn order off first problem and apply to the other.
2. Learn order off second problem and apply to the other.
3. For each hypothesis, determine the max score from each group and record its label.
'''

library(mass)
library(caret)

create.regression.problem = function(m, n, seed = 1) {
  set.seed(seed)
  Sigma = matrix(c(1,0,0,1), nrow = 2)
  mu = c(0, 0)
  true.null.covariates = mvrnorm(n = m, mu = mu, Sigma = Sigma)
  mu.false = c(2, 2)
  false.null.covariates = mvrnorm(n = n, mu = mu.false, Sigma = Sigma)
  true.null.labels = expand.grid(replicate(m, c(-1, 1), simplify = FALSE))
  false.null.labels = matrix(1, ncol = n)
  results = list(true.null.labels, false.null.labels, true.null.covariates, false.null.covariates)
  names(results) = c('true.null.label', 'false.null.label', 'true.null.covariates', 'false.null.covariates')
  return(results)
}

sample_size = (2^m)

problem = create.regression.problem(m = 3, n = 2)

for (i in 1:sample_size) {
  true.null.labels.1 = problem$true.null.label[i, ]
  false.null.labels.1 = problem$false.null.label[i, ]
  true.null.covariates.1 = problem$true.null.covariates
  false.null.covariates.1 = problem$false.null.covariates
  fold1 = cbind(matrix(true.null.labels.1, ncol = 1), as.matrix(true.null.covariates.1))
  fold1 = as.data.frame(fold1)
  names(fold1) = c('Label', 'V1', 'V2')
  fold1$Label = as.factor(fold1$Label)
  for (j in 1:sample_size) {
    true.null.labels.2 = problem$true.null.label[j, ]
    false.null.labels.2 = problem$false.null.label[j, ]
    true.null.covariates.2 = problem$true.null.covariates
    false.null.covariates.2 = problem$false.null.covariates
    fold2 = cbind(matrix(true.null.labels.2, ncol = 1), as.matrix(true.null.covariates.2))
    fold2 = as.data.frame(fold1)
    names(fold2) = c('Label', 'V1', 'V2')
    fold2$Label = as.factor(fold2$Label)
    
    # Estimate preprocessing parameters
    preproc.param.1 = fold1 %>% 
      preProcess(method = c("center", "scale"))
    preproc.param.2 = fold2 %>% 
      preProcess(method = c("center", "scale"))
    # Transform the data using the estimated parameters
    fold1.transformed = preproc.param.1 %>% predict(fold1)
    fold2.transformed = preproc.param.2 %>% predict(fold2)
    
    
  }
}