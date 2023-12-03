library(tidyverse)

data = read_csv('brandchoicesDataSet.csv')
beta = rep(1, 4)

p_t <- function(beta,t){
  x = data %>%
    filter(SET == t) %>%
    select(BRAND, FASH, QUAL, PRICE) %>%
    distinct() %>%
    select(-BRAND) %>%
    mutate(c = 1) %>%
    as.matrix()
  
  return(tibble(SET = t,
                p_t = 1 + sum(exp(x %*% beta))))
  
}

# tentativa falha

probs_llik <- function(beta){
  x <- data %>%
    select(BRAND, SET, FASH, QUAL, PRICE) %>%
    distinct() %>%
    select(-c(BRAND, SET)) %>%
    mutate(c = 1) %>%
    as.matrix()
  
  denom <- 1 + sum(exp(x %*% beta))
  
  probs <- data %>%
    mutate(p_tk = exp(FASH * beta[1] + QUAL * beta[2] + PRICE * beta[3] + beta[4])/denom)
  
  llik <- probs %>%
    mutate(llik_itk = CHOICE * log(p_tk)) %>%
    summarise(sum(llik_itk)) %>%
    as.numeric()
  return(list(probabilities = probs, llik = llik))
}

probs_llik <- function(beta){
  denom <- data %>%
    select(SET) %>%
    distinct() %>%
    {.$SET} %>%
    map_df(~ p_t(beta, .x))
  
  probs <- data %>%
    left_join(denom) %>%
    mutate(p_tk = ifelse(BRAND == 0,
                         1,
                         exp(FASH * beta[1] + QUAL * beta[2] + PRICE * beta[3] + beta[4]))/p_t)
  
   llik <- probs %>%
    mutate(llik_itk = CHOICE * log(p_tk)) %>%
    summarise(sum(llik_itk)) %>%
    as.numeric()
    return(list(probabilities = probs, llik = llik))
}

get_llik <- function(beta){
  return(-probs_llik(beta)$llik)
}

result <- optim(par = c(0,0,0,1),
                fn = get_llik,
                method= 'BFGS')

beta_ml <- result$par

probs_llik(beta_ml)$probabilities %>%
  select(BRAND, SET, p_tk) %>%
  arrange(SET, BRAND) %>%
  distinct() %>%
  pivot_wider(names_from = BRAND, values_from = p_tk)

data %>%
  group_by(SET, BRAND) %>%
  summarise(pbar = sum(CHOICE)/n()) %>%
  ungroup() %>%
  pivot_wider(names_from = BRAND, values_from = pbar)

# Item c

# Score function given beta, t, k and j
score_fun <- function(beta, t, k, j){
  x = data %>%
    filter(SET == t) %>%
    select(BRAND, FASH, QUAL, PRICE) %>%
    arrange(BRAND) %>%
    distinct() %>%
    select(-BRAND) %>%
    mutate(c = 1) %>%
    as.matrix()
  
  pt0 <- (1 + sum(exp(x %*% beta)))^(-1)
  
  if(k != 0){
    return(x[k+1,j] - exp(x[k+1, ] %*% beta) * score_fun(beta, t, 0, j))
  }else{
    return(- pt0 * (x[,j] %*% exp(x %*% beta)))
  }
    
}

# The equivalent for hessian
sscore_fun <- function(beta, t, k, j, l){
  x = data %>%
    filter(SET == t) %>%
    select(BRAND, FASH, QUAL, PRICE) %>%
    arrange(BRAND) %>%
    distinct() %>%
    select(-BRAND) %>%
    mutate(c = 1) %>%
    as.matrix()
  
  pt0 <- (1 + sum(exp(x %*% beta)))^(-1)
  
  if(k == 0){
    return(-pt0* (score_fun(beta, t, 0, l) * (x[,j] %*% exp(x %*% beta)) + ((x[,j] * x[,l]) %*% exp(x %*% beta))))
  }else{
    return(sscore_fun(beta, t, 0, j, l))
  }
}

# Score for each j
nabla <- function(beta, t, k){
  return(map(c(1:4), ~ score_fun(beta, t = t, k = k, j = .x)) %>%
           unlist)
}

# ... and the hessian
inner_fisher <- function(beta, t, k){
  H <- matrix(nrow = 4, ncol = 4)
  for(i in 1:4){
    for(j in 1:4){
      H[i,j] <- sscore_fun(beta, t, k, i, j)
    }
  }
  return(H)
}

# score score'
meat <- function(k,t){
  return(tibble(BRAND = k, SET = t) %>%
           mutate(scsc = list(nabla(beta_ml, t, k) %*% t(nabla(beta_ml, t, k)))))
}

# A function that stores the hessian for a pair (t,k) in a tibble structure
H_i <- function(k,t){
  return(tibble(BRAND = k, SET = t) %>%
           mutate(H = list(inner_fisher(beta_ml, t, k))))
}

# Now a tibble with matrices in its entries
H_kt <- pmap(
  data.frame(k = data$BRAND,
             t = data$SET) %>%
    distinct(),
  ~ H_i(k = ..1, t = ..2)
) %>%
  do.call(rbind, .)

# Fisher
J <- data %>%
  left_join(H_kt) %>%
  filter(CHOICE == 1) %>%
  {.$H} %>%
  Reduce('+', .) * (-1)

# Corretly specified estimate
Jinv <- solve(Ical)

# ... or just the slice of bread

# Analogously, a tibble with "meat" matrices
meat_kt <- pmap(
  data.frame(k = data$BRAND,
             t = data$SET) %>%
    distinct(),
  ~ meat(k = ..1, t = ..2)
) %>%
  do.call(rbind, .)

# The meat, or Score Score'
Meat <- data %>%
  left_join(meat_kt) %>%
  filter(CHOICE == 1) %>%
  {.$scsc} %>%
  Reduce('+', .)

# Misspecified 
Hinv %*% Meat %*% Hinv 
