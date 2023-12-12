library(tidyverse)
library(xtable)
library(ggcorrplot)
library(Hotelling)

source('resources.R')

data = read_csv('brandchoicesDataSet.csv')

# Question 1, item b

beta0 = rep(1, 4)/10

probs_llik <- function(beta){
  x = data %>%
    select(FASH, QUAL, PRICE) %>%
    mutate(c = 1) %>%
    as.matrix() %>%
    {exp(. %*% beta)} %>%
    `colnames<-`('p_tk')
  
  probs <- data %>%
    cbind(x) %>%
    group_by(SET, id) %>%
    mutate(p_tk = ifelse(BRAND == 0,
                         1,
                         p_tk)/(1+sum(p_tk))) %>%
    ungroup()
  
  llik <- probs %>%
    mutate(llik_itk = CHOICE * log(p_tk)) %>%
    dplyr::summarise(sum(llik_itk)) %>%
    as.numeric()
  return(list(probabilities = probs, llik = llik))
}

get_llik <- function(beta){
  return(-probs_llik(beta)$llik)
}

result <- optim(par = beta0,
                fn = get_llik,
                method= 'BFGS')

beta_ml <- result$par

optim_to_latex(result, label = "tab:mv",
               caption = 'Logit - Maximum Likelihood estimate',
               param_labels = c(paste0('$\\beta_', c(1:3,0),'$'))) %>%
  writeLines(., 'tables/mv.tex')

probs_llik(beta_ml)$probabilities %>%
  select(BRAND, SET, p_tk) %>%
  arrange(SET, BRAND) %>%
  distinct() %>%
  pivot_wider(names_from = BRAND, values_from = p_tk) %>%
  select(-SET) %>%
  as.matrix() %>%
  ggcorrplot(lab = T, lab_size = 8,
             colors = c(my_col[4], 'white', my_col[1])) +
  labs(x = 'Set', y = 'Brand', fill = '') +
  my_theme +
  guides(fill = 'none') +
  scale_x_continuous(breaks = c(1:8))
ggsave('figures/prob-mv.pdf', scale = 1.5, device = cairo_pdf)

pbar <- data %>%
  group_by(SET, BRAND) %>%
  dplyr::summarise(pbar = sum(CHOICE)/n()) %>%
  ungroup()

pbar %>%
  pivot_wider(names_from = BRAND, values_from = pbar) %>%
  select(-SET) %>%
  as.matrix() %>%
  ggcorrplot(lab = T, lab_size = 8,
             colors = c(my_col[4], 'white', my_col[1])) +
  labs(x = 'Set', y = 'Brand', fill = '') +
  my_theme +
  guides(fill = 'none') +
  scale_x_continuous(breaks = c(1:8))
ggsave('figures/prob-avg.pdf', scale = 1.5, device = cairo_pdf)

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
    return(-pt0* (score_fun(beta, t, 0, l) * (x[,j] %*% exp(x %*% beta)) +
                    ((x[,j] * x[,l]) %*% exp(x %*% beta))))
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
Jinv <- solve(J)

sqrt(diag(Jinv)) %>%
  t() %>%
  round(3) %>%
  `colnames<-`(c(paste0('$\\beta_', c(1:3,0),'$'))) %>%
  {print(xtable(., label = 'tab:normalAvar',
                caption = 'Estimated standard errors - Inverse Fisher'),
         type='latex', sanitize.text.function=identity,
         include.rownames = FALSE,
        file = "tables/normal.tex")}

# ... or just a slice of bread

# Analogously, a tibble with "meat" matrices
meat_kt <- pmap(
  data.frame(k = data$BRAND,
             t = data$SET) %>%
    distinct(),
  ~ meat(k = ..1, t = ..2)
) %>%
  do.call(rbind, .)

# The bread, or Score Score'
Meat <- data %>%
  left_join(meat_kt) %>%
  filter(CHOICE == 1) %>%
  {.$scsc} %>%
  Reduce('+', .)/400

# Misspecified 
sand <- Jinv %*% Meat %*% Jinv

sqrt(diag(sand)) %>%
  t() %>%
  round(3) %>%
  `colnames<-`(c(paste0('$\\beta_', c(1:3,0),'$'))) %>%
  {print(xtable(., label = 'tab:sandwich',
                caption = 'Estimated standard errors - Sandwich Estimator'),
         type='latex', sanitize.text.function=identity,
         include.rownames = FALSE,
         file = "tables/sandwich.tex")}

# Item d

probs_llik(beta_ml)$probabilities %>%
  select(BRAND, SET, p_tk) %>%
  arrange(SET, BRAND) %>%
  distinct() %>%
  group_by(SET) %>%
  mutate(pt1 = ifelse(BRAND == 1,
                      p_tk,
                      0)) %>%
  mutate(pt1 = sum(pt1)) %>%
  ungroup() %>%
  mutate(partial = ifelse(BRAND == 1, p_tk*(1-p_tk), -p_tk * pt1) *beta_ml[3]) %>%
  select(SET, BRAND, partial) %>%
  pivot_wider(names_from = BRAND, values_from = partial) %>%
  select(-SET) %>%
  as.matrix() %>%
  ggcorrplot(lab = T, lab_size = 8,
             colors = c(my_col[4], 'white', my_col[1])) +
  labs(x = 'Set', y = 'Brand', fill = '') +
  my_theme +
  guides(fill = 'none') +
  scale_x_continuous(breaks = c(1:8))
ggsave('figures/mg1d.pdf', scale = 1.5, device = cairo_pdf)

# Item e

data <- data %>%
  filter(BRAND != 0)

beta0_IIA = rep(1, 3)/10

probs_llik_IIA <- function(beta){
  x = data %>%
    select(FASH, QUAL, PRICE) %>%
    #mutate(c = 1) %>%
    as.matrix() %>%
    {exp(. %*% beta)} %>%
    `colnames<-`('p_tk')
  
  probs <- data %>%
    cbind(x) %>%
    group_by(SET, id) %>%
    mutate(p_tk = ifelse(BRAND == 0,
                         1,
                         p_tk)/(1+sum(p_tk))) %>%
    ungroup()
  
  llik <- probs %>%
    mutate(llik_itk = CHOICE * log(p_tk)) %>%
    dplyr::summarise(sum(llik_itk)) %>%
    as.numeric()
  return(list(probabilities = probs, llik = llik))
}

get_llik_IIA <- function(beta){
  return(-probs_llik_IIA(beta)$llik)
}

result_IIA <- optim(par = beta0_IIA,
                fn = get_llik_IIA,
                method= 'BFGS')

result_IIA_orig <- result_IIA

result_IIA$par <- c(result_IIA$par, 0)

optim_to_latex(result_IIA, label = "tab:mvIIA",
               caption = 'Logit - Independence of Irrelevant Alternative - Maximum Likelihood estimate',
               param_labels = c(paste0('$\\beta_', c(1:3,0),'$')))  %>%
  writeLines(., 'tables/logit_iia.tex')

result_IIA$par
beta_ml

# Question 2, item b
data = read_csv('brandchoicesDataSet.csv')

probs_llik_nested <- function(beta){
  x2 <- data %>%
    select(BRAND, SET, FASH, QUAL, PRICE) %>%
    filter(BRAND != 0) %>%
    distinct() %>%
    select(-c(BRAND, SET)) %>%
    mutate(c = 1) %>%
    as.matrix()
  
  x1 <- data %>%
    select(BRAND, SET, FASH, QUAL, PRICE) %>%
    filter(BRAND == 0) %>%
    distinct() %>%
    select(-c(BRAND, SET)) %>%
    mutate(c = 1) %>%
    as.matrix()
  
  denom <- sum(exp(x1 %*% beta[-5])) + (sum(exp(x2 %*% beta[-5] / beta[5])))^beta[5]
  
  probs <- data %>%
    mutate(p_tk = ifelse(BRAND == 0,
                         sum(exp(x1 %*% beta[-5])),
                         sum(exp(x2 %*% beta[-5] / beta[5]))^beta[5] *
                           exp((FASH * beta[1] + QUAL * beta[2] + PRICE * beta[3] + beta[4])/beta[5]) /
                           (sum(exp(x2 %*% beta[-5] / beta[5]))))/denom)
  
  llik <- probs %>%
    mutate(llik_itk = CHOICE * log(p_tk)) %>%
    dplyr::summarise(sum(llik_itk)) %>%
    as.numeric()
  return(list(probabilities = probs, llik = llik))
}

get_llik_nested <- function(beta){
  return(-probs_llik_nested(beta)$llik)
}

beta0_nested <- rep(1,5)
result_nested <- optim(par = beta0_nested,
                fn = get_llik_nested,
                method= 'BFGS')
beta_nested <-result_nested$par

optim_to_latex(result_nested, label = "tab:nestmv",
               caption = 'Nested Logit - Maximum Likelihood estimate',
               param_labels = c(paste0('$\\beta_', c(1:3,0),'$'), '$\\sigma$')) %>%
  writeLines(., 'tables/nestmv.tex')

# Item c

probs_llik_nested(beta_nested)$probabilities %>%
  select(BRAND, SET, p_tk) %>%
  arrange(SET, BRAND) %>%
  distinct() %>%
  filter(BRAND == 0) %>% # k = 0
  mutate(partial = - p_tk * (1-p_tk) * beta_nested[3]) %>%
  select(partial)  %>%
  as.matrix() %>% # Partial effects for each t
  `colnames<-`('0') %>%
  `rownames<-`(c(1:8)) %>%
  t() %>%
  `rownames<-`(c('Effect')) %>%
  {print(xtable(., label = 'tab:mg2c',
                caption = 'Marginal price increase - brand 1 on choice probability of brand 0'),
         type='latex', sanitize.text.function=identity,
         include.rownames = TRUE,
         file = "tables/Partial2c.tex")}
  
# Item d

min_dist <- function(beta){
  diag_pbar <- pbar %>%
    filter(BRAND != 0) %>%
    arrange(BRAND, SET) %>%
    {.$pbar %>% diag()}
  
  estim_prob <- probs_llik_nested(beta)$probabilities %>%
    select(BRAND, SET, p_tk) %>%
    filter(BRAND != 0) %>%
    arrange(BRAND, SET) %>%
    distinct() %>%
    {.$p_tk %>% diag}
  
  diff_prob <- (estim_prob - diag_pbar) %*% 
    (diag(1, nrow = 8) %>%
       rbind(.,.,.))
  
  map(1:3, ~
        data %>%
        select(BRAND, SET, FASH, QUAL, PRICE) %>%
        filter(BRAND == .x) %>%
        arrange(BRAND, SET) %>%
        distinct() %>%
        mutate(c = 1) %>%
        select(c, FASH, QUAL, PRICE) %>%
        as.matrix() %>% t() %>%
        `colnames<-`(1:8) %>%
        {1/8 * . %*% diff_prob[((.x-1) * 8 + 1):(.x * 8),]}) %>%
    do.call(rbind,.) %>%
    rowSums() %>%
    {t(.) %*% .} %>%
    as.numeric() %>%
    return()
  
}

result_gmm <- optim(par = beta0_nested,
                       fn = min_dist,
                       method= 'BFGS')
result_gmm

optim_to_latex(result_gmm, label = "tab:gmm",
               caption = 'Nested Logit - GMM estimate',
               param_labels = c(paste0('$\\beta_', c(1:3,0),'$'), '$\\sigma$')) %>%
  writeLines(., 'tables/nestgmm.tex')

# Question 3, item c/d

S = 1000

set.seed(1981)
eta <- rnorm(S)

data_expand <- data %>%
  #select(BRAND, FASH, QUAL, PRICE) %>%
  #distinct() %>%
  slice(rep(1:n(), each = 1000)) %>%
  mutate(Sim = rep(1:S, nrow(data)))

probs_llik_sim <- function(beta, eta){
  x = data_expand %>%
    select(FASH, QUAL, PRICE) %>%
    as.matrix() %>%
    {. %*% beta[-c(4:5)]} %>%
    `colnames<-`('p_tk')
  
  # Adding the "intercept"
  x <- exp(x + beta[4] + rep(eta, nrow(data)) * beta[5])
  
  probs <- data_expand %>%
    cbind(x) %>%
    group_by(SET, id, Sim) %>%
    mutate(p_tk = p_tk/(1+sum(p_tk))) %>%
    ungroup()
  
  probs_avg <- probs %>%
    group_by(id, BRAND, SET, CHOICE, FASH, QUAL, PRICE) %>%
    dplyr::summarise(Lis = mean(p_tk, na.rm = T)) %>%
    ungroup()
  
  llik <- probs_avg %>%
    mutate(llik_itk = CHOICE * log(Lis)) %>%
    dplyr::summarise(sum(llik_itk, na.rm = T)) %>%
    as.numeric()
  return(list(probabilities = probs_avg, llik = llik))
}

get_llik_sim <- function(beta){
  return(-probs_llik_sim(beta, eta)$llik)
}

result_mix_mv <- optim(par = beta0_nested,
                fn = get_llik_sim,
                method= 'BFGS')

optim_to_latex(result_mix_mv, label = "tab:mixed-mv",
               caption = 'Mixed Logit - Maximum Likelihood estimate',
               param_labels = c(paste0('$\\beta_', c(1:3,0),'$'), '$\\omega$')) %>%
  writeLines(., 'tables/mix-mv.tex')

# Item e

min_dist_mix <- function(beta){
  diag_pbar <- pbar %>%
    filter(BRAND != 0) %>%
    arrange(BRAND, SET) %>%
    {.$pbar %>% diag()}
  
  estim_prob <- probs_llik_sim(beta, eta)$probabilities %>%
    select(BRAND, SET, Lis) %>%
    filter(BRAND != 0) %>%
    arrange(BRAND, SET) %>%
    distinct() %>%
    {.$Lis %>% diag}
  
  diff_prob <- (estim_prob - diag_pbar) %*% 
    (diag(1, nrow = 8) %>%
       rbind(.,.,.))
  
  map(1:3, ~
        data %>%
        select(BRAND, SET, FASH, QUAL, PRICE) %>%
        filter(BRAND == .x) %>%
        arrange(BRAND, SET) %>%
        distinct() %>%
        mutate(c = 1) %>%
        select(c, FASH, QUAL, PRICE) %>%
        as.matrix() %>% t() %>%
        `colnames<-`(1:8) %>%
        {1/8 * . %*% diff_prob[((.x-1) * 8 + 1):(.x * 8),]}) %>%
    do.call(rbind,.) %>%
    rowSums() %>%
    {t(.) %*% .} %>%
    as.numeric() %>%
    return()
  
}

result_mix_gmm <- optim(par = beta0_nested,
                    fn = min_dist_mix,
                    method= 'BFGS')
result_mix_gmm

optim_to_latex(result_mix_gmm, label = "tab:mixed-gmm",
               caption = 'Mixed Logit - GMM estimate',
               param_labels = c(paste0('$\\beta_', c(1:3,0),'$'), '$\\omega$')) %>%
  writeLines(., 'tables/mix-gmm.tex')

# Question 4
orig_data <- data

B = 1000
Nb <- n_distinct(orig_data$id)
set.seed(1981)

ibs <- list()
betaB <- list()
betaB_IIA <- list()
samples <- list()

for(b in 1:B){
  print(paste0('b = ', b))
  ibs[[b]] <- sample(1:Nb, size = Nb, replace = T) # Item a
  
  # Item b
  
  # Subset the data frame based on the sampled IDs
  subset_df <- orig_data[orig_data$id %in% ibs[[b]], ]
  
  # Identify the counts of each ID in the sampled_ids vector
  id_counts <- table(ibs[[b]])
  
  # Create an empty data frame to store the final subset with duplicated entries
  data <- data.frame()
  
  # Loop through each unique ID in the sampled IDs vector
  for (idd in unique(ibs[[b]])) {
    # Subset the rows of subset_df that match the current ID
    subset_for_id <- subset_df[subset_df$id == idd, ]
    
    # Duplicate the rows based on the count of the current ID
    if (id_counts[as.character(idd)] > 1) {
      duplicated_entries <- subset_for_id[rep(1:nrow(subset_for_id),
                                              each = id_counts[as.character(idd)]), ]
      data <- rbind(data, duplicated_entries)
    } else {
      data <- rbind(data, subset_for_id)
    }
  }
  
  samples[[b]] <- data
  
  
  betaB[[b]] <- optim(par = beta0,
                 fn = get_llik,
                 method= 'BFGS')$par
  
  data <- data %>%
    filter(BRAND != 0)
  
  betaB_IIA[[b]] <- optim(par = beta0_IIA,
                     fn = get_llik_IIA,
                     method= 'BFGS')$par
}

betaB_IIA_tibble <- do.call(cbind, betaB_IIA) %>%
  t() %>%
  `colnames<-`(paste0('beta', c(1:3))) %>%
  data.frame() %>%
  tibble()

betaB_tibble <- do.call(cbind, betaB) %>%
  t() %>%
  `colnames<-`(paste0('beta', c(1:3,0))) %>%
  data.frame() %>%
  tibble()

hottest <- hotelling.test(betaB_tibble %>%
                            select(-beta0),
                         betaB_IIA_tibble)

hottest$stats

paste0("\\begin{table}[h]\n",
       "\\centering\n",
       "\\caption{Hotelling's $T^2$ test}\\label{tab:hotelling}\n",
       "\\begin{tabular}{lc}\n \\hline \n",
       'Test stat &', round(hottest$stats$statistic,3), '\\\\',
       '\nNumerator df & ', round(hottest$stats$df[1],3), '\\\\',
       '\nDenominator df & ', round(hottest$stats$df[2],3), '\\\\',
       '\nP-value & ', round(hottest$pval,3), '\\\\',
       "\\hline\\\\\n", "\\end{tabular}\n",
       "\\end{table}\n") %>%
  writeLines(., 'tables/hotelling.tex')


betaB_tibble %>%
  mutate(IIA = 0) %>%
  bind_rows(betaB_IIA_tibble %>%
              mutate(IIA = 1, beta0 = 0)) %>%
  mutate(IIA = factor(IIA)) %>%
  pivot_longer(beta1:beta0) %>%
  ggplot(aes(y = value, x = IIA)) + geom_boxplot() +
  facet_wrap(~name, scales = 'free') + my_theme
ggsave('figures/boxplotbetas.pdf', scale = 1.5, device = cairo_pdf)

# Item c

#i.

truePbar <- orig_data %>%
  group_by(SET, BRAND) %>%
  dplyr::summarise(pbar = sum(CHOICE)/n()) %>%
  ungroup() %>%
  filter(SET == 1 & BRAND == 1) %>%
  {.$pbar}

Pb <- map_df(samples,
             function(x) x %>%
               filter(SET == 1 & BRAND == 1) %>%
               dplyr::summarise(pbar = sum(CHOICE)/Nb)
               )

mean(Pb$pbar) # Estim
truePbar # "True"
var(Pb$pbar) # Estim
truePbar*(1-truePbar)/Nb # "True"

paste0("\\begin{table}[h!]\n",
       "\\centering\n",
       "\\caption{$\\bar{P}$ and $\\bar{P}^b$ mean and standard error}\\label{tab:pb-comp}\n",
       "\\begin{tabular}{lcc}\n \\hline \n",
       " & $\\bar{P}$ & $\\bar{P}^b$ \\\\ \n", '\\hline',
       '\nMean &', round(truePbar,3), " & ", round(mean(Pb$pbar),3), '\\\\',
       '\nStd. Dev. & ', round(sqrt(truePbar*(1-truePbar)/Nb),3), "&", round(sd(Pb$pbar),3), '\\\\',
       "\\hline\\\\\n", "\\end{tabular}\n",
       "\\end{table}\n") %>%
  writeLines(., 'tables/mean_var.tex')

# ii.
# alpha
set.seed(1981)
W <- map(1:B, function(.x) rexp(Nb, rate = 1))

# ... depois a gente termina

# beta

Pbtilde <- map_df(1:B,
                  function(x) samples[[x]] %>%
                    filter(SET == 1 & BRAND == 1) %>%
                    mutate(wCHOICE = W[[x]] * CHOICE) %>%
                    dplyr::summarise(pbar = sum(CHOICE)/sum(W[[x]])))

paste0("\\begin{table}[h!]\n",
       "\\centering\n",
       "\\caption{$\\bar{P}^b$ and $\\tilde{P}^b$ mean and standard error}\\label{tab:pb-comp}\n",
       "\\begin{tabular}{lcc}\n \\hline \n",
       " & $\\bar{P}^b$ & $\\tilde{P}^b$ \\\\ \n", '\\hline',
       '\nMean &', round(mean(Pb$pbar),3), " & ", round(mean(Pbtilde$pbar),3), '\\\\',
       '\nStd. Dev. & ', round(sd(Pb$pbar),3), "&", round(sd(Pbtilde$pbar),3), '\\\\',
       "\\hline\\\\\n", "\\end{tabular}\n",
       "\\end{table}\n") %>%
  writeLines(., 'tables/mean_var_tilde.tex')
