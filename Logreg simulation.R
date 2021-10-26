library(mice)
library(tidyverse)
library(data.table)

#################################

# IMPORT DATA

url <- "http://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data"

data <- read.csv(url, header=FALSE)

colnames(data) <- c(
  "age",
  "sex",# 0 = female, 1 = male
  "cp", # chest pain
  # 1 = typical angina,
  # 2 = atypical angina,
  # 3 = non-anginal pain,
  # 4 = asymptomatic
  "trestbps", # resting blood pressure (in mm Hg)
  "chol", # serum cholestoral in mg/dl
  "fbs",  # fasting blood sugar if less than 120 mg/dl, 1 = TRUE, 0 = FALSE
  "restecg", # resting electrocardiographic results
  # 1 = normal
  # 2 = having ST-T wave abnormality
  # 3 = showing probable or definite left ventricular hypertrophy
  "thalach", # maximum heart rate achieved
  "exang",   # exercise induced angina, 1 = yes, 0 = no
  "oldpeak", # ST depression induced by exercise relative to rest
  "slope", # the slope of the peak exercise ST segment
  # 1 = upsloping
  # 2 = flat
  # 3 = downsloping
  "ca", # number of major vessels (0-3) colored by fluoroscopy
  "thal", # this is short of thalium heart scan
  # 3 = normal (no cold spots)
  # 6 = fixed defect (cold spots during rest and exercise)
  # 7 = reversible defect (when cold spots only appear during exercise)
  "hd" # (the predicted attribute) - diagnosis of heart disease
  # 0 if less than or equal to 50% diameter narrowing
  # 1 if greater than 50% diameter narrowing
)

data[data == "?"] <- NA

data$hd <- ifelse(test=data$hd == 0, yes=0, no=1)

data <- na.omit(data)

#################################

# SELECT RELEVANT FEATURES

selected_features <- c("sex", "trestbps", "chol", "thalach", "exang", "hd")
data <- data %>% select(all_of(selected_features))

#################################

# SPLIT DATA INTO TRAIN SET & TEST SET

train_fraction <- .7
sample_size <- floor(train_fraction * nrow(data))
train_ids <- sample(seq_len(nrow(data)), size=sample_size)
train_data <- data[train_ids, ]
test_data <- data[-train_ids, ]

#################################

# INDUCE MISSINGNESS

# Note:
# In both MAR mechanisms in the paper the threshold for the 2 numeric variables as mechanism was approx. mean + 15% SD, so we use a similar value.
# Missingness application order is also described in the paper.
#
# Similar to the paper, we also have columns that would make no sense to miss:
# - Age and sex make no sense to be missing
# - trestbps and chol should never be missing
#
# Explanation of the missingness mechanic structure used in the code:
# Each missingness mechanic is described by a set (list of variables in the mechanic, missingness fraction of affected variables, set of affected variables)
# The equal_or_greater variable in the case of a numeric mechanic indicates whether the threshold is to be read as '>= threshold' -> missingness or the other way around, and for a categoric variable whether the threshold is 'mech_var == specific value' -> missingness.

ignored_cols <- c('hd') # For example, outcome col should not be made missing 'slightlyMAR' and 'veryMAR' do not use these percentages, their missingness percentages per mechanic are described in the mechanics.
continuous_mar_sd_percent <- 0.15

# Do not touch these names
mech_names <- c('mech_var', 'is_numeric', 'equal_or_greater', 'point_value')
list_names <- c('mechanics', 'missingness_frac', 'affected_by_mech')


apply_mech <- function(dataset, mech_structure)
{
  miss_mechanics <- mech_structure[['mechanics']]
  missing_frac <- mech_structure[['missingness_frac']]
  affected_vars <- mech_structure[['affected_by_mech']]
  nr_vars <- length(affected_vars)
  
  for (row_id in seq_len(nrow(dataset)))
  {
    mechanic_failed <- FALSE
    
    for (mechanic in miss_mechanics)
    {
      if (mechanic[['is_numeric']])
      {
        if ((mechanic[['equal_or_greater']] & 
             dataset[row_id, mechanic[['mech_var']]] < mechanic[['point_value']]) |
            (!mechanic[['equal_or_greater']] &
             dataset[row_id, mechanic[['mech_var']]] > mechanic[['point_value']]))
        {
          # If this is true, a mechanic on a numeric value is not fulfilled.
          # Do not apply this mechanism to this row, and skip to next row.
          # Due to double for loop structure, must be done via a secondary variable.
          mechanic_failed <- TRUE
          break
        }
      }
      else # Mechanic on categorical variable
      {
        if ((mechanic[['equal_or_greater']] & 
             dataset[row_id, mechanic[['mech_var']]] != mechanic[['point_value']]) |
            (!mechanic[['equal_or_greater']] &
             dataset[row_id, mechanic[['mech_var']]] == mechanic[['point_value']]))
        {
          # If this is true, a mechanic on a categorical value is not fulfilled.
          # Do not apply this mechanism to this row, and skip to next row.
          # Due to double for loop structure, must be done via a secondary variable.
          mechanic_failed <- TRUE
          break
        }
      }
    }
    
    # If any mechanism failed, this mechanic should not be applied to this row.
    if (mechanic_failed)
    {
      next
    }
    
    probs <- runif(nr_vars)
    
    for (i in seq_len(nr_vars))
    {
      if (probs[i] <= missing_frac)
      {
        # For example, with frac 0.6 and runif sampling 0-1 you want 60%
        #   of sampled values to trigger this. So yes, <= is right.
        dataset[row_id, affected_vars[i]] = NA
      }
    }
  }
  
  return(dataset)
}


induce_slight_mar <- function(full_dataset)
{
  # Original paper: for 4/5 predictors, mechanism of 1 cont. value >= or <=
  #   threshold -> set prob missing (age skipped, unreasonable to be missing)
  # For categorical values: one state contributes to missing %, other is observed
  
  # As we have 13 instead of 5 predictors, we apply a similar mechanism
  #   with one column giving the threshold except this mechanism is applied
  #   to a lot more columns. In the paper one column was used often, we act similarly.
  
  # Note that, in the case of numerical variables, thresholds are calculated based
  #   on the dataset in this function. Train/test set information thus should not leak.
  
  # Construct all used lower and upper thresholds used
  thalach_lower  <- mean(full_dataset$thalach) - continuous_mar_sd_percent * sd(full_dataset$thalach)
  trestbps_upper <- mean(full_dataset$trestbps) + continuous_mar_sd_percent * sd(full_dataset$trestbps)
  
  thalach_mech  <- list('thalach', TRUE, FALSE, thalach_lower)
  trestbps_mech <- list('trestbps', TRUE, TRUE, trestbps_upper)
  
  names(thalach_mech)  <- mech_names
  names(trestbps_mech) <- mech_names
  
  sex_mech <- list('sex', FALSE, TRUE, 0)
  
  names(sex_mech) <- mech_names
  
  chol_set           <- list(list(trestbps_mech), 0.4, c('chol'))
  exang_trestbps_set <- list(list(thalach_mech), 0.4, c('exang', 'trestbps'))
  thalach_set        <- list(list(sex_mech), 0.3, c('thalach'))
  
  names(chol_set)           <- list_names
  names(exang_trestbps_set) <- list_names
  names(thalach_set)        <- list_names
  
  missing_dataset <- full_dataset
  
  missing_dataset <- apply_mech(missing_dataset, chol_set)
  missing_dataset <- apply_mech(missing_dataset, exang_trestbps_set)
  missing_dataset <- apply_mech(missing_dataset, thalach_set)
  
  return(missing_dataset)
}


induce_very_mar <- function(full_dataset)
{
  # Original paper: for 4/5 predictors, mechanism of up to 4 other variables
  #   fulfilling their respective conditions -> set prob missing (age skipped,
  #   unreasonable to be missing)
  
  # As we have 13 instead of 5 predictors, we again apply a similar mechanism.
  
  # Note that, in the case of numerical variables, thresholds are calculated based
  #   on the dataset in this function. Train/test set information thus should not leak.
  
  # Construct all used lower and upper thresholds used
  trestbps_lower <- mean(full_dataset$trestbps) - continuous_mar_sd_percent * sd(full_dataset$trestbps)
  trestbps_upper <- mean(full_dataset$trestbps) + continuous_mar_sd_percent * sd(full_dataset$trestbps)
  chol_upper     <- mean(full_dataset$chol) + continuous_mar_sd_percent * sd(full_dataset$chol)
  
  trestbps_leq_mech <- list('trestbps', TRUE, FALSE, trestbps_lower)
  trestbps_geq_mech <- list('trestbps', TRUE, TRUE, trestbps_upper)
  chol_geq_mech     <- list('chol', TRUE, TRUE, chol_upper)
  
  names(trestbps_leq_mech) <- mech_names
  names(trestbps_geq_mech) <- mech_names
  names(chol_geq_mech) <- mech_names
  
  sex_m_mech     <- list('sex', FALSE, TRUE, 1)
  sex_f_mech     <- list('sex', FALSE, TRUE, 0)
  exang_yes_mech <- list('exang', FALSE, TRUE, 1)
  exang_no_mech  <- list('exang', FALSE, TRUE, 0)
  
  names(sex_m_mech)     <- mech_names
  names(sex_f_mech)     <- mech_names
  names(exang_yes_mech) <- mech_names
  names(exang_no_mech)  <- mech_names
  
  thalach_set     <- list(list(sex_f_mech, chol_geq_mech, exang_no_mech),
                          0.6, c('thalach'))
  chol_set        <- list(list(exang_no_mech, trestbps_geq_mech),
                          0.6, c('chol'))
  exang_set       <- list(list(trestbps_leq_mech, sex_m_mech),
                          0.5, c('exang'))
  trestbps_set    <- list(list(sex_f_mech),
                          0.5, c('trestbps'))
  
  names(thalach_set)  <- list_names
  names(chol_set)     <- list_names
  names(exang_set)    <- list_names
  names(trestbps_set) <- list_names
  
  missing_dataset <- full_dataset
  
  missing_dataset <- apply_mech(missing_dataset, thalach_set)
  missing_dataset <- apply_mech(missing_dataset, chol_set)
  missing_dataset <- apply_mech(missing_dataset, exang_set)
  missing_dataset <- apply_mech(missing_dataset, trestbps_set)
  
  return(missing_dataset)
}


induce_mcar_v2 <- function(full_dataset, miss_frac)
{
  # Particular missing "completely" at random pattern:
  # [miss_frac] of rows are missing everything except outcome, the
  #   rest has everything fully observed. This results in the correct amount
  #   of missingness per indicator, and also in total (excl. outcome, that is)
  missing_pattern <- c(0,0,0,0,0,1)
  amputation <- ampute(full_dataset, prop=miss_frac, patterns=missing_pattern, mech='MCAR')
  
  return(amputation$amp)
}

#################################

# BASELINE

calc_accu <- function(model, test_data)
{
  test_data$model_prob <- predict(model, test_data, type = "response")
  test_data <- test_data %>% mutate(model_pred = 1 * (model_prob > .5))
  test_data <- test_data %>% mutate(accurate = 1 * (model_pred == hd))
  return(sum(test_data$accurate)/nrow(test_data))
}

results <- list()

# calculate mean & standard deviation for each continuous feature
data_mean <- sapply(data, function(x) ifelse(is.numeric(x), mean(x), NA))
data_sd <- sapply(data, function(x) ifelse(is.numeric(x), sd(x), NA))

# do logistic regression
baseline_logit <- glm(hd ~ sex + trestbps + chol + thalach + exang, data = train_data, family = "binomial")

# calculate model accuracy
accu <- calc_accu(baseline_logit, test_data)

# calculate 95% confidence interval
ci <- confint.default(baseline_logit)

# save results
key <- "baseline"

results[[key]]$accu <- accu
results[[key]]$params <- cbind(mean=shift(data_mean), sd=shift(data_sd), summary(baseline_logit)$coef[, 1:2], ci)

#################################

# SIMULATION

# Note:
# Each iteration induces missingness to the train and test sets in a number of configurations: MCAR-0.15, MCAR-0.25, MCAR-0.5, MCAR-0.7, sMAR, vMAR.
# For each missingness configuration, we impute the missing data with and without outcome variable y (hd) using MICE.
# For each imputed train data set, we fit a logistic regression model and use this model on test set to measure the accuracy. On top of that, we collect the estimate, standard error and 95% confidence interval for each coefficient from the model.

simulation_runs = 100

missing_config <- c("mcar-0.15", "mcar-0.25", "mcar-0.5", "mcar-0.7", "smar", "vmar")
imputation_config <- c("cc", "mi+y", "mi-y")

print("Simulation started")

for (i in seq(simulation_runs))
{
  for (j in missing_config)
  {
    
    # induce missingness
    if (j == "smar")
    {
      missing_train_data <- induce_slight_mar(train_data)
      missing_test_data <- induce_slight_mar(test_data)
    }
    else if (j == "vmar")
    {
      missing_train_data <- induce_very_mar(train_data)
      missing_test_data <- induce_very_mar(test_data)
    }
    else # if the current missing mechanism is MCAR
    {
      miss_frac <- as.numeric(unlist(regmatches(j, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", j)))) # extract miss_frac from the current missing_config
      
      missing_train_data <- induce_mcar_v2(train_data, miss_frac)
      missing_test_data <- induce_mcar_v2(test_data, miss_frac)
    }
    
    for (k in imputation_config)
    {
      # impute missing values
      if (k == "cc")
      {
        imputed_train_data <- missing_train_data[complete.cases(missing_train_data), ]
        imputed_test_data <- missing_test_data[complete.cases(missing_test_data), ]
      }
      else
      {
        if (k == "mi-y")
        {
          drop <- c("hd")
          
          train_data_y <- missing_train_data$hd
          missing_train_data <- missing_train_data[ , !(names(missing_train_data) %in% drop)]
          
          test_data_y <- missing_test_data$hd
          missing_test_data <- missing_test_data[ , !(names(missing_test_data) %in% drop)]
        }
        
        imputed_train_data <- mice(missing_train_data, maxit=20, m=1, method="pmm", remove.collinear = FALSE, printFlag = FALSE)
        imputed_test_data <- mice(missing_test_data, maxit=20, m=1, method="pmm", remove.collinear = FALSE, printFlag = FALSE)
        
        imputed_train_data <- complete(imputed_train_data)
        imputed_test_data <- complete(imputed_test_data)
        
        if (k == "mi-y")
        {
          imputed_train_data$hd <- train_data_y
          imputed_test_data$hd <- test_data_y
        }
      }
      
      # do logistic regression
      sim_logit <- glm(hd ~ sex + trestbps + chol + thalach + exang, data = imputed_train_data, family = "binomial")
      
      # calculate model accuracy
      accu <- calc_accu(sim_logit, imputed_test_data)
      
      # calculate coverage of 95% confidence interval
      ci <- confint.default(sim_logit)
      baseline_coef <- results[["baseline"]]$params[ , 3]
      ci_coverage <- numeric(6)
      
      for (l in 1:length(ci_coverage))
        ci_coverage[l] <- between(baseline_coef[l], ci[l, ][1], ci[l, ][2])
      
      # prepare results
      key <- paste(j, k, sep = "_")
      
      if (is.null(results[[key]]$accu))
        results[[key]]$accu <- 0
      
      if (is.null(results[[key]]$params))
        results[[key]]$params <- matrix(0, nrow = 6, ncol = 3)
      
      # save accuracy to results as a sum
      results[[key]]$accu <- results[[key]]$accu + accu
      
      # save coefficients, standard error & coverage of 95% ci to results as a summed matrix
      results[[key]]$params <- results[[key]]$params + cbind(summary(sim_logit)$coef[ , 1:2], Coverage=ci_coverage)
    }
  }
  
  # calculate mean and coefficient bias after each iteration for all cases except the baseline
  iteration_results <- lapply(results[c(-1)], function(x)
  {
    x$accu <- x$accu / i
    x$params <- x$params / i
    
    bias <- baseline_coef - x$params[ , 1]
    x$params <- cbind(x$params, Bias=bias)
    x$params <- x$params[ , c("Estimate", "Bias", "Std. Error", "Coverage")]
    
    return(x)
  })
  
  print(iteration_results)
  print(paste("Above are results for iteration: ", i))

  # Update the final simulation results with the last iteration_results
  if (i == simulation_runs)
    results <- c(results[1], iteration_results)
}