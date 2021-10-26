# === Background information ===

# In both MAR mechanisms in the paper the threshold for the 2 numeric variables
#   as mechanism was approx. mean + 15% SD, so we use a similar value.
# Missingness application order is also described in the paper.

# Similar to the paper, we also have columns that would make no sense to miss:
# - Age and sex make no sense to be missing
# - trestbps and chol should never be missing

# Explanation of the missingness mechanic structure used in the code:
# Each missingness mechanic is described by a set (list of variables in the
#   mechanic, missingness fraction of affected variables, set of affected variables)
# The equal_or_greater variable in the case of a numeric mechanic indicates whether
#   the threshold is to be read as '>= threshold' -> missingness or the other way around,
#   and for a categoric variable whether the threshold is 'mech_var == specific value'
#   -> missingness.

library(tidyverse)
library(mice)

# Configurable paths and variables
in_folder <- 'C:/Users/20202596/Documents/R/data'
out_folder <- 'C:/Users/20202596/Documents/R/data/output'
main_filename <- 'processed-heart-disease.csv'
out_affix <- 'heartdisease.csv'

train_fraction <- 0.7
selection_cols <- c('exang', 'trestbps', 'chol',
                    'sex', 'thalach', 'hd')
ignored_cols <- c('hd') # For example, outcome col should not be made missing
# 'slightlyMAR' and 'veryMAR' do not use these percentages,
#   their missingness percentages per mechanic are described in the mechanics.
miss_percents <- c(0.15, 0.25, 0.5, 0.7) # 0.15 = original paper
miss_types <- c('MCAR')
continuous_mar_sd_percent <- 0.15

# Used for evaluation
writing_enabled = TRUE
simulation_wanted = TRUE
simulation_runs = 1000

# Previous process introduced a useless index col, remove it.
data <- read.csv(paste(in_folder, main_filename, sep='/'))
data <- data %>% select(-X) %>% select(all_of(selection_cols))
col_names <- colnames(data)

# Split into train & test sets
sample_size <- floor(train_fraction * nrow(data))
train_ids <- sample(seq_len(nrow(data)), size=sample_size)
train_data <- data[train_ids, ]
test_data <- data[-train_ids, ]

if (writing_enabled)
{
  write.csv(train_data, paste(out_folder, "truth_train.csv", sep='/'), row.names=FALSE)
  write.csv(test_data, paste(out_folder, "truth_test.csv", sep='/'), row.names=FALSE)
}

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



induce_mcar <- function(full_dataset, miss_frac)
{
  # Ampute does not let us ampute as much data in the way the original paper did,
  #   so implement MCAR amputation manually (error: "Proportion of missing cells
  #   is too large in combination with the desired number of missing variable)
  
  # For every column, produce a new set of ids for which rows will be made NA
  # to act in accordance with the original paper's "sample of about x% from
  #   each predictor was set to missing".
  missing_dataset <- full_dataset
  nr_missing <- floor(miss_frac * nrow(missing_dataset))
  
  for (i in seq(length(col_names)))
  {
    if (col_names[i] %in% ignored_cols)
    {
      # Obviously, 'ignored columns' should be ignored. Skip iteration.
      next
    }
    
    # For every column, sample and make a unique set of missing rows.
    miss_ids <- sample(seq_len(nrow(missing_dataset)), size=nr_missing)
    
    for (id in miss_ids)
    {
      missing_dataset[id, col_names[i]] = NA
    }
  }
  
  return(missing_dataset)
}



induce_mcar_v2 <- function(full_dataset, miss_frac)
{
  # Particular missing "completely" at random pattern:
  # [miss_frac] of rows are missing everything except outcome, the
  #   rest has everything fully observed. This results in the correct amount
  #   of missingness per indicator, and also in total (excl. outcome, that is)
  missing_pattern <- c(0,0,0,0,0,1)
  amputation <- ampute(full_dataset, prop=miss_frac,
                       patterns=missing_pattern, mech='MCAR')
  return(amputation$amp)
}



# Slightly MAR
# MICE's ampute() does not make possibly exactly what the paper did, sadly.
# Thus, we have implemented that concept ourselves.
# The variable selection and missingness patterns are covered in induce_slight_mar.
train_filename <- paste('slightlyMAR', 'train', out_affix, sep = '_')
test_filename <- paste('slightlyMAR', 'test', out_affix, sep = '_')
train_with_miss <- induce_slight_mar(train_data)
test_with_miss <- induce_slight_mar(test_data)
if (writing_enabled)
{
  write.csv(train_with_miss, paste(out_folder, train_filename, sep='/'), row.names=FALSE)
  write.csv(test_with_miss, paste(out_folder, test_filename, sep='/'), row.names=FALSE)
}



# Very MAR
# We keep the same amount of predictors needed for a mechanism to be active
#   (at most 4), but vary which columns these are more often as we have more.
train_filename <- paste('veryMAR', 'train', out_affix, sep = '_')
test_filename <- paste('veryMAR', 'test', out_affix, sep = '_')
train_with_miss <- induce_very_mar(train_data)
test_with_miss <- induce_very_mar(test_data)
if (writing_enabled)
{
  write.csv(train_with_miss, paste(out_folder, train_filename, sep='/'), row.names=FALSE)
  write.csv(test_with_miss, paste(out_folder, test_filename, sep='/'), row.names=FALSE)
}




# MCAR and optionally the non-paper MAR types
for (miss_frac in miss_percents)
{
  for (miss_type in miss_types)
  {
    train_filename <- paste(miss_type, miss_frac, 'train', out_affix, sep = '_')
    test_filename <- paste(miss_type, miss_frac, 'test', out_affix, sep = '_')
    
    # Generate both train and test missingness
    if (miss_type == 'MCAR')
    {
      #train_with_miss <- induce_mcar(train_data, miss_frac)
      #test_with_miss  <- induce_mcar(test_data, miss_frac)
      train_with_miss <- induce_mcar_v2(train_data, miss_frac)
      test_with_miss  <- induce_mcar_v2(test_data, miss_frac)
    }
    else
    {
      stop(paste("Missingness type requested but not implemented:", miss_type))
    }
      
    # Store both the create train and test miss sets
    if (writing_enabled)
    {
      write.csv(train_with_miss, paste(out_folder, train_filename, sep='/'), row.names=FALSE)
      write.csv(test_with_miss, paste(out_folder, test_filename, sep='/'), row.names=FALSE)
    }
  }
}

if (simulation_wanted)
{
  print("Initiating simulation runs, using full data (no train-test split here).")
  print(paste("Limiting MCAR tests to 15%, exact same function with different percentage",
              "should give the exact same SDs (and exp. mean == miss_frac obviously)"))
  
  slightMAR_missingness <- data.frame(matrix(nrow=simulation_runs, ncol=length(col_names)))
  veryMAR_missingness <- data.frame(matrix(nrow=simulation_runs, ncol=length(col_names)))
  MCAR_missingness <- data.frame(matrix(nrow=simulation_runs, ncol=length(col_names)))
  names(slightMAR_missingness) <- col_names
  names(veryMAR_missingness) <- col_names
  names(MCAR_missingness) <- col_names
  
  for(i in seq(simulation_runs))
  {
    if (i %% 100 == 0)
    {
      # Progress monitoring
      print(paste("Running simulation", i, "on all types"))
    }
    
    missing_data <- induce_slight_mar(data)
    missingness_percents <- map(missing_data, ~mean(is.na(.)))
    slightMAR_missingness[i,] <- missingness_percents
    
    missing_data <- induce_very_mar(data)
    missingness_percents <- map(missing_data, ~mean(is.na(.)))
    veryMAR_missingness[i,] <- missingness_percents
    
    missing_data <- induce_mcar_v2(data, 0.15)
    missingness_percents <- map(missing_data, ~mean(is.na(.)))
    MCAR_missingness[i,] <- missingness_percents
  }
  
  print("Simulations finished. Results:")
  
  simulation_summary_means <- data.frame(matrix(nrow=3, ncol=length(col_names)),
                                         row.names=c("Slight MAR mean % missing",
                                                     "Very MAR mean % missing",
                                                     "MCAR mean % missing"))
  names(simulation_summary_means) <- col_names
  simulation_summary_means[1,] = map(slightMAR_missingness, ~mean(.))
  simulation_summary_means[2,] = map(veryMAR_missingness, ~mean(.))
  simulation_summary_means[3,] = map(MCAR_missingness, ~mean(.))
  
  simulation_summary_sds <- data.frame(matrix(nrow=3, ncol=length(col_names)),
                                       row.names=c("Slight MAR % missing SD",
                                                   "Very MAR % missing SD",
                                                   "MCAR % missing SD"))
  names(simulation_summary_sds) <- col_names
  simulation_summary_sds[1,] = map(slightMAR_missingness, ~sd(.))
  simulation_summary_sds[2,] = map(veryMAR_missingness, ~sd(.))
  simulation_summary_sds[3,] = map(MCAR_missingness, ~sd(.))
  
  print("Missingness means:")
  print(simulation_summary_means)
  print("")
  print("Missingness SDs:")
  print(simulation_summary_sds)
}

