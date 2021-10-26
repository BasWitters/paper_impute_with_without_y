# calculate RB,PB,coverage rate
library(mice, warn.conflicts = FALSE)
library(tidyverse)

folder <- "C:/Users/20202596/Documents/R/mice_eval"
file_names <-list.files(folder)
output_path <- "C:/Users/20202596/Documents/R/data/with y"

truth_train <- read.csv("C:/Users/20202596/Documents/R/data/output/truth_train.csv")
truth_test <- read.csv("C:/Users/20202596/Documents/R/data/output/truth_test.csv")

truth_train_means <- data.frame(map(truth_train, ~mean(.))[-6])
truth_test_means <- data.frame(map(truth_test, ~mean(.))[-6])

results_list <- list()

simulate <- function(dataset, imp_method, truth_means, runs = 10) {
  mean_df <- data.frame(matrix(nrow=runs, ncol=dim(dataset)[2]-1))
  var_df  <- data.frame(matrix(nrow=runs, ncol=dim(dataset)[2]-1))
  rb_df   <- data.frame(matrix(nrow=runs, ncol=dim(dataset)[2]-1))
  pb_df   <- data.frame(matrix(nrow=runs, ncol=dim(dataset)[2]-1))
  names(mean_df) <- c('exang', 'trestbps', 'chol', 'sex', 'thalach')
  names(var_df) <- c('exang', 'trestbps', 'chol', 'sex', 'thalach')
  names(rb_df) <- c('exang', 'trestbps', 'chol', 'sex', 'thalach')
  names(pb_df) <- c('exang', 'trestbps', 'chol', 'sex', 'thalach')
  
  for (run in seq(runs))
  {
    if (run %% 10 == 0)
    {
      print(paste("Simulation run", run))
    }
    
    imp <- mice(dataset,maxit=20,m=1,remove.collinear = FALSE, method =imp_method, printFlag = FALSE)#"pmm" )
    imputed <- complete(imp)
    means <- imp$chainMean[,20,][-6] # -6 indexing to ignore hd column
    vars  <- imp$chainVar[,20,][-6]
    mean_df[run,] = means
    var_df[run,]  = vars
    rb_df[run,]   = means - truth_means
    pb_df[run,]   = 100*abs((means-truth_means)/truth_means)
  }
  
  return(list(mean_df, var_df, rb_df, pb_df))
}

for (name in file_names){
  import_path <- paste(folder,name,sep = "/")
  missing_data <- read.csv(import_path)
  if (dim(missing_data)[1] > 100)
  {
    # This is a set of training data, so use the training data means as truth
    truth_means <- truth_train_means
  }
  else
  {
    truth_means <- truth_test_means
  }
  
  mean_mat_list <- list()
  print(paste("Simulating results for", name))
  results_w_full_pmm <- simulate(missing_data, "pmm", truth_means, 100)
  pb_means <- data.frame(map(results_w_full_pmm[[4]], ~mean(.)))
  file_results <- list(data.frame(map(results_w_full_pmm[[1]], ~mean(.))),
                       data.frame(map(results_w_full_pmm[[2]], ~mean(.))),
                       data.frame(map(results_w_full_pmm[[3]], ~mean(.))),
                       data.frame(map(results_w_full_pmm[[4]], ~mean(.))))
  results_list <- append(results_list, list(file_results))
}

names(results_list) <- file_names
