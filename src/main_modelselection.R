library(cstmr)
library(stringr)
library(lubridate)
library(foreach)
library(doParallel)

path_out = "fiftyshadesofgrey/io/out/"
path_in = "fiftyshadesofgrey/io/in/"

## Set the working directory. Change this to the location of the example on the computer. Note that "/" is always used in R, also in Windows
setwd("fiftyshadesofgrey/src/")

# Source the scripts with functions in the "src" folder. Just a neat way of arranging helping functions in R
sapply(dir("allmodels", full.names=TRUE), source)
sapply(dir("utils", full.names=TRUE), source)

# RC models and forward selection scheme declaration
RC_models = list("Ti"=Ti, 
                 # 2nd ite
                 "TiTe"=TiTe, "TiTh"=TiTh, "TiTm"=TiTm, "TiTs" =TiTs,
                 # 2rd ite
                 "TiTmTe"=TiTmTe, "TiTmTs"=TiTmTs, "TiTmTh"=TiTmTh, "TiTeTs"=TiTeTs,
                 "TiTeTh"=TiTeTh, "TiThTs"=TiThTs, "TiTeAeRia"=TiTeAeRia,
                 # 4th ite
                 "TiTmTeAeRia"=TiTmTeAeRia, "TiTmTeTh"=TiTmTeTh, "TiTmTeTs"=TiTmTeTs, "TiTmThTs"=TiTmThTs,
                 "TiTeTsAeRia"=TiTeTsAeRia, "TiTeThTs"=TiTeThTs, "TiTeThAeRia"=TiTeThAeRia,
                 # 5th ite
                 "TiTeThTsAeRia"=TiTeThTsAeRia, "TiTmTeTsAeRia"=TiTmTeTsAeRia, "TiTmTeThAeRia"=TiTmTeThAeRia,
                 "TiTmTeThTs"=TiTmTeThTs,
                 # 6th ite
                 "TiTmTeThTsAeRia"=TiTmTeThTsAeRia)

RC_iterations = list("0" = list("Ti"), 
                     # 1st -> 2nd ite
                     "Ti" = list("TiTe", "TiTh", "TiTm", "TiTs"),
                     # 2nd -> 3rd ite
                     "TiTe" = list("TiTeAeRia", "TiTmTe", "TiTeTs", "TiTeTh"),
                     "TiTh" = list("TiTeTh", "TiThTs", "TiTmTh"),
                     "TiTm" = list("TiTmTe", "TiTmTs", "TiTmTh"),
                     "TiTs" = list("TiTeTs", "TiTmTs", "TiThTs"),
                     # 3rd -> 4th ite
                     "TiTmTe" = list("TiTmTeAeRia", "TiTmTeTh", "TiTmTeTs"),
                     "TiTmTs" = list("TiTmTeTs", "TiTmThTs"),
                     "TiTmTh" = list("TiTmTeTh", "TiTmThTs"),
                     "TiTeTs" = list("TiTeTsAeRia", "TiTeThTs", "TiTmTeTs"),
                     "TiTeTh" = list("TiTeThAeRia", "TiTeThTs", "TiTmTeTh"),
                     "TiTeAeRia" = list("TiTmTeAeRia", "TiTeTsAeRia", "TiTeThAeRia"),
                     "TiThTs" = list("TiTeThTs", "TiTmThTs"),
                     # 4th -> 5th ite
                     "TiTmTeAeRia" = list("TiTmTeTsAeRia", "TiTmTeThAeRia"),
                     "TiTmTeTh" = list("TiTmTeThAeRia", "TiTmTeThTs"),
                     "TiTmTeTs" = list("TiTmTeTsAeRia", "TiTmTeThTs"),
                     "TiTmThTs" = list("TiTmTeThTs"),
                     "TiTeTsAeRia" = list("TiTeThTsAeRia", "TiTmTeTsAeRia"),
                     "TiTeThTs" = list("TiTeThTsAeRia", "TiTmTeThTs"),
                     "TiTeThAeRia" = list("TiTeThTsAeRia", "TiTmTeThAeRia"),
                     # 5th -> 6th ite
                     "TiTeThTsAeRia" = list("TiTmTeThTsAeRia"),
                     "TiTmTeTsAeRia" = list("TiTmTeThTsAeRia"),
                     "TiTmTeThAeRia" = list("TiTmTeThTsAeRia"),
                     "TiTmTeThTs" = list("TiTmTeThTsAeRia")
                     )

# Initial parameter declaration
inl <- c()
inl[["Tm_ini"]] <- c(19, 19, 19, 20, 18)
inl[["Tm_lb"]] <- c(12, 12, 10, 5, 5)
inl[["Tm_ub"]] <- c(30, 30, 30, 35, 40)
inl[["Te_lb"]] <- c(-10, -10, 0, 0, 0)
inl[["Te_ub"]] <- c(40, 40, 40, 35, 35)
inl[["Ci_ini"]] <- c(5, 1, 3, 10, 15, 50)
inl[["Ci_ub"]] <- c(20, 30, 40, 60, 100, 300)
inl[["Ria_ini"]] <- c(5, 5, 5, 10, 1)
inl[["Rie_ini"]] <- c(5, 10, 10, 10, 1)
inl[["Rea_ini"]] <- c(5, 5, 10, 10, 5)
inl[["Rih_ini"]] <- c(5, 5, 10, 10, 1)
inl[["Rim_ini"]] <- c(1, 1, 3, 5, 1)
inl[["Ris_ini"]] <- c(1, 1, 3, 5, 1)
inl[["Rh_ini"]] <- c(0.5, 0.5,0.5, 1, 2)
inl[["p_ini"]] <- c(1, 1, 1, -1, 5)
inl[["p_lb"]] <- c(-40,-40,-40, -50, -20)
inl[["p_ub"]] <- c(20, 20, 20, 10, 40)
inl[["e_ini"]] <- c(-1, -1, 1, 1, 0.5)
inl[["e_lb"]] <- c(-50, -40, -40, -40, -20)
inl[["e_ub"]] <- c(20, 20, 20, 30, 40)


# Get all input files names
setwd(path_in)
list_file_names <- list.files(pattern = "blg_*")
# Save buildings with no converging model fit
list_nofit <- c()


### Parallel Looping over building input files
n.cores <- parallel::detectCores() - 1
#create the cluster
cl <- parallel::makeCluster(
  n.cores, 
  type = "FORK", #FORK, PSOCK
  outfile='log.txt'
  )
registerDoParallel(cl)

# PARALLEL LOOP
all_list_nofit <- foreach (file = list_file_names, 
                           .combine=cbind) %dopar% {
  # library(stringr)
  # library(lubridate)

  ### ---PREPROCESSING------------------------------------------------------------------------------------------------------
  ## Read the csv file
  df <- read.csv(paste(path_in,file,sep=""),sep=";",header=TRUE)
  # Get blg uuid from file name
  split_file_name <- str_split_fixed(file, "_", 2)
  uuid_filtering <- sapply(split_file_name, tail, 1)[[2]]
  uuid_filtering <- str_split_fixed(uuid_filtering, ".csv", 2)[1]

  # Convert to datetime
  df$timedate <- ymd_hms(df$t)
  ## df$t is now hours since start of the experiment.
  df$t <- seq(0, length(df$t)-1, by=1) / 4  # dt = 15 minutes
  
  # if (all(df$Th==0)){ # Testing boiler temperature data
  #   print('Elec heater detected')
  #   return(uuid_filtering)
  # }

  # Define input X
  X <- df[c("t", "Ps", "Ta", "Thm", "yTi", "timedate")]
  # Removing NA values
  X = X[complete.cases(X), ]
  
  # Automated check of value initialization from the data set
  inl[["Th_ini"]] <- c(X$Thm[1], X$Thm[1], X$Thm[1], 40, 45)
  inl[["Th_ub"]] <- c(90, 90, 90, 90, 90)
  
  if (X$yTi[1] - X$Ta[1] <= 0) { #testing initial conditions 1.1
    inl[["Te_ini"]] <- c(1, 1, 1, 10, 15)
  } else if (X$yTi[1] - X$Ta[1] >= 35) { #testing initial conditions 1.2
    inl[["Te_ini"]] <- c(15, 15, 15, 10, 15)
  } else{
    inl[["Te_ini"]] <- c(X$yTi[1] - X$Ta[1], X$yTi[1] - X$Ta[1], X$yTi[1] - X$Ta[1], 10, 15)
  }
  if ( all(X$Thm < 1) & all(X$Thm > 0) ) { #testing initial conditions 2.1
    print('This boiler is binary')
    print(uuid_filtering)
    inl[["Th_ini"]] <- c(0.1, 0.1, 0.1, 0.1, 0.1)
    inl[["Th_ub"]] <- c(1, 1, 1, 1, 1)
  } else if (X$Thm[1] < 10) { #testing initial conditions 2.2
    inl[["Ti_lb"]] <- c(5, 5, 5, 15, 15)
    inl[["Ti_ub"]] <- c(30, 30, 25, 35, 35)
  } else if (X$Thm[1] > 35) { #testing initial conditions 2.3
    inl[["Ti_lb"]] <- c(10, 10, 10, 15, 15)
    inl[["Ti_ub"]] <- c(40, 40, 40, 45, 45)
  } else {
    inl[["Ti_lb"]] <- c(10, 10, 10, 15, 15)
    inl[["Ti_ub"]] <- c(30, 30, 25, 35, 35)
  }
  
  
  ### ---- GREYBOX MODEL FITTING ------------------------------------------------------------
  while_condition <- TRUE
  ite <- "0"

  ### Model selection
  while (while_condition) {
    no_fit_condition <- TRUE
    
    # Structure to save max loglikelihood
    loglik_ite <- c()
    
    ### Looping over models in each iteration
    for (m in RC_iterations[[ite]]) { # m - model loop
      
      # Select model to fit
      RC_model_considered <- RC_models[[m]]
      
      # Looping over multiple initial values
      for (i in (1:5)) { # i - intial parameter values loop
        # Parameter dictionary (list)
        input_param <- input_ini(X, inl, i)
        ## ----executeTiTe,results="hide"------------------------------------------
        try(fited_model <- RC_model_considered(X, input_param))
        # Results extraction - appending loglik_ite list with the current model's log-likelihood
        try(loglik_ite <- c(loglik_ite, fited_model$loglik))
        # Saving best current model - if its loglikelihood is the maximum observed yet
        try(
          if (tail(loglik_ite, n=1) == max(loglik_ite)) { # max likelyhood test
            fited_model$model_name <- m
            best_fit <- fited_model
            best_model_name <- m
            no_fit_condition <- FALSE
        })
      } # End of [ite][model][initial values] fitting loop
      
    } # End of [ite][model] fitting loop
    
    
    # Model did not converge
    if(no_fit_condition){
      while_condition <- FALSE
      return(uuid_filtering)
    
    # Model did converge - 1st ite save
    } else if (ite == "0") {
      # Update the model list to loop over for next iteration
      ite <- best_model_name
      best_fit$best_path <- c(best_model_name)
      # Simply save the best model (max likelyhood selection)
      save(best_fit, file=paste(path_out,toString(uuid_filtering),'_fit.rda', sep="")) 
      
    # Model did converge - best ite save
    } else {
      # Saving the current best model from this iteration
      best_fit_new <- best_fit
      # Comparing this model with the previous one
      load(paste(path_out,toString(uuid_filtering),'_fit.rda', sep=""))
      ## Calculate lambda
      chisqStat <- -2 * (max(best_fit$loglik) - max(best_fit_new$loglik))
      ## If this gives a p-value smaller than confidence limit, i.e. 5\%, then the
      ## new model is significantly better than the previous one
      prmDiff <- best_fit_new$model$NPARAM - best_fit$model$NPARAM
      ## Calculating the p-value of the test
      p_val <- 1 - pchisq(chisqStat, prmDiff)
      
      if (p_val < 0.05) { # If the p-value shows significance - we save the new model
          # New model is significantly better than the last
          # Saving best model
          best_newpath <- c(best_fit$best_path, best_model_name)
          best_fit <- best_fit_new  # changing the naming convention of the best_fit model
          best_fit$best_path <- best_newpath
          save(best_fit, file=paste(path_out,toString(uuid_filtering),'_fit.rda', sep=""))
          # Update the model list to loop over for next iteration
          ite <- best_model_name
        } else {
          # New model shows no significant improvement of the previous model
          # We can terminate the model iteration
          while_condition <- FALSE
          return(NA)
        }
    }
    
    
  } # End of while loop
  
  list_nofit
} # End of [uuid] loop

stopCluster(cl)


save(all_list_nofit, file=paste(path_out,'all_list_nofit.csv', sep=""))