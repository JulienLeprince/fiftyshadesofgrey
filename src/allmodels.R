Ti <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps + 1/(Ci)*Rh*Thm )*dt + exp(p11)*dw1)
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTe <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) )*dt + exp(p22)*dw2)
  model$addSystem(dTi ~ ( 1/(Ci*Rie)*(Te-Ti) + Aw/Ci*Ps + 1/(Ci)*Rh*Thm )*dt + exp(p11)*dw1)
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTeAeRia <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) + Ae/Ce*Ps )*dt + exp(p22)*dw2)
  model$addSystem(dTi ~ ( 1/(Ci*Rie)*(Te-Ti) + 1/(Ci)*Rh*Thm + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps )*dt + exp(p11)*dw1)
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter(  Ae = c(init=P[["Ae_ini"]]  ,lb=P[["Ae_lb"]]  ,ub=P[["Ae_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTh <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rih)*(Th-Ti) + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTh ~ ( 1/(Ch*Rih)*(Ti-Th) + 1/(Ch)*Rh*Thm)*dt + exp(p22)*dw2 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Th = c(init=P[["Th_ini"]]  ,lb=P[["Th_lb"]]  ,ub=P[["Th_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Ch = c(init=P[["Ch_ini"]]  ,lb=P[["Ch_lb"]]  ,ub=P[["Ch_ub"]])  )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Rih = c(init=P[["Rih_ini"]] ,lb=P[["Rih_lb"]] ,ub=P[["Rih_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTm <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rim)*(Tm-Ti) + 1/(Ci*Ria)*(Ta-Ti) + 1/(Ci)*Rh*Thm + Aw/Ci*Ps )*dt + exp(p11)*dw1 )
  model$addSystem(dTm ~ ( 1/(Cm*Rim)*(Ti-Tm) )*dt + exp(p22)*dw2 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Tm = c(init=P[["Tm_ini"]]  ,lb=P[["Tm_lb"]]  ,ub=P[["Tm_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Cm = c(init=P[["Cm_ini"]]  ,lb=P[["Cm_lb"]]  ,ub=P[["Cm_ub"]])  )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Rim = c(init=P[["Rim_ini"]] ,lb=P[["Rim_lb"]] ,ub=P[["Rim_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTs <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Ris)*(Ts-Ti) + 1/(Ci*Rh)*Thm + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTs ~ ( 1/(Cs*Ris)*(Ti-Ts))*dt + exp(p22)*dw2 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ts is the state, yTi is the measured output
  model$addObs(yTi ~ Ts)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Ts = c(init=P[["Ts_ini"]]  ,lb=P[["Ts_lb"]]  ,ub=P[["Ts_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Cs = c(init=P[["Cs_ini"]]  ,lb=P[["Cs_lb"]]  ,ub=P[["Cs_ub"]])  )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Ris = c(init=P[["Ris_ini"]] ,lb=P[["Ris_lb"]] ,ub=P[["Ris_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTeTh <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rih)*(Th-Ti) + 1/(Ci*Rie)*(Te-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) )*dt + exp(p22)*dw2 )
  model$addSystem(dTh ~ ( 1/(Ch*Rih)*(Ti-Th) + Rh/Ch*Thm)*dt + exp(p33)*dw3 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Th = c(init=P[["Th_ini"]]  ,lb=P[["Th_lb"]]  ,ub=P[["Th_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter(  Ch = c(init=P[["Ch_ini"]]  ,lb=P[["Ch_lb"]]  ,ub=P[["Ch_ub"]])  )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Rih = c(init=P[["Rih_ini"]] ,lb=P[["Rih_lb"]] ,ub=P[["Rih_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTeThAeRia <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rih)*(Th-Ti) + 1/(Ci*Rie)*(Te-Ti) + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) + Ae/Ce*Ps )*dt + exp(p22)*dw2 )
  model$addSystem(dTh ~ ( 1/(Ch*Rih)*(Ti-Th) + Rh/Ch*Thm)*dt + exp(p33)*dw3 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Th = c(init=P[["Th_ini"]]  ,lb=P[["Th_lb"]]  ,ub=P[["Th_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter(  Ch = c(init=P[["Ch_ini"]]  ,lb=P[["Ch_lb"]]  ,ub=P[["Ch_ub"]])  )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Rih = c(init=P[["Rih_ini"]] ,lb=P[["Rih_lb"]] ,ub=P[["Rih_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter(  Ae = c(init=P[["Ae_ini"]]  ,lb=P[["Ae_lb"]]  ,ub=P[["Ae_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTmTe <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rim)*(Tm-Ti) + 1/Ci*Rh*Thm + 1/(Ci*Rie)*(Te-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) )*dt + exp(p22)*dw2 )
  model$addSystem(dTm ~ ( 1/(Cm*Rim)*(Ti-Tm))*dt + exp(p33)*dw3 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Tm = c(init=P[["Tm_ini"]]  ,lb=P[["Tm_lb"]]  ,ub=P[["Tm_ub"]])  )
  ## Set the initial value for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter(  Cm = c(init=P[["Cm_ini"]]  ,lb=P[["Cm_lb"]]  ,ub=P[["Cm_ub"]])  )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Rim = c(init=P[["Rim_ini"]] ,lb=P[["Rim_lb"]] ,ub=P[["Rim_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTmTeAeRia <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rim)*(Tm-Ti) + 1/Ci*Rh*Thm + 1/(Ci*Rie)*(Te-Ti) + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) + Ae/Ce*Ps )*dt + exp(p22)*dw2 )
  model$addSystem(dTm ~ ( 1/(Cm*Rim)*(Ti-Tm))*dt + exp(p33)*dw3 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Tm = c(init=P[["Tm_ini"]]  ,lb=P[["Tm_lb"]]  ,ub=P[["Tm_ub"]])  )
  ## Set the initial value for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter(  Cm = c(init=P[["Cm_ini"]]  ,lb=P[["Cm_lb"]]  ,ub=P[["Cm_ub"]])  )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Rim = c(init=P[["Rim_ini"]] ,lb=P[["Rim_lb"]] ,ub=P[["Rim_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter(  Ae = c(init=P[["Ae_ini"]]  ,lb=P[["Ae_lb"]]  ,ub=P[["Ae_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiThTs <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Ris)*(Ts-Ti) + 1/(Ci*Rih)*(Th-Ti) + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTh ~ ( 1/(Ch*Rih)*(Ti-Th) + 1/(Ch)*Rh*Thm)*dt + exp(p33)*dw3 )
  model$addSystem(dTs ~ ( 1/(Cs*Ris)*(Ti-Ts))*dt + exp(p22)*dw2 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ts is the state, yTi is the measured output
  model$addObs(yTi ~ Ts)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Th = c(init=P[["Th_ini"]]  ,lb=P[["Th_lb"]]  ,ub=P[["Th_ub"]])  )
  model$setParameter(  Ts = c(init=P[["Ts_ini"]]  ,lb=P[["Ts_lb"]]  ,ub=P[["Ts_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Ch = c(init=P[["Ch_ini"]]  ,lb=P[["Ch_lb"]]  ,ub=P[["Ch_ub"]])  )
  model$setParameter(  Cs = c(init=P[["Cs_ini"]]  ,lb=P[["Cs_lb"]]  ,ub=P[["Cs_ub"]])  )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Ris = c(init=P[["Ris_ini"]] ,lb=P[["Ris_lb"]] ,ub=P[["Ris_ub"]]) )
  model$setParameter( Rih = c(init=P[["Rih_ini"]] ,lb=P[["Rih_lb"]] ,ub=P[["Rih_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}
                  
TiTmTh <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rim)*(Tm-Ti) + 1/(Ci*Rih)*(Th-Ti) + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTh ~ ( 1/(Ch*Rih)*(Ti-Th) + 1/(Ch)*Rh*Thm)*dt + exp(p33)*dw3 )
  model$addSystem(dTm ~ ( 1/(Cm*Rim)*(Ti-Tm))*dt + exp(p22)*dw2 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Th = c(init=P[["Th_ini"]]  ,lb=P[["Th_lb"]]  ,ub=P[["Th_ub"]])  )
  model$setParameter(  Tm = c(init=P[["Tm_ini"]]  ,lb=P[["Tm_lb"]]  ,ub=P[["Tm_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Ch = c(init=P[["Ch_ini"]]  ,lb=P[["Ch_lb"]]  ,ub=P[["Ch_ub"]])  )
  model$setParameter(  Cm = c(init=P[["Cm_ini"]]  ,lb=P[["Cm_lb"]]  ,ub=P[["Cm_ub"]])  )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Rim = c(init=P[["Rim_ini"]] ,lb=P[["Rim_lb"]] ,ub=P[["Rim_ub"]]) )
  model$setParameter( Rih = c(init=P[["Rih_ini"]] ,lb=P[["Rih_lb"]] ,ub=P[["Rih_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTmTs <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Ris)*(Ts-Ti) + 1/(Ci*Rim)*(Tm-Ti) + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps + (1/Ci)*Rh*Thm )*dt + exp(p11)*dw1 )
  model$addSystem(dTm ~ ( 1/(Cm*Rim)*(Ti-Tm) )*dt + exp(p33)*dw3 )
  model$addSystem(dTs ~ ( 1/(Cs*Ris)*(Ti-Ts))*dt + exp(p22)*dw2 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ts is the state, yTi is the measured output
  model$addObs(yTi ~ Ts)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Tm = c(init=P[["Tm_ini"]]  ,lb=P[["Tm_lb"]]  ,ub=P[["Tm_ub"]])  )
  model$setParameter(  Ts = c(init=P[["Ts_ini"]]  ,lb=P[["Ts_lb"]]  ,ub=P[["Ts_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Cm = c(init=P[["Cm_ini"]]  ,lb=P[["Cm_lb"]]  ,ub=P[["Cm_ub"]])  )
  model$setParameter(  Cs = c(init=P[["Cs_ini"]]  ,lb=P[["Cs_lb"]]  ,ub=P[["Cs_ub"]])  )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Ris = c(init=P[["Ris_ini"]] ,lb=P[["Ris_lb"]] ,ub=P[["Ris_ub"]]) )
  model$setParameter( Rim = c(init=P[["Rim_ini"]] ,lb=P[["Rim_lb"]] ,ub=P[["Rim_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTeTs <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Ris)*(Ts-Ti) + 1/Ci*Rh*Thm + 1/(Ci*Rie)*(Te-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) )*dt + exp(p22)*dw2 )
  model$addSystem(dTs ~ ( 1/(Cs*Ris)*(Ti-Ts))*dt + exp(p33)*dw3 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## Set the observation equation: Ts is the state, yTi is the measured output
  model$addObs(yTi ~ Ts)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Ts = c(init=P[["Ts_ini"]]  ,lb=P[["Ts_lb"]]  ,ub=P[["Ts_ub"]])  )
  ## Set the initial value for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter(  Cs = c(init=P[["Cs_ini"]]  ,lb=P[["Cs_lb"]]  ,ub=P[["Cs_ub"]])  )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Ris = c(init=P[["Ris_ini"]] ,lb=P[["Ris_lb"]] ,ub=P[["Ris_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTeTsAeRia <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Ris)*(Ts-Ti) + 1/Ci*Rh*Thm + 1/(Ci*Rie)*(Te-Ti) + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) + Ae/Ce*Ps )*dt + exp(p22)*dw2 )
  model$addSystem(dTs ~ ( 1/(Cs*Ris)*(Ti-Ts))*dt + exp(p33)*dw3 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## Set the observation equation: Ts is the state, yTi is the measured output
  model$addObs(yTi ~ Ts)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Ts = c(init=P[["Ts_ini"]]  ,lb=P[["Ts_lb"]]  ,ub=P[["Ts_ub"]])  )
  ## Set the initial value for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter(  Cs = c(init=P[["Cs_ini"]]  ,lb=P[["Cs_lb"]]  ,ub=P[["Cs_ub"]])  )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Ris = c(init=P[["Ris_ini"]] ,lb=P[["Ris_lb"]] ,ub=P[["Ris_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter(  Ae = c(init=P[["Ae_ini"]]  ,lb=P[["Ae_lb"]]  ,ub=P[["Ae_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}


TiTmTeTh <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rim)*(Tm-Ti) + 1/(Ci*Rih)*(Th-Ti) + 1/(Ci*Rie)*(Te-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) )*dt + exp(p22)*dw2 )
  model$addSystem(dTh ~ ( 1/(Ch*Rih)*(Ti-Th) + 1/Ch*Rh*Thm)*dt + exp(p33)*dw3 )
  model$addSystem(dTm ~ ( 1/(Cm*Rim)*(Ti-Tm))*dt + exp(p44)*dw4 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Th = c(init=P[["Th_ini"]]  ,lb=P[["Th_lb"]]  ,ub=P[["Th_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Tm = c(init=P[["Tm_ini"]]  ,lb=P[["Tm_lb"]]  ,ub=P[["Tm_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Cm = c(init=P[["Cm_ini"]]  ,lb=P[["Cm_lb"]]  ,ub=P[["Cm_ub"]])  )
  model$setParameter(  Ch = c(init=P[["Ch_ini"]]  ,lb=P[["Ch_lb"]]  ,ub=P[["Ch_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter( Rim = c(init=P[["Rim_ini"]] ,lb=P[["Rim_lb"]] ,ub=P[["Rim_ub"]]) )
  model$setParameter( Rih = c(init=P[["Rih_ini"]] ,lb=P[["Rih_lb"]] ,ub=P[["Rih_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( p44 = c(init=P[["p44_ini"]] ,lb=P[["p44_lb"]] ,ub=P[["p44_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTmTeThAeRia <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rim)*(Tm-Ti) + 1/(Ci*Rih)*(Th-Ti) + 1/(Ci*Rie)*(Te-Ti) + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) + Ae/Ce*Ps )*dt + exp(p22)*dw2 )
  model$addSystem(dTh ~ ( 1/(Ch*Rih)*(Ti-Th) + 1/Ch*Rh*Thm)*dt + exp(p33)*dw3 )
  model$addSystem(dTm ~ ( 1/(Cm*Rim)*(Ti-Tm))*dt + exp(p44)*dw4 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Th = c(init=P[["Th_ini"]]  ,lb=P[["Th_lb"]]  ,ub=P[["Th_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Tm = c(init=P[["Tm_ini"]]  ,lb=P[["Tm_lb"]]  ,ub=P[["Tm_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Cm = c(init=P[["Cm_ini"]]  ,lb=P[["Cm_lb"]]  ,ub=P[["Cm_ub"]])  )
  model$setParameter(  Ch = c(init=P[["Ch_ini"]]  ,lb=P[["Ch_lb"]]  ,ub=P[["Ch_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter( Rim = c(init=P[["Rim_ini"]] ,lb=P[["Rim_lb"]] ,ub=P[["Rim_ub"]]) )
  model$setParameter( Rih = c(init=P[["Rih_ini"]] ,lb=P[["Rih_lb"]] ,ub=P[["Rih_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter(  Ae = c(init=P[["Ae_ini"]]  ,lb=P[["Ae_lb"]]  ,ub=P[["Ae_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( p44 = c(init=P[["p44_ini"]] ,lb=P[["p44_lb"]] ,ub=P[["p44_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTeThTs <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Ris)*(Ts-Ti) + 1/(Ci*Rih)*(Th-Ti) + 1/(Ci*Rie)*(Te-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) )*dt + exp(p22)*dw2 )
  model$addSystem(dTh ~ ( 1/(Ch*Rih)*(Ti-Th) + 1/Ch*Rh*Thm)*dt + exp(p33)*dw3 )
  model$addSystem(dTs ~ ( 1/(Cs*Ris)*(Ti-Ts))*dt + exp(p44)*dw4 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ts is the state, yTi is the measured output
  model$addObs(yTi ~ Ts)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Th = c(init=P[["Th_ini"]]  ,lb=P[["Th_lb"]]  ,ub=P[["Th_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Ts = c(init=P[["Ts_ini"]]  ,lb=P[["Ts_lb"]]  ,ub=P[["Ts_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Ch = c(init=P[["Ch_ini"]]  ,lb=P[["Ch_lb"]]  ,ub=P[["Ch_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter(  Cs = c(init=P[["Cs_ini"]]  ,lb=P[["Cs_lb"]]  ,ub=P[["Cs_ub"]])  )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter( Ris = c(init=P[["Ris_ini"]] ,lb=P[["Ris_lb"]] ,ub=P[["Ris_ub"]]) )
  model$setParameter( Rih = c(init=P[["Rih_ini"]] ,lb=P[["Rih_lb"]] ,ub=P[["Rih_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( p44 = c(init=P[["p44_ini"]] ,lb=P[["p44_lb"]] ,ub=P[["p44_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTeThTsAeRia <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Ris)*(Ts-Ti) + 1/(Ci*Rih)*(Th-Ti) + 1/(Ci*Rie)*(Te-Ti) + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) + Ae/Ce*Ps)*dt + exp(p22)*dw2 )
  model$addSystem(dTh ~ ( 1/(Ch*Rih)*(Ti-Th) + 1/Ch*Rh*Thm)*dt + exp(p33)*dw3 )
  model$addSystem(dTs ~ ( 1/(Cs*Ris)*(Ti-Ts))*dt + exp(p44)*dw4 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ts is the state, yTi is the measured output
  model$addObs(yTi ~ Ts)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Th = c(init=P[["Th_ini"]]  ,lb=P[["Th_lb"]]  ,ub=P[["Th_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Ts = c(init=P[["Ts_ini"]]  ,lb=P[["Ts_lb"]]  ,ub=P[["Ts_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Ch = c(init=P[["Ch_ini"]]  ,lb=P[["Ch_lb"]]  ,ub=P[["Ch_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter(  Cs = c(init=P[["Cs_ini"]]  ,lb=P[["Cs_lb"]]  ,ub=P[["Cs_ub"]])  )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter( Ris = c(init=P[["Ris_ini"]] ,lb=P[["Ris_lb"]] ,ub=P[["Ris_ub"]]) )
  model$setParameter( Rih = c(init=P[["Rih_ini"]] ,lb=P[["Rih_lb"]] ,ub=P[["Rih_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter(  Ae = c(init=P[["Ae_ini"]]  ,lb=P[["Ae_lb"]]  ,ub=P[["Ae_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( p44 = c(init=P[["p44_ini"]] ,lb=P[["p44_lb"]] ,ub=P[["p44_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTmThTs <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rim)*(Tm-Ti) + 1/(Ci*Ris)*(Ts-Ti) + 1/(Ci*Rih)*(Th-Ti) + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTh ~ ( 1/(Ch*Rih)*(Ti-Th) + 1/Ch*Rh*Thm)*dt + exp(p22)*dw2 )
  model$addSystem(dTm ~ ( 1/(Cm*Rim)*(Ti-Tm))*dt + exp(p33)*dw3 )
  model$addSystem(dTs ~ ( 1/(Cs*Ris)*(Ti-Ts))*dt + exp(p44)*dw4 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ts is the state, yTi is the measured output
  model$addObs(yTi ~ Ts)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Th = c(init=P[["Th_ini"]]  ,lb=P[["Th_lb"]]  ,ub=P[["Th_ub"]])  )
  model$setParameter(  Tm = c(init=P[["Tm_ini"]]  ,lb=P[["Tm_lb"]]  ,ub=P[["Tm_ub"]])  )
  model$setParameter(  Ts = c(init=P[["Ts_ini"]]  ,lb=P[["Ts_lb"]]  ,ub=P[["Ts_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Cm = c(init=P[["Cm_ini"]]  ,lb=P[["Cm_lb"]]  ,ub=P[["Cm_ub"]])  )
  model$setParameter(  Ch = c(init=P[["Ch_ini"]]  ,lb=P[["Ch_lb"]]  ,ub=P[["Ch_ub"]])  )
  model$setParameter(  Cs = c(init=P[["Cs_ini"]]  ,lb=P[["Cs_lb"]]  ,ub=P[["Cs_ub"]])  )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Ris = c(init=P[["Ris_ini"]] ,lb=P[["Ris_lb"]] ,ub=P[["Ris_ub"]]) )
  model$setParameter( Rim = c(init=P[["Rim_ini"]] ,lb=P[["Rim_lb"]] ,ub=P[["Rim_ub"]]) )
  model$setParameter( Rih = c(init=P[["Rih_ini"]] ,lb=P[["Rih_lb"]] ,ub=P[["Rih_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( p44 = c(init=P[["p44_ini"]] ,lb=P[["p44_lb"]] ,ub=P[["p44_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTmTeTs <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rim)*(Tm-Ti) + 1/(Ci*Ris)*(Ts-Ti) + 1/(Ci*Rie)*(Te-Ti) + 1/Ci*Rh*Thm + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) )*dt + exp(p22)*dw2 )
  model$addSystem(dTm ~ ( 1/(Cm*Rim)*(Ti-Tm))*dt + exp(p33)*dw3 )
  model$addSystem(dTs ~ ( 1/(Cs*Ris)*(Ti-Ts))*dt + exp(p44)*dw4 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ts is the state, yTi is the measured output
  model$addObs(yTi ~ Ts)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Tm = c(init=P[["Tm_ini"]]  ,lb=P[["Tm_lb"]]  ,ub=P[["Tm_ub"]])  )
  model$setParameter(  Ts = c(init=P[["Ts_ini"]]  ,lb=P[["Ts_lb"]]  ,ub=P[["Ts_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Cm = c(init=P[["Cm_ini"]]  ,lb=P[["Cm_lb"]]  ,ub=P[["Cm_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter(  Cs = c(init=P[["Cs_ini"]]  ,lb=P[["Cs_lb"]]  ,ub=P[["Cs_ub"]])  )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter( Ris = c(init=P[["Ris_ini"]] ,lb=P[["Ris_lb"]] ,ub=P[["Ris_ub"]]) )
  model$setParameter( Rim = c(init=P[["Rim_ini"]] ,lb=P[["Rim_lb"]] ,ub=P[["Rim_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( p44 = c(init=P[["p44_ini"]] ,lb=P[["p44_lb"]] ,ub=P[["p44_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTmTeTsAeRia <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rim)*(Tm-Ti) + 1/(Ci*Ris)*(Ts-Ti) + 1/(Ci*Rie)*(Te-Ti) + 1/Ci*Rh*Thm + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) + Ae/Ce*Ps )*dt + exp(p22)*dw2 )
  model$addSystem(dTm ~ ( 1/(Cm*Rim)*(Ti-Tm))*dt + exp(p33)*dw3 )
  model$addSystem(dTs ~ ( 1/(Cs*Ris)*(Ti-Ts))*dt + exp(p44)*dw4 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ts is the state, yTi is the measured output
  model$addObs(yTi ~ Ts)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Tm = c(init=P[["Tm_ini"]]  ,lb=P[["Tm_lb"]]  ,ub=P[["Tm_ub"]])  )
  model$setParameter(  Ts = c(init=P[["Ts_ini"]]  ,lb=P[["Ts_lb"]]  ,ub=P[["Ts_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Cm = c(init=P[["Cm_ini"]]  ,lb=P[["Cm_lb"]]  ,ub=P[["Cm_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter(  Cs = c(init=P[["Cs_ini"]]  ,lb=P[["Cs_lb"]]  ,ub=P[["Cs_ub"]])  )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter( Ris = c(init=P[["Ris_ini"]] ,lb=P[["Ris_lb"]] ,ub=P[["Ris_ub"]]) )
  model$setParameter( Rim = c(init=P[["Rim_ini"]] ,lb=P[["Rim_lb"]] ,ub=P[["Rim_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter(  Ae = c(init=P[["Ae_ini"]]  ,lb=P[["Ae_lb"]]  ,ub=P[["Ae_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( p44 = c(init=P[["p44_ini"]] ,lb=P[["p44_lb"]] ,ub=P[["p44_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}



TiTmTeThTs <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rim)*(Tm-Ti) + 1/(Ci*Ris)*(Ts-Ti) + 1/(Ci*Rih)*(Th-Ti) + 1/(Ci*Rie)*(Te-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) )*dt + exp(p22)*dw2 )
  model$addSystem(dTh ~ ( 1/(Ch*Rih)*(Ti-Th) + 1/Ch*Rh*Thm)*dt + exp(p33)*dw3 )
  model$addSystem(dTm ~ ( 1/(Cm*Rim)*(Ti-Tm))*dt + exp(p44)*dw4 )
  model$addSystem(dTs ~ ( 1/(Cs*Ris)*(Ti-Ts))*dt + exp(p55)*dw5 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ts is the state, yTi is the measured output
  model$addObs(yTi ~ Ts)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Th = c(init=P[["Th_ini"]]  ,lb=P[["Th_lb"]]  ,ub=P[["Th_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Tm = c(init=P[["Tm_ini"]]  ,lb=P[["Tm_lb"]]  ,ub=P[["Tm_ub"]])  )
  model$setParameter(  Ts = c(init=P[["Ts_ini"]]  ,lb=P[["Ts_lb"]]  ,ub=P[["Ts_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Cm = c(init=P[["Cm_ini"]]  ,lb=P[["Cm_lb"]]  ,ub=P[["Cm_ub"]])  )
  model$setParameter(  Ch = c(init=P[["Ch_ini"]]  ,lb=P[["Ch_lb"]]  ,ub=P[["Ch_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter(  Cs = c(init=P[["Cs_ini"]]  ,lb=P[["Cs_lb"]]  ,ub=P[["Cs_ub"]])  )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter( Ris = c(init=P[["Ris_ini"]] ,lb=P[["Ris_lb"]] ,ub=P[["Ris_ub"]]) )
  model$setParameter( Rim = c(init=P[["Rim_ini"]] ,lb=P[["Rim_lb"]] ,ub=P[["Rim_ub"]]) )
  model$setParameter( Rih = c(init=P[["Rih_ini"]] ,lb=P[["Rih_lb"]] ,ub=P[["Rih_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( p44 = c(init=P[["p44_ini"]] ,lb=P[["p44_lb"]] ,ub=P[["p44_ub"]]) )
  model$setParameter( p55 = c(init=P[["p55_ini"]] ,lb=P[["p55_lb"]] ,ub=P[["p55_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}

TiTmTeThTsAeRia <- function(Dat, P){
  ## Generate a new object of class ctsm
  model <- ctsm$new()
  # Model tolerance for numerical ODE solution
  model$options$odeeps <- 1E-5
  model$options$iEKFeps <- 1E-5
  model$options$eps <- 1E-5
  model$options$svdEps <- 1E-5
  ## Add system equations and thereby also states
  model$addSystem(dTi ~ ( 1/(Ci*Rim)*(Tm-Ti) + 1/(Ci*Ris)*(Ts-Ti) + 1/(Ci*Rih)*(Th-Ti) + 1/(Ci*Rie)*(Te-Ti) + 1/(Ci*Ria)*(Ta-Ti) + Aw/Ci*Ps  )*dt + exp(p11)*dw1 )
  model$addSystem(dTe ~ ( 1/(Ce*Rie)*(Ti-Te) + 1/(Ce*Rea)*(Ta-Te) + Ae/Ce*Ps )*dt + exp(p22)*dw2 )
  model$addSystem(dTh ~ ( 1/(Ch*Rih)*(Ti-Th) + 1/Ch*Rh*Thm)*dt + exp(p33)*dw3 )
  model$addSystem(dTm ~ ( 1/(Cm*Rim)*(Ti-Tm))*dt + exp(p44)*dw4 )
  model$addSystem(dTs ~ ( 1/(Cs*Ris)*(Ti-Ts))*dt + exp(p55)*dw5 )
  ## Set the names of the inputs
  model$addInput(Ta,Ps,Thm)
  ## ----addObs--------------------------------------------------------------
  ## Set the observation equation: Ts is the state, yTi is the measured output
  model$addObs(yTi ~ Ts)
  ## Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  ## ----initialValues,tidy=FALSE--------------------------------------------
  ## Set the initial value (for the optimization) of the states at the start time point
  model$setParameter(  Ti = c(init=P[["Ti_ini"]]  ,lb=P[["Ti_lb"]]  ,ub=P[["Ti_ub"]])  )
  model$setParameter(  Th = c(init=P[["Th_ini"]]  ,lb=P[["Th_lb"]]  ,ub=P[["Th_ub"]])  )
  model$setParameter(  Te = c(init=P[["Te_ini"]]  ,lb=P[["Te_lb"]]  ,ub=P[["Te_ub"]])  )
  model$setParameter(  Tm = c(init=P[["Tm_ini"]]  ,lb=P[["Tm_lb"]]  ,ub=P[["Tm_ub"]])  )
  model$setParameter(  Ts = c(init=P[["Ts_ini"]]  ,lb=P[["Ts_lb"]]  ,ub=P[["Ts_ub"]])  )
  ## Set the initial value of the parameters for the optimization
  model$setParameter(  Ci = c(init=P[["Ci_ini"]]  ,lb=P[["Ci_lb"]]  ,ub=P[["Ci_ub"]])  )
  model$setParameter(  Cm = c(init=P[["Cm_ini"]]  ,lb=P[["Cm_lb"]]  ,ub=P[["Cm_ub"]])  )
  model$setParameter(  Ch = c(init=P[["Ch_ini"]]  ,lb=P[["Ch_lb"]]  ,ub=P[["Ch_ub"]])  )
  model$setParameter(  Ce = c(init=P[["Ce_ini"]]  ,lb=P[["Ce_lb"]]  ,ub=P[["Ce_ub"]])  )
  model$setParameter(  Cs = c(init=P[["Cs_ini"]]  ,lb=P[["Cs_lb"]]  ,ub=P[["Cs_ub"]])  )
  model$setParameter( Ria = c(init=P[["Ria_ini"]] ,lb=P[["Ria_lb"]] ,ub=P[["Ria_ub"]]) )
  model$setParameter( Rie = c(init=P[["Rie_ini"]] ,lb=P[["Rie_lb"]] ,ub=P[["Rie_ub"]]) )
  model$setParameter( Rea = c(init=P[["Rea_ini"]] ,lb=P[["Rea_lb"]] ,ub=P[["Rea_ub"]]) )
  model$setParameter( Ris = c(init=P[["Ris_ini"]] ,lb=P[["Ris_lb"]] ,ub=P[["Ris_ub"]]) )
  model$setParameter( Rim = c(init=P[["Rim_ini"]] ,lb=P[["Rim_lb"]] ,ub=P[["Rim_ub"]]) )
  model$setParameter( Rih = c(init=P[["Rih_ini"]] ,lb=P[["Rih_lb"]] ,ub=P[["Rih_ub"]]) )
  model$setParameter( Rh  = c(init=P[["Rh_ini"]]  ,lb=P[["Rh_lb"]]  ,ub=P[["Rh_ub"]])  )
  model$setParameter(  Aw = c(init=P[["Aw_ini"]]  ,lb=P[["Aw_lb"]]  ,ub=P[["Aw_ub"]])  )
  model$setParameter(  Ae = c(init=P[["Ae_ini"]]  ,lb=P[["Ae_lb"]]  ,ub=P[["Ae_ub"]])  )
  model$setParameter( p11 = c(init=P[["p11_ini"]] ,lb=P[["p11_lb"]] ,ub=P[["p11_ub"]]) )
  model$setParameter( p22 = c(init=P[["p22_ini"]] ,lb=P[["p22_lb"]] ,ub=P[["p22_ub"]]) )
  model$setParameter( p33 = c(init=P[["p33_ini"]] ,lb=P[["p33_lb"]] ,ub=P[["p33_ub"]]) )
  model$setParameter( p44 = c(init=P[["p44_ini"]] ,lb=P[["p44_lb"]] ,ub=P[["p44_ub"]]) )
  model$setParameter( p55 = c(init=P[["p55_ini"]] ,lb=P[["p55_lb"]] ,ub=P[["p55_ub"]]) )
  model$setParameter( e11 = c(init=P[["e11_ini"]] ,lb=P[["e11_lb"]] ,ub=P[["e11_ub"]]) )
  ## Run the parameter optimization
  fit <- model$estimate(Dat,threads=prm$threads)
  ## Replace the data to have all series available for analysis
  fit$data[[1]] <- Dat
  ## Return the fit
  return(fit)
}