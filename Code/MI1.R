################# IPAW with MI #################
# IPAW function
MI.ipaw.func1 <- function(m.mi = 10,                      # number of imputations
                          iter = 1,                       # simulation index
                          K = 60,                         # Max follow-up time (in months)
                          n = 1000,                       # subject per arm
                          m = 24,                         # gap in measurement (in months)
                          imp.meth = "2l", 
                                                          # Possible options: 
                                                          # "l" = pmm, 
                                                          # "2l" = group pmm, 
                                                          # "m" = mean, 
                                                          # "2m" = group mean, 
                                                          # "r" = regression, 
                                                          # "br" = bayesian regression, 
                                                          # "rf" = random forest
                          sim.new.data = TRUE,            # whether to generate data or read old data
                          save.new.sim.data = TRUE,       # whether to save generated data
                          mainPath = "~/GitHub/LOCF_MI/", # main directory to save
                          scenario.code = 1,              # simulation scenario indicator
                          sigma = 2,                      # SD of L1  
                          theta2 = 0,                     # true effect 
                          B.unmeasured = FALSE,           # if U is unmeasured
                          RD = FALSE,                     # If effect estimate is RD or OR 
                          trunc = 0.01,                   # truncation level
                          beta1_0=0,                      # Intercept in group 0 in generating L
                          beta2_0=-5,                     # Intercept in group 1 in generating L
                          beta1_1=6,                      # slope for U in group 0 in generating L
                          beta2_1=3,                      # slope for U in group 1 in generating L 
                          theta0=-11,                     # intercept in the outcome model
                          theta1=8                        # slope for U in the outcome model
                          ){ 
  impdata1 <- vector("list", m.mi) 
  ipw.mean <- ipw.sd <- ipw.min <- ipw.max <- ipw.q1 <- ipw.q3 <- vector("list", m.mi)
  ipw.mean.trunc <- ipw.sd.trunc <- ipw.min.trunc <- ipw.max.trunc <- ipw.q1.trunc <- ipw.q3.trunc <- vector("list", m.mi)
  
  fit.ipaw <- fit.ipaw.trunc <- vector("list", m.mi)
  cn1 <- cd1 <- cn0 <- cd0 <- cf0 <- cft <- vector("list", m.mi)

  if (sim.new.data == TRUE){
    simdat <- produce.data(iter = iter, K = K, n = n, m = m, save.data = save.new.sim.data, 
                           sigma = sigma, theta2 = theta2, mainPath = mainPath,
                           scenario.code = scenario.code,
                           beta1_0=beta1_0, beta2_0=beta2_0, beta1_1=beta1_1, beta2_1=beta2_1, 
                           theta0=theta0, theta1=theta1)
    cat("Data", iter, "created for MI\n")
  } else {
    scenario.code1 <- as.numeric(as.character(substr(scenario.code, 1, 2)))
    simdat <- readRDS(file = paste0(mainPath, paste0("simulations/scenario", scenario.code1,"/m",m,"/Data", 
                                                     iter,"sigma",sigma, ".RDS")))
    cat("Data", iter, "read from local repo for MI\n")
  }
  dat <- simdat
  
  if (B.unmeasured == TRUE) dat <- dat[,-c("B", "cavgL1", "sumuse", "Ause", "avguse", "wgt_temp", "Auselag")]
  if (B.unmeasured == FALSE) dat <- dat[,-c("cavgL1", "sumuse", "Ause", "avguse", "wgt_temp", "Auselag")]
  
  dat$Nalen <- nelsonaalen(dat, t0, Y)
  
  # Step 0: Set imputation model
  ini <- mice(data = dat, maxit = 0, print = FALSE)
  pred <- ini$pred
  
  # Prediction model
  pred[,"id"] <- pred["id",] <- 0
  if (imp.meth == "2l") pred[,"id"] <- -2
  if (imp.meth == "2r") pred[,"id"] <- -2
  if (imp.meth == "2m") pred[,"id"] <- -2
  pred[,"Z"] <- pred["Z",] <- 0
  pred[,"L2"] <- pred["L2",] <- 0
  pred[,"A"] <- pred["A",] <- 0

  # Set imputation method
  meth <- ini$meth
  if (imp.meth == "2l") {
    meth["L1"] <- "2l.pmm"
  } 
  if (imp.meth == "l") {
    meth["L1"] <- "pmm"
  }
  if (imp.meth == "rf") {
    meth["L1"] <- "rfcont"
  }
  if (imp.meth == "r") {
    meth["L1"] <- "norm.nob"
  }
  if (imp.meth == "2r") {
    meth["L1"] <- "2L.norm"
  }
  if (imp.meth == "br") {
    meth["L1"] <- "norm"
  }
  if (imp.meth == "2m") {
    meth["L1"] <- "2l.groupmean"
  }
  if (imp.meth == "m") {
    meth["L1"] <- "mean"
  }

  imputation <- mice(data = dat, 
                     seed = iter, 
                     predictorMatrix = pred,
                     method = meth, 
                     m = m.mi, 
                     maxit = 1, 
                     print = FALSE)
  impdata <- mice::complete(imputation, action="long")
  
  #Remove .id variable from the model as it was created in an intermediate step
  impdata$.id <- NULL
  
  for (i in 1:m.mi) {
    impdata1[[i]] <- impdata[impdata$.imp==i,] %>% 
      group_by(id) %>% 
      mutate(cavgL1 = cumsum(L1)/(t0+1)) %>% 
      ungroup()
  }
  
  for (M in 1:m.mi) {
    dat <- impdata1[[M]]
    dat$wgt_temp <- 0.0
    ## For records with Z=1
    simdatAC1 <- dat[dat$Z==1 & dat$Alag1==1,]
    
    # Calculate numerator alpha_hat coeffs for stabilized IP weights
    if (B.unmeasured == FALSE) {
      numprobA1 <- glm(A ~ B + t0, data = simdatAC1, family = binomial(link = "logit"))
    } else {
      numprobA1 <- glm(A ~ t0, data = simdatAC1, family = binomial(link = "logit"))
    }
    
    # Calculate denominator alpha_hat coeffs for stabilized IP weights
    if (B.unmeasured == FALSE) {
      denomprobA1 <- glm(A ~ B + cavgL1 + L2lag + t0, data = simdatAC1, 
                       family = binomial(link = "logit"))
    } else {
      denomprobA1 <- glm(A ~ cavgL1 + L2lag + t0, data = simdatAC1, 
                         family = binomial(link = "logit"))
    }
    
    #weight contributions that are non-zero (records with A=1)
    simdatAC1[simdatAC1$A==1,]$wgt_temp <-
      predict(numprobA1, simdatAC1[simdatAC1$A==1,], type = "response")/
      predict(denomprobA1, simdatAC1[simdatAC1$A==1,], type = "response")
    
    #weight contributions that are zero (records with A=0)
    simdatAC1[simdatAC1$A==0,]$wgt_temp <- 0
    
    # Calc. final IP weights for complete data
    ipweights1 <- numeric()
    for(i in 1:n){
      data_subset <- subset(simdatAC1, id==i)
      # Cum. product of temp. weights for each indiv. over time
      ipweights1 <- c(ipweights1, cumprod(data_subset$wgt_temp))
    }
    
    # Append weights to data frame
    simdatAC1$ipaw <- ipweights1
    simdatAC1$ipaw.trunc <- simdatAC1$ipaw
    simdatAC1$ipaw.trunc[simdatAC1$ipaw <= 
                           quantile.NA.rm(simdatAC1$ipaw, trunc)] <- quantile.NA.rm(simdatAC1$ipaw, trunc)
    simdatAC1$ipaw.trunc[simdatAC1$ipaw >  
                           quantile.NA.rm(simdatAC1$ipaw, 1-trunc)] <- quantile.NA.rm(simdatAC1$ipaw, 1-trunc)
    
    ## For records with Z=0
    simdatAC0 <- dat[dat$Z==0 & dat$Alag1==0,]
    
    # Calculate numerator alpha_hat coeffs for stabilized IP weights
    if (B.unmeasured == FALSE) {
      numprobA0 <- glm(A ~ B + t0, data=simdatAC0, family = binomial(link = "logit"))
    } else {
      numprobA0 <- glm(A ~ t0, data=simdatAC0, family = binomial(link = "logit"))
    }
    
    # Calculate denominator alpha_hat coeffs for stabilized IP weights
    if (B.unmeasured == FALSE) {
      denomprobA0 <- glm(A ~ B + cavgL1 + L2lag + t0, data=simdatAC0, 
                       family=binomial(link = "logit"))
    } else {
      denomprobA0 <- glm(A ~ cavgL1 + L2lag + t0, data=simdatAC0, 
                         family=binomial(link = "logit"))
    }
    
    #weight contributions that are non-zero (records with A=0)
    simdatAC0[simdatAC0$A==0,]$wgt_temp <-
      (1-predict(numprobA0, simdatAC0[simdatAC0$A==0,],type = "response"))/
      (1-predict(denomprobA0, simdatAC0[simdatAC0$A==0,], type = "response"))
    
    #weight contributions that are zero (records with A=1)
    simdatAC0[simdatAC0$A==1,]$wgt_temp <- 0
    
    # Calc. final IP weights for complete data
    ipweights0 <- numeric()
    ipwuweights0 <- numeric()
    for(i in (n+1):(2*n)){
      data_subset <- subset(simdatAC0, id==i)
      # Cum. product of temp. weights for each indiv. over time
      ipweights0 <- c(ipweights0, cumprod(data_subset$wgt_temp))
    }
    # Append weights to data frame
    simdatAC0$ipaw <- ipweights0
    simdatAC0$ipaw.trunc <- simdatAC0$ipaw
    simdatAC0$ipaw.trunc[simdatAC0$ipaw <= 
                           quantile.NA.rm(simdatAC0$ipaw, trunc)] <- quantile.NA.rm(simdatAC0$ipaw, trunc)
    simdatAC0$ipaw.trunc[simdatAC0$ipaw >  
                           quantile.NA.rm(simdatAC0$ipaw, 1-trunc)] <- quantile.NA.rm(simdatAC0$ipaw, 1-trunc)
    
    dat.ipaw <- rbind(simdatAC1, simdatAC0)
    dat.ipaw <- dat.ipaw[dat.ipaw$Z == dat.ipaw$Alag1, ]
    
    # Weight summary
    ipw.mean[[M]] = with(dat.ipaw, mean(ipaw, na.rm = FALSE))
    ipw.sd[[M]] = with(dat.ipaw, sd(ipaw, na.rm = FALSE))
    ipw.min[[M]] = with(dat.ipaw, min(ipaw, na.rm = FALSE))
    ipw.max[[M]] = with(dat.ipaw, max(ipaw, na.rm = FALSE))
    ipw.q1[[M]] = with(dat.ipaw, quantile.NA.rm(ipaw, 0.25))
    ipw.q3[[M]] = with(dat.ipaw, quantile.NA.rm(ipaw, 0.75))
    
    ipw.mean.trunc[[M]] = with(dat.ipaw, mean(ipaw.trunc, na.rm = FALSE))
    ipw.sd.trunc[[M]] = with(dat.ipaw, sd(ipaw.trunc, na.rm = FALSE))
    ipw.min.trunc[[M]] = with(dat.ipaw, min(ipaw.trunc, na.rm = FALSE))
    ipw.max.trunc[[M]] = with(dat.ipaw, max(ipaw.trunc, na.rm = FALSE))
    ipw.q1.trunc[[M]] = with(dat.ipaw, quantile.NA.rm(ipaw.trunc, 0.25))
    ipw.q3.trunc[[M]] = with(dat.ipaw, quantile.NA.rm(ipaw.trunc, 0.75))
    
    ## IPAW of stabilized weight
    if (RD == FALSE){
      if (B.unmeasured == FALSE) {
        fit.ipaw[[M]] <- coeftest(glm2(Y ~ A + B + t0, data = dat.ipaw, weights = ipaw, 
                                       family = binomial("logit")), vcov. = sandwich,
                                  cluster = "id", type = "HC1")
      } else {
        fit.ipaw[[M]] <- coeftest(glm2(Y ~ A + t0, data = dat.ipaw, weights = ipaw, 
                                       family = binomial("logit")), vcov. = sandwich,
                                  cluster = "id", type = "HC1")
      }
      
      if (B.unmeasured == FALSE) {
        fit.ipaw.trunc[[M]] <- coeftest(glm2(Y ~ A + B + t0, data = dat.ipaw, weights = ipaw.trunc, 
                                             family = binomial("logit")), vcov. = sandwich,
                                        cluster = "id", type = "HC1")
      } else {
        fit.ipaw.trunc[[M]] <- coeftest(glm2(Y ~ A + t0, data = dat.ipaw, weights = ipaw.trunc, 
                                             family = binomial("logit")), vcov. = sandwich,
                                        cluster = "id", type = "HC1")
      }
    }
    if (RD == TRUE){
      if (B.unmeasured == FALSE) {
        fit.ipaw[[M]] <- coeftest(lm(Y ~ A + B + t0, data = dat.ipaw, weights = ipaw), 
                                  vcov. = sandwich,
                                  cluster = "id", type = "HC1")
      } else {
        fit.ipaw[[M]] <- coeftest(lm(Y ~ A + t0, data = dat.ipaw, weights = ipaw), 
                                  vcov. = sandwich,
                                  cluster = "id", type = "HC1")
      }
      
      if (B.unmeasured == FALSE) {
        fit.ipaw.trunc[[M]] <- coeftest(lm(Y ~ A + B + t0, data = dat.ipaw, weights = ipaw.trunc), 
                                        vcov. = sandwich,
                                        cluster = "id", type = "HC1")
      } else {
        fit.ipaw.trunc[[M]] <- coeftest(lm(Y ~ A + t0, data = dat.ipaw, weights = ipaw.trunc), 
                                        vcov. = sandwich,
                                        cluster = "id", type = "HC1")
      }
    }
    
    cn1[[M]] <- numprobA1$converged == FALSE
    cd1[[M]] <- denomprobA1$converged == FALSE
    cn0[[M]] <- numprobA0$converged == FALSE
    cd0[[M]] <- denomprobA0$converged == FALSE
    cf0[[M]] <- fit.ipaw$converged == FALSE
    cft[[M]] <- fit.ipaw.trunc$converged == FALSE
    
  }
  
  # Pooling the Results 
  est <- summary(pool(fit.ipaw))
  ipaw.m.logHR <- est$estimate[est$term=="A"]
  ipaw.m.HR <- exp(ipaw.m.logHR)
  ipaw.SE <- est$std.error[est$term=="A"]
  ipaw.LB <- ipaw.m.logHR - 1.96*ipaw.SE
  ipaw.UB <- ipaw.m.logHR + 1.96*ipaw.SE
  ipaw.p <- est$p.value[est$term=="A"]
  
  est.trunc <- summary(pool(fit.ipaw.trunc))
  ipaw.m.logHR.trunc <- est.trunc$estimate[est.trunc$term=="A"]
  ipaw.m.HR.trunc <- exp(ipaw.m.logHR.trunc)
  ipaw.SE.trunc <- est.trunc$std.error[est.trunc$term=="A"]
  ipaw.LB.trunc <- ipaw.m.logHR.trunc - 1.96*ipaw.SE.trunc
  ipaw.UB.trunc <- ipaw.m.logHR.trunc + 1.96*ipaw.SE.trunc
  ipaw.p.trunc <- est.trunc$p.value[est.trunc$term=="A"]

  ipaw.mean <- mean(unlist(ipw.mean))
  ipaw.sd <- mean(unlist(ipw.sd))
  ipaw.min <- mean(unlist(ipw.min))
  ipaw.max <- mean(unlist(ipw.max))
  ipaw.q1 <- mean(unlist(ipw.q1))
  ipaw.q3 <- mean(unlist(ipw.q3))
  
  ipaw.mean.trunc <- mean(unlist(ipw.mean.trunc))
  ipaw.sd.trunc <- mean(unlist(ipw.sd.trunc))
  ipaw.min.trunc <- mean(unlist(ipw.min.trunc))
  ipaw.max.trunc <- mean(unlist(ipw.max.trunc))
  ipaw.q1.trunc <- mean(unlist(ipw.q1.trunc))
  ipaw.q3.trunc <- mean(unlist(ipw.q3.trunc))
  
  cn1x <- sum(unlist(cn1), na.rm = FALSE)
  cd1x <- sum(unlist(cd1), na.rm = FALSE)
  cn0x <- sum(unlist(cn0), na.rm = FALSE)
  cd0x <- sum(unlist(cd0), na.rm = FALSE)
  cf0x <- sum(unlist(cf0), na.rm = FALSE)
  cftx <- sum(unlist(cft), na.rm = FALSE)
  
  res.est <- c(ipaw.m.logHR, ipaw.SE, ipaw.LB, ipaw.UB,ipaw.p, 
               ipaw.mean, ipaw.sd, ipaw.min, ipaw.max, ipaw.q1, ipaw.q3,
               ipaw.m.logHR.trunc, ipaw.SE.trunc, ipaw.LB.trunc, ipaw.UB.trunc,ipaw.p.trunc, 
               ipaw.mean.trunc, ipaw.sd.trunc, ipaw.min.trunc, ipaw.max.trunc, ipaw.q1.trunc, ipaw.q3.trunc,
               cn1x, cd1x, cn0x, cd0x, cf0x, cftx)
  return(res.est)
}

# MI.ipaw.func1(m.mi = 3, iter = 2,  m = 12, sim.new.data = TRUE, save.new.sim.data = FALSE, mainPath = mainPath, scenario.code = 99,imp.meth = "r" )
# MI.ipaw.func1(m.mi = 5, iter = 2,  m = 3, sim.new.data = FALSE, imp.meth = "r")
