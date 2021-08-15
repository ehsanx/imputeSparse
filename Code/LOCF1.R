################# IPAW with LOCF #################
# IPAW function
LOCF.ipaw.func1 <- function(iter = 1,                       # simulation index
                            K = 60,                         # Max follow-up time (in months)
                            n = 1000,                       # subject per arm
                            m = 3,                          # gap in measurement (in months)
                            sim.new.data = TRUE,            # whether to generate data or read old data
                            save.new.sim.data = TRUE,       # whether to save generated data
                            mainPath = "~/GitHub/LOCF_MI/", # main directory to save
                            scenario.code = 1,              # simulation scenario indicator
                            sigma = 2,                      # SD of L1  
                            theta2 = 0,                     # true effect 
                            B.unmeasured = FALSE,           # if U is unmeasured
                            RD = FALSE,                     # If effect estimate is RD or OR 
                            trunc = 0.01,                   # truncation level
                            beta1_0=0,                      # Parameters in data generation
                            beta2_0=-5,                     # Parameters in data generation 
                            beta1_1=6,                      # Parameters in data generation
                            beta2_1=3,                      # Parameters in data generation
                            theta0=-11,                     # Parameters in data generation
                            theta1=8                        # Parameters in data generation
                            ){ 
  if (sim.new.data == TRUE){
    simdat <- produce.data(iter = iter, K = K, n = n, m = m, save.data = save.new.sim.data, 
                           sigma = sigma, theta2 = theta2, mainPath = mainPath, 
                           scenario.code = scenario.code, 
                           beta1_0=beta1_0, beta2_0=beta2_0, beta1_1=beta1_1, beta2_1=beta2_1, 
                           theta0=theta0, theta1=theta1)
    cat("Data", iter, "created for LOCF\n")
  } else {
    scenario.code1 <- as.numeric(as.character(substr(scenario.code, 1, 2)))
    simdat <- readRDS(file = paste0(mainPath, paste0("simulations/scenario", scenario.code1,"/m",m,"/Data", 
                                                     iter,"sigma",sigma, ".RDS")))
    cat("Data", iter, "read from local repo for LOCF\n")
  }
  dat0 <- simdat
  dat <- dat0 
  dat$L1 <- na.locf(dat0$L1) 
  dat <- dat %>% 
    group_by(id) %>%
    mutate(cavgL1 = cumsum(L1)/(t0+1)) 
  
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
  simdatAC1$ipaw.trunc[simdatAC1$ipaw <= quantile.NA.rm(simdatAC1$ipaw, trunc)] <- quantile.NA.rm(simdatAC1$ipaw, trunc)
  simdatAC1$ipaw.trunc[simdatAC1$ipaw >  quantile.NA.rm(simdatAC1$ipaw, 1-trunc)] <- quantile.NA.rm(simdatAC1$ipaw, 1-trunc)

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
  simdatAC0$ipaw.trunc[simdatAC0$ipaw <= quantile.NA.rm(simdatAC0$ipaw, trunc)] <- quantile.NA.rm(simdatAC0$ipaw, trunc)
  simdatAC0$ipaw.trunc[simdatAC0$ipaw >  quantile.NA.rm(simdatAC0$ipaw, 1-trunc)] <- quantile.NA.rm(simdatAC0$ipaw, 1-trunc)
  
  dat.ipaw <- rbind(simdatAC1, simdatAC0)
  dat.ipaw <- dat.ipaw[dat.ipaw$Z == dat.ipaw$Alag1, ]
  
  # Weight summary
  ipaw.mean = with(dat.ipaw, mean(ipaw, na.rm = FALSE))
  ipaw.sd = with(dat.ipaw, sd(ipaw, na.rm = FALSE))
  ipaw.min = with(dat.ipaw, min(ipaw, na.rm = FALSE))
  ipaw.max = with(dat.ipaw, max(ipaw, na.rm = FALSE))
  ipaw.q1 = with(dat.ipaw, quantile.NA.rm(ipaw, 0.25))
  ipaw.q3 = with(dat.ipaw, quantile.NA.rm(ipaw, 0.75))
  
  ipaw.mean.trunc = with(dat.ipaw, mean(ipaw.trunc, na.rm = FALSE))
  ipaw.sd.trunc = with(dat.ipaw, sd(ipaw.trunc, na.rm = FALSE))
  ipaw.min.trunc = with(dat.ipaw, min(ipaw.trunc, na.rm = FALSE))
  ipaw.max.trunc = with(dat.ipaw, max(ipaw.trunc, na.rm = FALSE))
  ipaw.q1.trunc = with(dat.ipaw, quantile.NA.rm(ipaw.trunc, 0.25))
  ipaw.q3.trunc = with(dat.ipaw, quantile.NA.rm(ipaw.trunc, 0.75))
  
  ## IPAW of stabilized weight - sandwich SE
  if (RD == FALSE){
    if (B.unmeasured == FALSE) {
      fit.ipaw <- glm2(Y ~ A + B + t0, data = dat.ipaw, weights = ipaw, 
                       family = binomial("logit"))
    } else {
      fit.ipaw <- glm2(Y ~ A + t0, data = dat.ipaw, weights = ipaw, 
                       family = binomial("logit"))
    }
  }
  if (RD == TRUE){
    if (B.unmeasured == FALSE) {
      fit.ipaw <- lm(Y ~ A + B + t0, data = dat.ipaw, weights = ipaw)
    } else {
      fit.ipaw <- lm(Y ~ A + t0, data = dat.ipaw, weights = ipaw)
    }
  }
  cov.m1 <- vcovCL(fit.ipaw, type = "HC1", cluster = dat.ipaw$id)
  co.test <- coeftest(fit.ipaw, vcov = cov.m1)
  co.CI <- coefci(fit.ipaw, vcov = cov.m1, level = 0.95)
  
  # Summary
  ipaw.m.logHR <- coef(co.test)[["A"]]
  ipaw.m.HR <- exp(ipaw.m.logHR)
  ipaw.SE <- sqrt(diag(cov.m1))[["A"]]
  ipaw.LB <- co.CI["A","2.5 %"]
  ipaw.UB <- co.CI["A","97.5 %"]
  if (RD == FALSE) ipaw.p <- co.test["A","Pr(>|z|)"]
  if (RD == TRUE) ipaw.p <- co.test["A","Pr(>|t|)"]
  
  ## IPAW of truncated stabilized weight - sandwich SE
  if (RD == FALSE){
    if (B.unmeasured == FALSE) {
      fit.ipaw.trunc <- glm2(Y ~ A + B + t0, data = dat.ipaw, weights = ipaw.trunc, 
                             family = binomial("logit"))
    } else {
      fit.ipaw.trunc <- glm2(Y ~ A + t0, data = dat.ipaw, weights = ipaw.trunc, 
                             family = binomial("logit"))
    }
  }
  if (RD == TRUE){
    if (B.unmeasured == FALSE) {
      fit.ipaw.trunc <- lm(Y ~ A + B + t0, data = dat.ipaw, weights = ipaw.trunc)
    } else {
      fit.ipaw.trunc <- lm(Y ~ A + t0, data = dat.ipaw, weights = ipaw.trunc)
    }
  }
  
  cov.m1.trunc <- vcovCL(fit.ipaw.trunc, type = "HC1", cluster = dat.ipaw$id)
  co.test.trunc <- coeftest(fit.ipaw.trunc, vcov = cov.m1.trunc)
  co.CI.trunc <- coefci(fit.ipaw.trunc, vcov = cov.m1.trunc, level = 0.95)
  
  # Summary for truncated
  ipaw.m.logHR.trunc <- coef(co.test.trunc)[["A"]]
  ipaw.m.HR.trunc <- exp(ipaw.m.logHR.trunc)
  ipaw.SE.trunc <- sqrt(diag(cov.m1.trunc))[["A"]]
  ipaw.LB.trunc <- co.CI.trunc["A","2.5 %"]
  ipaw.UB.trunc <- co.CI.trunc["A","97.5 %"]
  if (RD == FALSE) ipaw.p.trunc <- co.test.trunc["A","Pr(>|z|)"]
  if (RD == TRUE) ipaw.p.trunc <- co.test.trunc["A","Pr(>|t|)"]
  
  cn1 <- numprobA1$converged == FALSE
  cd1 <- denomprobA1$converged == FALSE
  cn0 <- numprobA0$converged == FALSE
  cd0<- denomprobA0$converged == FALSE
  if (RD == FALSE) cf0 <- fit.ipaw$converged == FALSE
  if (RD == TRUE) cf0 <- TRUE
  if (RD == FALSE) cft <- fit.ipaw.trunc$converged == FALSE
  if (RD == TRUE) cft <- TRUE
  
  res.est <- c(ipaw.m.logHR, ipaw.SE, ipaw.LB, ipaw.UB,ipaw.p,
               ipaw.mean, ipaw.sd, ipaw.min, ipaw.max, ipaw.q1, ipaw.q3,
               ipaw.m.logHR.trunc, ipaw.SE.trunc, ipaw.LB.trunc, ipaw.UB.trunc,ipaw.p.trunc,
               ipaw.mean.trunc, ipaw.sd.trunc, ipaw.min.trunc, ipaw.max.trunc, ipaw.q1.trunc, ipaw.q3.trunc,
               cn1, cd1, cn0, cd0, cf0, cft)
  return(res.est)
}

# LOCF.ipaw.func1(iter = 2, sim.new.data = TRUE, save.new.sim.data = TRUE, mainPath = mainPath)
# LOCF.ipaw.func1(iter = 2, sim.new.data = FALSE)
