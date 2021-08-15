require(data.table) 
require(zoo) 
require(dplyr) 
require(glm2) 
require(sandwich)
require(lmtest)
require(mice)
require(miceadds)
require(CALIBERrfimpute)
require(lme4)

source(paste0(mainPath, "codes/datagen.R"))
source(paste0(mainPath, "codes/LOCF1.R"))
source(paste0(mainPath, "codes/MI1.R"))

quantile.NA.rm <- function(x, probs = trunc, ...){
  quantile(x, probs, na.rm = TRUE, ...)
}

est.both <- function(m.mi = 10, iter = 1, K = 60, n = 1000, m = 24, #imp.meth = "2l",
                     #sim.new.data = TRUE, save.new.sim.data = FALSE, 
                     sigma = 2,theta2=0,
                     mainPath = "~/GitHub/LOCF_MI/", scenario.code = 1, B.unmeasured = FALSE,
                     beta1_0=0, beta2_0=-5, beta1_1=6, beta2_1=3, theta0=-11, theta1=8, RD = FALSE){
  simdat <- produce.data(iter = iter, K = K, n = n, m = m, save.data = TRUE, 
                         sigma = sigma, theta2 = theta2, mainPath = mainPath, 
                         scenario.code = scenario.code, 
                         beta1_0=beta1_0, beta2_0=beta2_0, beta1_1=beta1_1, beta2_1=beta2_1, 
                         theta0=theta0, theta1=theta1)
  names.l1 <- c("ipaw.m.logHR", "ipaw.SE", "ipaw.LB", "ipaw.UB","ipaw.p",
                "ipaw.mean", "ipaw.sd", "ipaw.min", "ipaw.max", "ipaw.q1", "ipaw.q3",
                "ipaw.m.logHR.trunc", "ipaw.SE.trunc", "ipaw.LB.trunc", "ipaw.UB.trunc","ipaw.p.trunc",
                "ipaw.mean.trunc", "ipaw.sd.trunc", "ipaw.min.trunc", "ipaw.max.trunc", "ipaw.q1.trunc", 
                "ipaw.q3.trunc","conv.n1", "conv.d1", "conv.n0", "conv.d0", "conv.fit", "conv.fit.tr")
  l1 <- tryCatch(LOCF.ipaw.func1(iter = iter, K = K, n = n, m = m, 
                                 sim.new.data = FALSE, save.new.sim.data = FALSE, 
                                 mainPath = mainPath, scenario.code = scenario.code, 
                                 sigma = sigma,theta2=theta2, B.unmeasured = B.unmeasured,
                                 beta1_0=beta1_0, beta2_0=beta2_0, beta1_1=beta1_1, beta2_1=beta2_1, 
                                 theta0=theta0, theta1=theta1, RD = RD),
                 error = function(e) { rep(NA, length(names.l1)) } )
  m1 <- tryCatch(MI.ipaw.func1(m.mi = m.mi, iter = iter, K = K, n = n, m = m, imp.meth = "l",
                               sim.new.data = FALSE, save.new.sim.data = FALSE, 
                               # generating twice, but if not, the saving is not fast enough from locf 
                               mainPath = mainPath, scenario.code = scenario.code, 
                               sigma = sigma,theta2=theta2, B.unmeasured = B.unmeasured,
                               beta1_0=beta1_0, beta2_0=beta2_0, beta1_1=beta1_1, beta2_1=beta2_1, 
                               theta0=theta0, theta1=theta1, RD = RD),
                 error = function(e) { rep(NA, length(names.l1)) } )  
  
  m2 <- tryCatch(MI.ipaw.func1(m.mi = m.mi, iter = iter, K = K, n = n, m = m, imp.meth = "2l",
                               sim.new.data = FALSE, save.new.sim.data = FALSE, 
                               # generating twice, but if not, the saving is not fast enough from locf 
                               mainPath = mainPath, scenario.code = scenario.code, 
                               sigma = sigma,theta2=theta2, B.unmeasured = B.unmeasured,
                               beta1_0=beta1_0, beta2_0=beta2_0, beta1_1=beta1_1, beta2_1=beta2_1, 
                               theta0=theta0, theta1=theta1, RD = RD),
                 error = function(e) { rep(NA, length(names.l1)) } )  
  
  m3 <- tryCatch(MI.ipaw.func1(m.mi = m.mi, iter = iter, K = K, n = n, m = m, imp.meth = "rf",
                               sim.new.data = FALSE, save.new.sim.data = FALSE, 
                               # generating twice, but if not, the saving is not fast enough from locf 
                               mainPath = mainPath, scenario.code = scenario.code, 
                               sigma = sigma,theta2=theta2, B.unmeasured = B.unmeasured,
                               beta1_0=beta1_0, beta2_0=beta2_0, beta1_1=beta1_1, beta2_1=beta2_1, 
                               theta0=theta0, theta1=theta1, RD = RD),
                 error = function(e) { rep(NA, length(names.l1)) } )  
  
  res <- c(l1,m1,m2,m3)
  names.m1 <- paste0(names.l1,".MI")
  names.m2 <- paste0(names.l1,".2MI")
  names.m3 <- paste0(names.l1,".rMI")
  
  names(res) <- c(names.l1,names.m1,names.m2,names.m3)
  return(res)
}


est.both.r <- function(m.mi = 10, iter = 1, K = 60, n = 1000, m = 24, #imp.meth = "2l",
                     #sim.new.data = TRUE, save.new.sim.data = FALSE, 
                     sigma = 2,theta2=0,
                     mainPath = "~/GitHub/LOCF_MI/", scenario.code = 1, B.unmeasured = FALSE,
                     beta1_0=0, beta2_0=-5, beta1_1=6, beta2_1=3, theta0=-11, theta1=8, RD = FALSE){
  names.l1 <- c("ipaw.m.logHR", "ipaw.SE", "ipaw.LB", "ipaw.UB","ipaw.p",
                "ipaw.mean", "ipaw.sd", "ipaw.min", "ipaw.max", "ipaw.q1", "ipaw.q3",
                "ipaw.m.logHR.trunc", "ipaw.SE.trunc", "ipaw.LB.trunc", "ipaw.UB.trunc","ipaw.p.trunc",
                "ipaw.mean.trunc", "ipaw.sd.trunc", "ipaw.min.trunc", "ipaw.max.trunc", "ipaw.q1.trunc", 
                "ipaw.q3.trunc","conv.n1", "conv.d1", "conv.n0", "conv.d0", "conv.fit", "conv.fit.tr")
  m1 <- tryCatch(MI.ipaw.func1(m.mi = m.mi, iter = iter, K = K, n = n, m = m, imp.meth = "r",
                               sim.new.data = FALSE, save.new.sim.data = FALSE, 
                               # generating twice, but if not, the saving is not fast enough from locf 
                               mainPath = mainPath, scenario.code = scenario.code, 
                               sigma = sigma,theta2=theta2, B.unmeasured = B.unmeasured,
                               beta1_0=beta1_0, beta2_0=beta2_0, beta1_1=beta1_1, beta2_1=beta2_1, 
                               theta0=theta0, theta1=theta1, RD = RD),
                 error = function(e) { rep(NA, length(names.l1)) } )  
  
  m2 <- tryCatch(MI.ipaw.func1(m.mi = m.mi, iter = iter, K = K, n = n, m = m, imp.meth = "2r",
                               sim.new.data = FALSE, save.new.sim.data = FALSE, 
                               # generating twice, but if not, the saving is not fast enough from locf 
                               mainPath = mainPath, scenario.code = scenario.code, 
                               sigma = sigma,theta2=theta2, B.unmeasured = B.unmeasured,
                               beta1_0=beta1_0, beta2_0=beta2_0, beta1_1=beta1_1, beta2_1=beta2_1, 
                               theta0=theta0, theta1=theta1, RD = RD),
                 error = function(e) { rep(NA, length(names.l1)) } )  
  
  m3 <- tryCatch(MI.ipaw.func1(m.mi = m.mi, iter = iter, K = K, n = n, m = m, imp.meth = "br",
                               sim.new.data = FALSE, save.new.sim.data = FALSE, 
                               # generating twice, but if not, the saving is not fast enough from locf 
                               mainPath = mainPath, scenario.code = scenario.code, 
                               sigma = sigma,theta2=theta2, B.unmeasured = B.unmeasured,
                               beta1_0=beta1_0, beta2_0=beta2_0, beta1_1=beta1_1, beta2_1=beta2_1, 
                               theta0=theta0, theta1=theta1, RD = RD),
                 error = function(e) { rep(NA, length(names.l1)) } )
  m4 <- tryCatch(MI.ipaw.func1(m.mi = m.mi, iter = iter, K = K, n = n, m = m, imp.meth = "2m",
                               sim.new.data = FALSE, save.new.sim.data = FALSE, 
                               # generating twice, but if not, the saving is not fast enough from locf 
                               mainPath = mainPath, scenario.code = scenario.code, 
                               sigma = sigma,theta2=theta2, B.unmeasured = B.unmeasured,
                               beta1_0=beta1_0, beta2_0=beta2_0, beta1_1=beta1_1, beta2_1=beta2_1, 
                               theta0=theta0, theta1=theta1, RD = RD),
                 error = function(e) { rep(NA, length(names.l1)) } )
  m5 <- tryCatch(MI.ipaw.func1(m.mi = m.mi, iter = iter, K = K, n = n, m = m, imp.meth = "m",
                               sim.new.data = FALSE, save.new.sim.data = FALSE, 
                               # generating twice, but if not, the saving is not fast enough from locf 
                               mainPath = mainPath, scenario.code = scenario.code, 
                               sigma = sigma,theta2=theta2, B.unmeasured = B.unmeasured,
                               beta1_0=beta1_0, beta2_0=beta2_0, beta1_1=beta1_1, beta2_1=beta2_1, 
                               theta0=theta0, theta1=theta1, RD = RD),
                 error = function(e) { rep(NA, length(names.l1)) } )
  
  res <- c(m1,m2,m3,m4,m5)
  names.m1 <- paste0(names.l1,".lrMI")
  names.m2 <- paste0(names.l1,".2rMI")
  names.m3 <- paste0(names.l1,".brMI")
  names.m4 <- paste0(names.l1,".2mMI")
  names.m5 <- paste0(names.l1,".mMI")
  
  names(res) <- c(names.m1,names.m2,names.m3,names.m4,names.m5)
  return(res)
}