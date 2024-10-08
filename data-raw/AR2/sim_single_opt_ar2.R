#' Simulates n data sets and performs bootstrapping
#'
#' Performs n tests
#'
#' For convenience the intercept, beta_0 is set to zero
#'
#'@param beta Numeric.  Values for beta_1 (slope/trend parameter)
#'@param rho Numeric  Values for autoregressive parameter
#'@param sigma Numeric.  Values of standard deviation of noise
#'@param nT Numeric.  Values for the length of time series to simulate
#'@param nSims Numeric.  Number of time series to simulate
#'@param nBootSims Numeric.  Number of bootstrap data sets
#'
#'@examples
#'\dontrun{
#'sim_single_opt(beta = 0.5,rho =0.5, sigma = 0.25, nT = 10,
#' nSims = 200, nBootSims = 500)
#'}
#'
#'@importFrom foreach %dopar%
#'
#'@export


sim_single_opt_ar2 <- function(beta = 0,
                       rho1 = 0,
                       rho2 = 0,
                       sigma = 0.25,
                       nT = 10,
                       nSims = 200,
                       nBootSims = 500) {

  # nC <- parallel::detectCores()
  # cl <- parallel::makeCluster(nC)
  # doParallel::registerDoParallel(cl)

  # allocate memory to a bunch of vectors
  LR <- vector(mode="numeric",length=nBootSims-1) # likelihood ratio statistic
  pVal_boot <- vector(mode="numeric",length=nSims) # pvalue for bootstrap
  pValChi2 <- vector(mode="numeric",length=nSims) # pvalue for chi sq

  for (isim in 1:nSims) {


    # simulate data set
    data <- simulate_ar2(alpha=0,beta=beta,sigma,c(rho1,rho2),nT)
    # fit under the null and alternative
    null <- fit_ar2_opt(data,rho=c(rho1,rho2),hypothesis="null")
    alt <- fit_ar2_opt(data,rho=c(rho1,rho2),hypothesis="alt")

    # likelihood statistic
    LRstatObs <- -2*(null$likelihood-alt$likelihood)
    # pvalue using chi square approximation
    pValChi2[isim] <- 1-pchisq(LRstatObs,1) # uses distributional theory



    # bootstrapping in parallel
    bootStats <- foreach::foreach(iboot = 2:nBootSims,.combine='c',.packages = "arfit") %dopar% {
    # simulate under Null
      bootdata <- simulate_ar2(alpha=null$betaEst,beta=0,null$sigmaEst,null$rhoEst,nT)
      # fit under null and alt
      nullBoot <- fit_ar2_opt(bootdata,rho=null$rhoEst,hypothesis="null")
      altBoot <- fit_ar2_opt(bootdata,rho=null$rhoEst,hypothesis="alt")
      # statisicic
      LR[iboot-1] <- -2*(nullBoot$likelihood-altBoot$likelihood)
    } # end bootstrap



    # now we can calculate the p-value based on the bootstrapping
    LRstat <- c(LRstatObs,bootStats)
    pVal_boot[isim] <- sum(LRstat >= LRstat[1])/nBootSims
    #print(paste0("pval_chi = ",pValChi2[isim],". pval_boot = ",pVal_boot[isim]))
  } # end sims
  # summarises the proportion of times the test rejects the null
  pValue <- sum(pVal_boot <= 0.05)/nSims
  pvChi <- sum(pValChi2 <= 0.05) /nSims

  #parallel::stopCluster(cl)
  return(list(pvChi=pvChi,pValue=pValue))

}


