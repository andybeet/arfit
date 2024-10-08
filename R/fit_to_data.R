#' Run test on user data set
#'
#' Fits a linear model to data set assuming AR1 error
#'
#' @param dataSet Data frame. Two columns (x,y)
#' @param nBootSims Numeric scalar. Number of bootstrap samples to perform
#'
#'
#'@export

fit_to_data <- function(dataSet=NULL,nBootSims=999) {

  dataValidation <- check_data_validation(dataSet)
  dataSet <- dataValidation$dataSet
  missingValues <- dataValidation$missingValues

  nT <- nrow(dataSet)
  data <- data.frame(x=c(1:nT),y=dataSet$y)

  # fit under the null and alternative
  # null <- fit_ar1_grid(data,hypothesis ="null")
  # alt <- fit_ar1_grid(data,hypothesis="alt")
   null <- fit_ar1_opt(data,rho=0.5,hypothesis ="null")
   alt <- fit_ar1_opt(data,null$rhoEst,hypothesis="alt")
  # preallocate likelihood ratio statistic vector
  LRstat <- vector(mode="numeric",length=nBootSims+1)
  # LR stat for data
  LRstat[1] <- -2*(null$likelihood-alt$likelihood)
  print(paste0("LR stat = ",LRstat[1]))
  # pvalue using chi square approximation
  pValChi2 <- 1-pchisq(LRstat[1],1) # uses distributional theory

  # Perform bootstrapping
  for (iboot in 2:(nBootSims+1)) {
    # simulate under Null
    bootdata <- simulate_ar1(alpha=null$betaEst,beta=0,null$sigmaEst,null$rhoEst,nT,missingValues = missingValues)
    # fit under null and alt
    # nullBoot <- fit_ar1_grid(bootdata,hypothesis="null")
    # altBoot <- fit_ar1_grid(bootdata,hypothesis="alt")
     nullBoot <- fit_ar1_opt(bootdata,null$rhoEst,hypothesis="null")
     altBoot <- fit_ar1_opt(bootdata,null$rhoEst,hypothesis="alt")
    # statisicic
    LRstat[iboot] <- -2*(nullBoot$likelihood-altBoot$likelihood)
  } # end bootstrap

  # now we can calculate the p-value based on the bootstrapping
  pVal_boot <- sum(LRstat >= LRstat[1])/(nBootSims+1)
  print(paste0("pval_boot = ",pVal_boot))

  # plot fits under the null and alternative
  #png("Figure7.png",width=900,height=600,units="px")
  par(mai=c(1,1.5,0,0),oma=c(0,0,1,1))
  plot(dataSet$x,dataSet$y,type="l",xlab="Year",ylab="Response",
       cex.lab=2.5,cex.axis=2,lwd=2)
  lines(dataSet$x,rep(null$betaEst,nT),col="black",lty=2,lwd=2)
  lines(dataSet$x,alt$betaEst[1]+alt$betaEst[2]*c(1:nT),col="black",lty=3,lwd=2)
  #dev.off()
  return(list(null=null, alt=alt,pValChi2=pValChi2))
}
