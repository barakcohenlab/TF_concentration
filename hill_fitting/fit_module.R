library(boot)



AIC_from_fit <- function(nls_fit){
  res <- nls_fit$fvec
  # RSS
  RSS <- sum(res^2)
  # Sample size
  n <- length(res)
  # Number of parameters
  k <- length(nls_fit$par)
  AIC_value <- n * log(RSS/n) + 2 * k
  BIC_value <- n * log(RSS/n) + k * log(n)
  
  return(AIC_value)
  
}

fit_function_order_loop_coefstd <- function(conc_tf,lib_3rMean,r1_hsp_matrix,r2_hsp_matrix,r3_hsp_matrix,Kd_start_,n_start,lower_Kd,upper_Kd,range_x,seq_name,printName=F){
  #conc_tf <-conc_foxa1*conc_hnf4a # conc_foxa1 # conc_dox
  if(printName){
  #  print(seq_name)
  }
  observed_mean_tmp <- unlist(lib_3rMean[seq_name,])
  observed_mean <- observed_mean_tmp[!is.na(observed_mean_tmp)]
  conc_tf_noNA <-conc_tf[!is.na(observed_mean_tmp)]
  
  observed_tmp <- c(as.numeric(r1_hsp_matrix[seq_name,]),as.numeric(r2_hsp_matrix[seq_name,]),as.numeric(r3_hsp_matrix[seq_name,]))
  observed <- observed_tmp[!is.na(observed_tmp)]
  x <- rep(conc_tf,3)[!is.na(observed_tmp)]
  
  lower_basal = min(observed)
  
  hill_equation <- function(parS,xx){parS$basal + (parS$Fmax-parS$basal)/(1+(parS$Kd^(parS$n))/(xx^parS$n))}
  
  residFun <- function(p, observed, xx){((observed)^1-(hill_equation(p, xx)^1))} ## or log

  lower_basal = min(observed)
  if(seq_name == "v2_basal_hsp68"){
    Kd_start <- 0.1
  }else{
    Kd_start <- Kd_start_
  }

  ### r1 r2 r3
  observed_1c_tmp <- as.numeric(r1_hsp_matrix[seq_name,])
  observed_1c <- observed_1c_tmp[!is.na(observed_1c_tmp)]
  x_1c <- conc_tf[!is.na(observed_1c_tmp)]
  nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Fmax=observed_mean[length(observed_mean)],basal=observed_mean[1]), fn=residFun, observed = observed_1c, xx=x_1c, control = nls.lm.control(nprint=0,ftol=0.05),lower=c(lower_Kd,0.1,lower_basal,lower_basal),upper = c(upper_Kd,10,Inf,Inf))
  #  points(x_1c, observed_1c,  pch=16, cex=1,col="purple")
  #  lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col="purple")
  r1_coef <- as.list(coef(nls_sol_order))
  
  observed_1c_tmp <- as.numeric(r2_hsp_matrix[seq_name,])
  observed_1c <- observed_1c_tmp[!is.na(observed_1c_tmp)]
  x_1c <- conc_tf[!is.na(observed_1c_tmp)]
  nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Fmax=observed_mean[length(observed_mean)],basal=observed_mean[1]), fn=residFun, observed = observed_1c, xx=x_1c, control = nls.lm.control(nprint=0,ftol=0.05),lower=c(lower_Kd,0.1,lower_basal,lower_basal),upper = c(upper_Kd,10,Inf,Inf))
  #  points(x_1c, observed_1c,  pch=16, cex=1,col="green")
  #  lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col="green")
  r2_coef <- as.list(coef(nls_sol_order))
  
  observed_1c_tmp <- as.numeric(r3_hsp_matrix[seq_name,])
  observed_1c <- observed_1c_tmp[!is.na(observed_1c_tmp)]
  x_1c <- conc_tf[!is.na(observed_1c_tmp)]
  nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Fmax=observed_mean[length(observed_mean)],basal=observed_mean[1]), fn=residFun, observed = observed_1c, xx=x_1c, control = nls.lm.control(nprint=0,ftol=0.05),lower=c(lower_Kd,0.1,lower_basal,lower_basal),upper = c(upper_Kd,10,Inf,Inf))
  # points(x_1c, observed_1c,  pch=16, cex=1,col="blue")
  # lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col="blue")
  r3_coef <- as.list(coef(nls_sol_order))
  
  par_range_low <- c(min(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd,na.rm = TRUE),
                 min(r1_coef$n,r2_coef$n,r3_coef$n,na.rm = TRUE),
                 min(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax,na.rm = TRUE),
                 min(r1_coef$basal,r2_coef$basal,r3_coef$basal,na.rm = TRUE))
  par_range_low <- c(lower_Kd,0.1,lower_basal,lower_basal)
  
  par_range_high<- c(max(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd,na.rm = TRUE),
                     max(r1_coef$n,r2_coef$n,r3_coef$n,na.rm = TRUE),
                     max(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax,na.rm = TRUE),
                     max(r1_coef$basal,r2_coef$basal,r3_coef$basal,na.rm = TRUE))
  par_range_high <- c(upper_Kd,10,Inf,Inf)
  
  # par_start_0 <- list(Kd=median(c(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd),na.rm = TRUE),
  #                   n=median(c(r1_coef$n,r2_coef$n,r3_coef$n),na.rm = TRUE),
  #                   Fmax=median(c(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax),na.rm = TRUE),
  #                   basal=median(c(r1_coef$basal,r2_coef$basal,r3_coef$basal),na.rm = TRUE))
  # 
  Kd_start_test <- median(c(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd),na.rm = TRUE)/(c(0.2,0.5,1,2,5))
  r2_best <- -Inf
  best_i <- median(c(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd),na.rm = TRUE)
  for (i in Kd_start_test) {
    par_start <- list(Kd=i,
                      n=median(c(r1_coef$n,r2_coef$n,r3_coef$n),na.rm = TRUE),
                      Fmax=median(c(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax),na.rm = TRUE),
                      basal=median(c(r1_coef$basal,r2_coef$basal,r3_coef$basal),na.rm = TRUE))
    nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0,ftol=0.05),lower=par_range_low,upper = par_range_high)
    tmp_r2 <- cor(observed_mean, hill_equation(as.list(coef(nls_sol_order)), conc_tf_noNA),use="pairwise.complete.obs")^2
    #tmp_r2 <- AIC_from_fit(nls_sol_order)
    if((!is.na(tmp_r2)) && tmp_r2>r2_best){best_i <- i;r2_best <- tmp_r2}
  }
  
  par_start <- list(Kd=best_i,
                      n=median(c(r1_coef$n,r2_coef$n,r3_coef$n),na.rm = TRUE),
                      Fmax=median(c(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax),na.rm = TRUE),
                      basal=median(c(r1_coef$basal,r2_coef$basal,r3_coef$basal),na.rm = TRUE))
  
  #nls_sol_order <- nls.lm(par = list(Kd=1.5/(kd_score_F_table[seq_name,]$kas*kd_score_H_table[seq_name,]$kas),n=1,Fmax=observed_mean[8],basal=observed_mean[1]), fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0),lower=c(0,0.1,lower_basal,lower_basal),upper = c(Inf,10,Inf,Inf))
  #nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Fmax=observed_mean[length(observed_mean)],basal=observed_mean[1]), fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0,ftol=0.05),lower=c(lower_Kd,0.1,lower_basal,lower_basal),upper = c(upper_Kd,10,Inf,Inf))
  nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0,ftol=0.05),lower=par_range_low,upper = par_range_high)
  
return_aic <- AIC_from_fit(nls_sol_order)
  range_x <- range_x
 print(seq_name)
 plot(c(range_x,conc_tf_noNA), c(hill_equation(as.list(coef(nls_sol_order)), range_x),observed_mean), col="white",type="l",log="x", xlab="[TF]", ylab = "Activity",ylim=c(min(observed),max(observed)))
 points(conc_tf_noNA, observed_mean,  pch=16, cex=1)
 lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col=2)

  to_return_x <- cor(observed, hill_equation(as.list(coef(nls_sol_order)), x))^2
  to_return_x <- cor(observed_mean, hill_equation(as.list(coef(nls_sol_order)), conc_tf_noNA))^2
  
  to_return_coef <- as.list(coef(nls_sol_order))
  a <- to_return_coef$Fmax
  b <- to_return_coef$basal
  to_return_coef$regP <- log2(a/b)

  result <- try(vcov(nls_sol_order), silent=TRUE)
  if (inherits(result, "try-error")) {
    #stderr <- lapply(as.list(to_return_coef), function(x) replace(x, TRUE, NA)) 
    #stderr$regP <- NA
    pvalue <- lapply(as.list(to_return_coef), function(x) replace(x, TRUE, NA)) 
  boot_fun <- function(observed,index) {
    nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed[index], xx=x[index], control = nls.lm.control(nprint=0,ftol=0.05),lower=par_range_low,upper = par_range_high)
    boot_coef <- as.list(coef(nls_sol_order))
    boot_coef$regP <- log2(boot_coef$Fmax/boot_coef$basal)
    as.numeric(boot_coef)
  }
  result <- boot(data = observed, statistic = boot_fun, R = 500)
  stderr <- as.list(apply(result$t, 2, sd))
  names(stderr) <- c("Kd","n","Fmax","basal","regP")

  }else{
    stderr <- as.list(summary(nls_sol_order)$coefficients[,"Std. Error"])
    pvalue <- as.list(summary(nls_sol_order)$coefficients[,"Pr(>|t|)"])
    SE_a <- stderr$Fmax
    SE_b <- stderr$basal
    cov_ab <- vcov(nls_sol_order)["Fmax","basal"]
    # Compute variance and SE
    var_theta <- (SE_a/a)^2 + (SE_b/b)^2 - 2 * cov_ab / (a * b)
    stderr$regP <- sqrt(var_theta)/log(2)
    
  }
  bio_stderr <- list(Kd=sd(c(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd)),n=sd(c(r1_coef$n,r2_coef$n,r3_coef$n)),Fmax=sd(c(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax)),basal=sd(c(r1_coef$basal,r2_coef$basal,r3_coef$basal)),regP=sd(c(log2(r1_coef$Fmax/r1_coef$basal),log2(r2_coef$Fmax/r2_coef$basal),log2(r3_coef$Fmax/r3_coef$basal))))
  
    #list(Kd=std.error(c(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd,to_return_coef$Kd)),n=std.error(c(r1_coef$n,r2_coef$n,r3_coef$n,to_return_coef$n)),Fmax=std.error(c(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax,to_return_coef$Fmax)),basal=std.error(c(r1_coef$basal,r2_coef$basal,r3_coef$basal,to_return_coef$basal)))
  
  return(list(r2=to_return_x,coef=to_return_coef,std=stderr,pvalue=pvalue,coef1=r1_coef,coef3=r3_coef,coef2=r2_coef,biostd=bio_stderr,aic=return_aic))

}


fit_2r_function_order_loop_coefstd <- function(conc_tf,lib_3rMean,r1_hsp_matrix,r2_hsp_matrix,Kd_start_,n_start,lower_Kd,upper_Kd,range_x,seq_name,printName=F){
  #conc_tf <-conc_foxa1*conc_hnf4a # conc_foxa1 # conc_dox # r3_hsp_matrix
  if(printName){
      print(seq_name)
  }
  observed_mean_tmp <- unlist(lib_3rMean[seq_name,])
  observed_mean <- observed_mean_tmp[!is.na(observed_mean_tmp)]
  conc_tf_noNA <-conc_tf[!is.na(observed_mean_tmp)]
  
  observed_tmp <- c(as.numeric(r1_hsp_matrix[seq_name,]),as.numeric(r2_hsp_matrix[seq_name,]))
  observed <- observed_tmp[!is.na(observed_tmp)]
  x <- rep(conc_tf,2)[!is.na(observed_tmp)]
  
  lower_basal = min(observed)
  
  hill_equation <- function(parS,xx){parS$basal + (parS$Fmax-parS$basal)/(1+(parS$Kd^(parS$n))/(xx^parS$n))}
  
  residFun <- function(p, observed, xx){((observed)^1-(hill_equation(p, xx)^1))} ## or log
  
  lower_basal = min(observed)
  if(seq_name == "v2_basal_hsp68"){
    Kd_start <- 0.1
  }else{
    Kd_start <- Kd_start_
  }
  

  ### r1 r2 r3
  observed_1c_tmp <- as.numeric(r1_hsp_matrix[seq_name,])
  observed_1c <- observed_1c_tmp[!is.na(observed_1c_tmp)]
  x_1c <- conc_tf[!is.na(observed_1c_tmp)]
  nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Fmax=observed_mean[1],basal=observed_mean[length(observed_mean)]), fn=residFun, observed = observed_1c, xx=x_1c, control = nls.lm.control(nprint=0,ftol=0.05),lower=c(lower_Kd,0.1,lower_basal,lower_basal),upper = c(upper_Kd,10,Inf,Inf))
  plot(c(range_x,conc_tf_noNA), c(hill_equation(as.list(coef(nls_sol_order)), range_x),observed_mean), col="white",type="l",log="x", xlab="[TF]", ylab = "Activity",ylim=c(min(observed),max(observed)))
  points(x_1c, observed_1c,  pch=16, cex=1,col="purple")
    lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col="purple")
  r1_coef <- as.list(coef(nls_sol_order))
  
  observed_1c_tmp <- as.numeric(r2_hsp_matrix[seq_name,])
  observed_1c <- observed_1c_tmp[!is.na(observed_1c_tmp)]
  x_1c <- conc_tf[!is.na(observed_1c_tmp)]
  nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Fmax=observed_mean[1],basal=observed_mean[length(observed_mean)]), fn=residFun, observed = observed_1c, xx=x_1c, control = nls.lm.control(nprint=0,ftol=0.05),lower=c(lower_Kd,0.1,lower_basal,lower_basal),upper = c(upper_Kd,10,Inf,Inf))
    points(x_1c, observed_1c,  pch=16, cex=1,col="green")
    lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col="green")
  r2_coef <- as.list(coef(nls_sol_order))
  
  
    par_range_low <- c(lower_Kd,0.1,lower_basal,lower_basal)
  
    par_range_high <- c(upper_Kd,10,Inf,Inf)
  
    Kd_start_test <- median(c(r1_coef$Kd,r2_coef$Kd),na.rm = TRUE)/(c(0.2,0.5,1,2,5))
    r2_best <- -Inf
    best_i <- median(c(r1_coef$Kd,r2_coef$Kd),na.rm = TRUE)
    for (i in Kd_start_test) {
      par_start <- list(Kd=i,
                        n=median(c(r1_coef$n,r2_coef$n),na.rm = TRUE),
                        Fmax=median(c(r1_coef$Fmax,r2_coef$Fmax),na.rm = TRUE),
                        basal=median(c(r1_coef$basal,r2_coef$basal),na.rm = TRUE))
      nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0,ftol=0.05),lower=par_range_low,upper = par_range_high)
      tmp_r2 <- cor(observed_mean, hill_equation(as.list(coef(nls_sol_order)), conc_tf_noNA),use="pairwise.complete.obs")^2
      #tmp_r2 <- AIC_from_fit(nls_sol_order)
      if((!is.na(tmp_r2)) && tmp_r2>r2_best){best_i <- i;r2_best <- tmp_r2}
    }
    
  par_start <- list(Kd=best_i,
                    n=median(c(r1_coef$n,r2_coef$n),na.rm = TRUE),
                    Fmax=median(c(r1_coef$Fmax,r2_coef$Fmax),na.rm = TRUE),
                    basal=median(c(r1_coef$basal,r2_coef$basal),na.rm = TRUE))
  
  #nls_sol_order <- nls.lm(par = list(Kd=1.5/(kd_score_F_table[seq_name,]$kas*kd_score_H_table[seq_name,]$kas),n=1,Fmax=observed_mean[8],basal=observed_mean[1]), fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0),lower=c(0,0.1,lower_basal,lower_basal),upper = c(Inf,10,Inf,Inf))
  #nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Fmax=observed_mean[length(observed_mean)],basal=observed_mean[1]), fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0,ftol=0.05),lower=c(lower_Kd,0.1,lower_basal,lower_basal),upper = c(upper_Kd,10,Inf,Inf))
  nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0,ftol=0.05),lower=par_range_low,upper = par_range_high)
  
  
  #plot(c(range_x,conc_tf_noNA), c(hill_equation(as.list(coef(nls_sol_order)), range_x),observed_mean), col="white",type="l",log="x", xlab="[TF]", ylab = "Activity")
  points(conc_tf_noNA, observed_mean,  pch=16, cex=1)
  lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col=2)
  
  to_return_x <- cor(observed, hill_equation(as.list(coef(nls_sol_order)), x))^2
  to_return_x <- cor(observed_mean, hill_equation(as.list(coef(nls_sol_order)), conc_tf_noNA))^2
  
  to_return_coef <- as.list(coef(nls_sol_order))
  a <- to_return_coef$Fmax
  b <- to_return_coef$basal
  to_return_coef$regP <- log2((1+a)/(b+1))
  
  result <- try(vcov(nls_sol_order), silent=TRUE)
  if (inherits(result, "try-error")) {
    #stderr <- lapply(as.list(to_return_coef), function(x) replace(x, TRUE, NA)) 
    #stderr$regP <- NA
    pvalue <- lapply(as.list(to_return_coef), function(x) replace(x, TRUE, NA)) 
    boot_fun <- function(observed,index) {
      nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed[index], xx=x[index], control = nls.lm.control(nprint=0,ftol=0.05),lower=par_range_low,upper = par_range_high)
      boot_coef <- as.list(coef(nls_sol_order))
      boot_coef$regP <- log2(boot_coef$Fmax/boot_coef$basal)
      as.numeric(boot_coef)
    }
    result <- boot(data = observed, statistic = boot_fun, R = 500)
    stderr <- as.list(apply(result$t, 2, sd))
    names(stderr) <- c("Kd","n","Fmax","basal","regP")
    
  }else{
    stderr <- as.list(summary(nls_sol_order)$coefficients[,"Std. Error"])
    pvalue <- as.list(summary(nls_sol_order)$coefficients[,"Pr(>|t|)"])
    SE_a <- stderr$Fmax
    SE_b <- stderr$basal
    cov_ab <- vcov(nls_sol_order)["Fmax","basal"]
    # Compute variance and SE
    var_theta <- (SE_a/a)^2 + (SE_b/b)^2 - 2 * cov_ab / (a * b)
    stderr$regP <- sqrt(var_theta)/log(2)
    
  }
  bio_stderr <- list(Kd=sd(c(r1_coef$Kd,r2_coef$Kd)),n=sd(c(r1_coef$n,r2_coef$n)),Fmax=sd(c(r1_coef$Fmax,r2_coef$Fmax)),basal=sd(c(r1_coef$basal,r2_coef$basal)),regP=sd(c(log2(r1_coef$Fmax/r1_coef$basal),log2(r2_coef$Fmax/r2_coef$basal))))
  
  #list(Kd=std.error(c(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd,to_return_coef$Kd)),n=std.error(c(r1_coef$n,r2_coef$n,r3_coef$n,to_return_coef$n)),Fmax=std.error(c(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax,to_return_coef$Fmax)),basal=std.error(c(r1_coef$basal,r2_coef$basal,r3_coef$basal,to_return_coef$basal)))
  
  return(list(r2=to_return_x,coef=to_return_coef,std=stderr,pvalue=pvalue,coef1=r1_coef,coef2=r2_coef,biostd=bio_stderr))
  
}



fit_function_order_loop <- function(conc_tf,lib_3rMean,r1_hsp_matrix,r2_hsp_matrix,r3_hsp_matrix,Kd_start_,n_start,lower_Kd,upper_Kd,range_x,seq_name,printName=F){
  #conc_tf <-conc_foxa1*conc_hnf4a # conc_foxa1 # conc_dox
  if(printName){
  #  print(seq_name)
  }
  observed_mean_tmp <- unlist(lib_3rMean[seq_name,])
  observed_mean <- observed_mean_tmp[!is.na(observed_mean_tmp)]
  conc_tf_noNA <-conc_tf[!is.na(observed_mean_tmp)]
  
  observed_tmp <- c(as.numeric(r1_hsp_matrix[seq_name,]),as.numeric(r2_hsp_matrix[seq_name,]),as.numeric(r3_hsp_matrix[seq_name,]))
  observed <- observed_tmp[!is.na(observed_tmp)]
  x <- rep(conc_tf,3)[!is.na(observed_tmp)]
  
  #x2 <- rep(conc_tf2,3)
  lower_basal = min(observed)
  
  hill_equation <- function(parS,xx){parS$basal + (parS$Fmax-parS$basal)/(1+(parS$Kd^(parS$n))/(xx^parS$n))}
  
  residFun <- function(p, observed, xx){((observed)^1-(hill_equation(p, xx)^1))} ## or log
  # 1.5/(kd_score_F_table[seq_name,]$kas*kd_score_H_table[seq_name,]$kas)
  lower_basal = 0
  lower_basal = min(observed)
  if(seq_name == "v2_basal_hsp68"){
    Kd_start <- 0.1
  }else{
    Kd_start <- Kd_start_
  }
  
  ### r1 r2 r3
  observed_1c_tmp <- as.numeric(r1_hsp_matrix[seq_name,])
  observed_1c <- observed_1c_tmp[!is.na(observed_1c_tmp)]
  x_1c <- conc_tf[!is.na(observed_1c_tmp)]
  nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Fmax=observed_mean[length(observed_mean)],basal=observed_mean[1]), fn=residFun, observed = observed_1c, xx=x_1c, control = nls.lm.control(nprint=0,ftol=0.05),lower=c(lower_Kd,0.1,lower_basal,lower_basal),upper = c(upper_Kd,10,Inf,Inf))
  #  points(x_1c, observed_1c,  pch=16, cex=1,col="purple")
  #  lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col="purple")
  r1_coef <- as.list(coef(nls_sol_order))
  
  observed_1c_tmp <- as.numeric(r2_hsp_matrix[seq_name,])
  observed_1c <- observed_1c_tmp[!is.na(observed_1c_tmp)]
  x_1c <- conc_tf[!is.na(observed_1c_tmp)]
  nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Fmax=observed_mean[length(observed_mean)],basal=observed_mean[1]), fn=residFun, observed = observed_1c, xx=x_1c, control = nls.lm.control(nprint=0,ftol=0.05),lower=c(lower_Kd,0.1,lower_basal,lower_basal),upper = c(upper_Kd,10,Inf,Inf))
  #  points(x_1c, observed_1c,  pch=16, cex=1,col="green")
  #  lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col="green")
  r2_coef <- as.list(coef(nls_sol_order))
  
  observed_1c_tmp <- as.numeric(r3_hsp_matrix[seq_name,])
  observed_1c <- observed_1c_tmp[!is.na(observed_1c_tmp)]
  x_1c <- conc_tf[!is.na(observed_1c_tmp)]
  nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Fmax=observed_mean[length(observed_mean)],basal=observed_mean[1]), fn=residFun, observed = observed_1c, xx=x_1c, control = nls.lm.control(nprint=0,ftol=0.05),lower=c(lower_Kd,0.1,lower_basal,lower_basal),upper = c(upper_Kd,10,Inf,Inf))
  # points(x_1c, observed_1c,  pch=16, cex=1,col="blue")
  # lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col="blue")
  r3_coef <- as.list(coef(nls_sol_order))
  
  
  par_range_low <- c(min(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd,na.rm = TRUE),
                     min(r1_coef$n,r2_coef$n,r3_coef$n,na.rm = TRUE),
                     min(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax,na.rm = TRUE),
                     min(r1_coef$basal,r2_coef$basal,r3_coef$basal,na.rm = TRUE))
  
  par_range_high<- c(max(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd,na.rm = TRUE),
                     max(r1_coef$n,r2_coef$n,r3_coef$n,na.rm = TRUE),
                     max(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax,na.rm = TRUE),
                     max(r1_coef$basal,r2_coef$basal,r3_coef$basal,na.rm = TRUE))
  
  par_start <- list(Kd=median(c(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd),na.rm = TRUE),
                    n=median(c(r1_coef$n,r2_coef$n,r3_coef$n),na.rm = TRUE),
                    Fmax=median(c(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax),na.rm = TRUE),
                    basal=median(c(r1_coef$basal,r2_coef$basal,r3_coef$basal),na.rm = TRUE))
  
  #nls_sol_order <- nls.lm(par = list(Kd=1.5/(kd_score_F_table[seq_name,]$kas*kd_score_H_table[seq_name,]$kas),n=1,Fmax=observed_mean[8],basal=observed_mean[1]), fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0),lower=c(0,0.1,lower_basal,lower_basal),upper = c(Inf,10,Inf,Inf))
  #nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Fmax=observed_mean[length(observed_mean)],basal=observed_mean[1]), fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0,ftol=0.05),lower=c(lower_Kd,0.1,lower_basal,lower_basal),upper = c(upper_Kd,10,Inf,Inf))
  nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0,ftol=0.05),lower=par_range_low,upper = par_range_high)
  

  
  range_x <- range_x
  # plot(c(range_x,conc_tf_noNA), c(hill_equation(as.list(coef(nls_sol_order)), range_x),observed_mean), col="white",type="l",log="x", xlab="[TF]", ylab = "Activity")
  # points(conc_tf_noNA, observed_mean,  pch=16, cex=1)
  # lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col=2)
  # 
  to_return_x <- cor(observed, hill_equation(as.list(coef(nls_sol_order)), x))^2
  to_return_x <- cor(observed_mean, hill_equation(as.list(coef(nls_sol_order)), conc_tf_noNA))^2
  
  to_return_coef <- as.list(coef(nls_sol_order))
  a <- to_return_coef$Fmax
  b <- to_return_coef$basal
  to_return_coef$regP <- log2(a/b)
  

  stderr <- list(Kd=sd(c(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd,to_return_coef$Kd)),n=sd(c(r1_coef$n,r2_coef$n,r3_coef$n,to_return_coef$n)),Fmax=sd(c(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax,to_return_coef$Fmax)),basal=sd(c(r1_coef$basal,r2_coef$basal,r3_coef$basal,to_return_coef$basal)))
  
  return(list(r2=to_return_x,coef=to_return_coef,coef1=r1_coef,coef3=r3_coef,coef2=r2_coef,std=as.list(stderr)))

}





multi_hill_loop_coefstd <- function(conc_tf,lib_3rMean,r1_hsp_matrix,r2_hsp_matrix,r3_hsp_matrix,Kd_start_,n_start,lower_Kd,upper_Kd,range_x,seq_name,printName=F){
  #conc_tf <-conc_foxa1*conc_hnf4a # conc_foxa1 # conc_dox
  if(printName){
    #  print(seq_name)
  }
  observed_mean_tmp <- unlist(lib_3rMean[seq_name,])
  observed_mean <- observed_mean_tmp[!is.na(observed_mean_tmp)]
  conc_tf_noNA <-conc_tf[!is.na(observed_mean_tmp)]
  
  observed_tmp <- c(as.numeric(r1_hsp_matrix[seq_name,]),as.numeric(r2_hsp_matrix[seq_name,]),as.numeric(r3_hsp_matrix[seq_name,]))
  observed <- observed_tmp[!is.na(observed_tmp)]
  x <- rep(conc_tf,3)[!is.na(observed_tmp)]
  
  lower_basal = min(observed)
  
  hill_equation <- function(parS,xx){parS$basal + (parS$Fmax-parS$basal)*(1/(1+(parS$Kd^(parS$n))/(xx^parS$n))+1/(1+(parS$Kd2^(parS$n2))/(xx^parS$n2)))}
  
  residFun <- function(p, observed, xx){((observed)^1-(hill_equation(p, xx)^1))} ## or log
  
  lower_basal = min(observed)
  if(seq_name == "v2_basal_hsp68"){
    Kd_start <- 0.1
  }else{
    Kd_start <- Kd_start_
  }
  
  ### r1 r2 r3
  rep_lower <- c(lower_Kd,0.1,lower_Kd,0.1,lower_basal,lower_basal)
  rep_upper = c(upper_Kd,10,upper_Kd,10,Inf,Inf)
  observed_1c_tmp <- as.numeric(r1_hsp_matrix[seq_name,])
  observed_1c <- observed_1c_tmp[!is.na(observed_1c_tmp)]
  x_1c <- conc_tf[!is.na(observed_1c_tmp)]
  nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Kd2=Kd_start*2,n2=1/n_start,Fmax=observed_mean[length(observed_mean)],basal=observed_mean[1]), fn=residFun, observed = observed_1c, xx=x_1c, control = nls.lm.control(nprint=0,ftol=0.05),lower=rep_lower,upper = rep_upper)
  #  points(x_1c, observed_1c,  pch=16, cex=1,col="purple")
  #  lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col="purple")
  r1_coef <- as.list(coef(nls_sol_order))
  
  observed_1c_tmp <- as.numeric(r2_hsp_matrix[seq_name,])
  observed_1c <- observed_1c_tmp[!is.na(observed_1c_tmp)]
  x_1c <- conc_tf[!is.na(observed_1c_tmp)]
  nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Kd2=Kd_start*2,n2=1/n_start,Fmax=observed_mean[length(observed_mean)],basal=observed_mean[1]), fn=residFun, observed = observed_1c, xx=x_1c, control = nls.lm.control(nprint=0,ftol=0.05),lower=rep_lower,upper = rep_upper)
  #  points(x_1c, observed_1c,  pch=16, cex=1,col="green")
  #  lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col="green")
  r2_coef <- as.list(coef(nls_sol_order))
  
  observed_1c_tmp <- as.numeric(r3_hsp_matrix[seq_name,])
  observed_1c <- observed_1c_tmp[!is.na(observed_1c_tmp)]
  x_1c <- conc_tf[!is.na(observed_1c_tmp)]
  nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Kd2=Kd_start*2,n2=1/n_start,Fmax=observed_mean[length(observed_mean)],basal=observed_mean[1]), fn=residFun, observed = observed_1c, xx=x_1c, control = nls.lm.control(nprint=0,ftol=0.05),lower=rep_lower,upper = rep_upper)
  # points(x_1c, observed_1c,  pch=16, cex=1,col="blue")
  # lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col="blue")
  r3_coef <- as.list(coef(nls_sol_order))
  
    par_range_low <- c(lower_Kd,0.1,lower_Kd,0.1,lower_basal,lower_basal)
  
   par_range_high <- c(upper_Kd,10,upper_Kd,1,Inf,Inf)
  
  par_start <- list(Kd=median(c(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd),na.rm = TRUE),
                    n=median(c(r1_coef$n,r2_coef$n,r3_coef$n),na.rm = TRUE),
                    Kd2=median(c(r1_coef$Kd2,r2_coef$Kd2,r3_coef$Kd2),na.rm = TRUE),
                    n2=median(c(r1_coef$n2,r2_coef$n2,r3_coef$n2),na.rm = TRUE),
                    Fmax=median(c(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax),na.rm = TRUE),
                    basal=median(c(r1_coef$basal,r2_coef$basal,r3_coef$basal),na.rm = TRUE))
  
  #nls_sol_order <- nls.lm(par = list(Kd=1.5/(kd_score_F_table[seq_name,]$kas*kd_score_H_table[seq_name,]$kas),n=1,Fmax=observed_mean[8],basal=observed_mean[1]), fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0),lower=c(0,0.1,lower_basal,lower_basal),upper = c(Inf,10,Inf,Inf))
  #nls_sol_order <- nls.lm(par = list(Kd=Kd_start,n=n_start,Kd2=Kd_start,n2=n_start,Fmax=observed_mean[length(observed_mean)],basal=observed_mean[1]), fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0,ftol=0.05),lower=c(lower_Kd,0.1,lower_basal,lower_basal),upper = c(upper_Kd,10,Inf,Inf))
  nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0,ftol=0.001),lower=par_range_low,upper = par_range_high)
  
  
  range_x <- range_x
  print(seq_name)
  plot(c(range_x,conc_tf_noNA), c(hill_equation(as.list(coef(nls_sol_order)), range_x),observed_mean), col="white",type="l",log="x", xlab="[TF]", ylab = "Activity",ylim=c(min(observed),max(observed)))
  points(conc_tf_noNA, observed_mean,  pch=16, cex=1)
  lines(range_x, hill_equation(as.list(coef(nls_sol_order)), range_x),col=2)
  
  to_return_x <- cor(observed, hill_equation(as.list(coef(nls_sol_order)), x))^2
  to_return_x <- cor(observed_mean, hill_equation(as.list(coef(nls_sol_order)), conc_tf_noNA))^2
  
  to_return_coef <- as.list(coef(nls_sol_order))
  a <- to_return_coef$Fmax
  b <- to_return_coef$basal
  to_return_coef$regP <- log2(a/b)
  
  result <- try(vcov(nls_sol_order), silent=TRUE)
  if (inherits(result, "try-error")) {
    #stderr <- lapply(as.list(to_return_coef), function(x) replace(x, TRUE, NA)) 
    #stderr$regP <- NA
    pvalue <- lapply(as.list(to_return_coef), function(x) replace(x, TRUE, NA)) 
    boot_fun <- function(observed,index) {
      nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed[index], xx=x[index], control = nls.lm.control(nprint=0,ftol=0.05),lower=par_range_low,upper = par_range_high)
      boot_coef <- as.list(coef(nls_sol_order))
      boot_coef$regP <- log2(boot_coef$Fmax/boot_coef$basal)
      as.numeric(boot_coef)
    }
    result <- boot(data = observed, statistic = boot_fun, R = 5)
    stderr <- as.list(apply(result$t, 2, sd))
    names(stderr) <- c("Kd","n","Kd2","n2","Fmax","basal","regP")
    
  }else{
    stderr <- as.list(summary(nls_sol_order)$coefficients[,"Std. Error"])
    pvalue <- as.list(summary(nls_sol_order)$coefficients[,"Pr(>|t|)"])
    SE_a <- stderr$Fmax
    SE_b <- stderr$basal
    cov_ab <- vcov(nls_sol_order)["Fmax","basal"]
    # Compute variance and SE
    var_theta <- (SE_a/a)^2 + (SE_b/b)^2 - 2 * cov_ab / (a * b)
    stderr$regP <- sqrt(var_theta)/log(2)
    
  }
  bio_stderr <- list(Kd=sd(c(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd)),n=sd(c(r1_coef$n,r2_coef$n,r3_coef$n)),Fmax=sd(c(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax)),basal=sd(c(r1_coef$basal,r2_coef$basal,r3_coef$basal)),regP=sd(c(log2(r1_coef$Fmax/r1_coef$basal),log2(r2_coef$Fmax/r2_coef$basal),log2(r3_coef$Fmax/r3_coef$basal))))
  
  #list(Kd=std.error(c(r1_coef$Kd,r2_coef$Kd,r3_coef$Kd,to_return_coef$Kd)),n=std.error(c(r1_coef$n,r2_coef$n,r3_coef$n,to_return_coef$n)),Fmax=std.error(c(r1_coef$Fmax,r2_coef$Fmax,r3_coef$Fmax,to_return_coef$Fmax)),basal=std.error(c(r1_coef$basal,r2_coef$basal,r3_coef$basal,to_return_coef$basal)))
  
  return(list(r2=to_return_x,coef=to_return_coef,std=stderr,pvalue=pvalue,coef1=r1_coef,coef3=r3_coef,coef2=r2_coef,biostd=bio_stderr))
  
}




thermo_fit_1 <- function(conc_tf,lib_3rMean,r1_hsp_matrix,r2_hsp_matrix,r3_hsp_matrix,Kd_tf,Kd_start_,c_start=1.5,lower_Kd=0,upper_Kd=Inf,range_x,seq_name,printName=F){
  if(printName){
    #  print(seq_name)
  }
  observed_mean_tmp <- unlist(lib_3rMean[seq_name,])
  observed_mean <- observed_mean_tmp[!is.na(observed_mean_tmp)]
  conc_tf_noNA <-conc_tf[!is.na(observed_mean_tmp)]
  
  observed_tmp <- c(as.numeric(r1_hsp_matrix[seq_name,]),as.numeric(r2_hsp_matrix[seq_name,]),as.numeric(r3_hsp_matrix[seq_name,]))
  observed <- observed_tmp[!is.na(observed_tmp)]
  x <- rep(conc_tf,3)[!is.na(observed_tmp)]
  
  #  Kd_tf <- Kd_tf
  
  #thermo_equation <- function(parS,xx){Fmax1 * parS$basal * (1+parS$c1*Kd_tf*xx) * (1+parS$c2*parS$Kdr*xx) / ( (1+parS$basal)*(1+Kd_tf*xx)*(1+parS$Kdr*xx) )}
  thermo_equation <- function(parS,xx){ parS$basal * (1+parS$c1*Kd_tf*xx) * (1+parS$c2*parS$Kdr*xx) / ( (1+Kd_tf*xx)*(1+parS$Kdr*xx) )}
  
  residFun <- function(p, observed, xx){((observed)^1-(thermo_equation(p, xx)^1))} ## or log
  
  lower_basal = 0
  # if(seq_name == "v2_basal_hsp68"){
  #   Kd_start <- 0.1
  # }else{
  #   Kd_start <- Kd_start_
  # }
  # Kd_start  <- kd_score_F_table[seq_name,]$kas
  brp_start <- lib_3rMean[seq_name,]$s5 #(1/(1-lib_3rMean[seq_name,]$s5/Fmax1)-1)
  par_range_low <- c(lower_basal,0,0,0)
  par_range_high <- c(Inf,Inf,Inf,Inf)
  
  kdr_start <- Kd_tf/(c(.5,1,2.5,5,10,20))
  r2_best <- 0
  best_i <- NA
  for (i in kdr_start) {
    par_start <- list(basal=brp_start,
                      c1=c_start,
                      c2=0.5,
                      Kdr= i) #min(i,0.0001))
    nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0),lower=par_range_low,upper = par_range_high)
    tmp_r2 <- cor(observed_mean, thermo_equation(as.list(coef(nls_sol_order)), conc_tf_noNA))^2
    if((!is.na(tmp_r2)) && tmp_r2>r2_best){best_i <- i}
    
  }
  par_start <- list(basal=brp_start,
                    c1=c_start,
                    c2=0.5,
                    Kdr=best_i)
  
  nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0),lower=par_range_low,upper = par_range_high)
  range_x <- range_x
  print(seq_name)
  plot(c(range_x,conc_tf_noNA), c(thermo_equation(as.list(coef(nls_sol_order)), range_x),observed_mean), col="white",type="l",log="x", xlab="[TF]", ylab = "Activity",ylim=c(min(observed),max(observed)))
  points(conc_tf_noNA, observed_mean,  pch=16, cex=1)
  lines(range_x, thermo_equation(as.list(coef(nls_sol_order)), range_x),col=2)
  
  #to_return_x <- cor(observed, thermo_equation(as.list(coef(nls_sol_order)), x))^2
  to_return_x <- cor(observed_mean, thermo_equation(as.list(coef(nls_sol_order)), conc_tf_noNA))^2
  
  to_return_coef <- as.list(coef(nls_sol_order))
  aic <- AIC_from_fit(nls_sol_order)
  return(list(r2=to_return_x,coef=to_return_coef,aic=aic))
  
}


thermo_fit <- function(conc_tf,lib_3rMean,r1_hsp_matrix,r2_hsp_matrix,r3_hsp_matrix,Kd_tf,Kd_start_,c_start=1.5,lower_Kd=0,upper_Kd=Inf,range_x,seq_name,printName=F){
  if(printName){
    #  print(seq_name)
  }
  observed_mean_tmp <- unlist(lib_3rMean[seq_name,])
  observed_mean <- observed_mean_tmp[!is.na(observed_mean_tmp)]
  conc_tf_noNA <-conc_tf[!is.na(observed_mean_tmp)]
  
  observed_tmp <- c(as.numeric(r1_hsp_matrix[seq_name,]),as.numeric(r2_hsp_matrix[seq_name,]),as.numeric(r3_hsp_matrix[seq_name,]))
  observed <- observed_tmp[!is.na(observed_tmp)]
  x <- rep(conc_tf,3)[!is.na(observed_tmp)]
  
  #  Kd_tf <- Kd_tf
  
  #thermo_equation <- function(parS,xx){Fmax1 * parS$basal * (1+parS$c1*Kd_tf*xx) * (1+parS$c2*parS$Kdr*xx) / ( (1+parS$basal)*(1+Kd_tf*xx)*(1+parS$Kdr*xx) )}
  thermo_equation <- function(parS,xx){ parS$basal * (1+parS$c1*parS$Kdtf*xx) * (1+parS$c2*parS$Kdr*xx) / ( (1+parS$Kdtf*xx)*(1+parS$Kdr*xx) )}
  
  residFun <- function(p, observed, xx){((observed)^1-(thermo_equation(p, xx)^1))} ## or log
  
  lower_basal = 0
  # if(seq_name == "v2_basal_hsp68"){
  #   Kd_start <- 0.1
  # }else{
  #   Kd_start <- Kd_start_
  # }
  # Kd_start  <- kd_score_F_table[seq_name,]$kas
  brp_start <- lib_3rMean[seq_name,]$s5 #(1/(1-lib_3rMean[seq_name,]$s5/Fmax1)-1)
  par_range_low <- c(lower_basal,0,0,0,0)
  par_range_high <- c(Inf,Inf,Inf,Inf,Inf)
  
  kdr_start <- Kd_tf/(c(.5,1,2.5,5,10,20))
  r2_best <- Inf
  best_i <- NA
  for (i in kdr_start) {
    par_start <- list(basal=brp_start,
                      c1=c_start,
                      c2=0.5,
                      Kdr= i,
                      Kdtf=Kd_start_) #min(i,0.0001))
    nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0),lower=par_range_low,upper = par_range_high)
    #tmp_r2 <- cor(observed_mean, thermo_equation(as.list(coef(nls_sol_order)), conc_tf_noNA))^2
    tmp_r2 <- AIC_from_fit(nls_sol_order)
    if((!is.na(tmp_r2)) && tmp_r2<r2_best){best_i <- i;r2_best <- tmp_r2}
    
  }
  kdtf_start <- Kd_start_ *(c(.02,.1,.5,1,5,20))
  r2_best <- Inf
  best_i2 <- NA
  for (i in kdtf_start) {
    par_start <- list(basal=brp_start,
                      c1=c_start,
                      c2=0.5,
                      Kdr= best_i,
                      Kdtf=i) #min(i,0.0001))
    nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0),lower=par_range_low,upper = par_range_high)
    #tmp_r2 <- cor(observed_mean, thermo_equation(as.list(coef(nls_sol_order)), conc_tf_noNA))^2
    tmp_r2 <- AIC_from_fit(nls_sol_order)
    if((!is.na(tmp_r2)) && tmp_r2<r2_best){best_i2 <- i;r2_best <- tmp_r2}
    
  }
  par_start <- list(basal=brp_start,
                    c1=c_start,
                    c2=0.5,
                    Kdr=best_i,
                    Kdtf=best_i2)
  
  nls_sol_order <- nls.lm(par = par_start, fn=residFun, observed = observed, xx=x, control = nls.lm.control(nprint=0),lower=par_range_low,upper = par_range_high)
  range_x <- range_x
  print(seq_name)
  plot(c(range_x,conc_tf_noNA), c(thermo_equation(as.list(coef(nls_sol_order)), range_x),observed_mean), col="white",type="l",log="x", xlab="[TF]", ylab = "Activity",ylim=c(min(observed),max(observed)))
  points(conc_tf_noNA, observed_mean,  pch=16, cex=1)
  lines(range_x, thermo_equation(as.list(coef(nls_sol_order)), range_x),col=2)
  
  #to_return_x <- cor(observed, thermo_equation(as.list(coef(nls_sol_order)), x))^2
  to_return_x <- cor(observed_mean, thermo_equation(as.list(coef(nls_sol_order)), conc_tf_noNA))^2
  
  to_return_coef <- as.list(coef(nls_sol_order))
  aic <- AIC_from_fit(nls_sol_order)
  return(list(r2=to_return_x,coef=to_return_coef,aic=aic))
  
}


