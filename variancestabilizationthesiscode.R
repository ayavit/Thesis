# library for NB glm
library(MASS)

# library for boxplots
library(fields)





find_model.nb.data = function(curr_data,curr_model_names){
  theta <- tryCatch({
    a = optimize(nb.glm.theta.data, data=curr_data,curr_model_names = curr_model_names,
                 interval=c(0.049, 11), maximum=TRUE)   # 0.05 and 100 are the boundaries of the search area of theta - they can be changed!
    theta = as.numeric(a[1])
  }, warning = function(war) {
    a = optimize(nb.glm.theta.data, data=curr_data,curr_model_names = curr_model_names,
                 interval=c(0.049, 11), maximum=TRUE)  # 0.05 and 100 are the boundaries of the search area of theta - they can be changed!
    theta = as.numeric(a[1])
    return(theta)
  }, error = function(err) {
    print(paste("MY_ERROR find_model.nb.data:  ",err))
    theta = -1
    return(theta)
  }, finally = {
  }) # END tryCatch
  
  if(theta!=-1){
    model.nb = glm(curr_model_names,
                   family = negative.binomial(theta),data = curr_data,maxit = 1000)
  }else{model.nb = 0}
  
  return(model.nb)
}

# curr_model_names is the input of a "standard" glm - for examples y~. here is how you write it as a formula:
# model_formula = as.formula("y","~"," . -exposure + offset(log(exposure))"))
# model_formula = as.formula(colnames(data)[1],"~"," . -exposure + offset(log(exposure))"))
# 
# You can change the limits of the shape paramter by changing the values 0.05 and 100.
# Here is an example of calling the function:
# model.nb = find_model.nb.data(data,model_formula)


get_theta = function(model.nb){
  tmp = strsplit(model.nb$family[[1]],"")
  tmp2 = tmp[[1]][19:(length(tmp[[1]])-1)]
  theta = as.numeric(paste(tmp2,collapse = ""))
  return(theta)
}




# explanatory variable: a vector of length 25
# x_values <- seq(from = 0, to = 4.8, by = 0.2)

# values of explanatory variable (exp(mean gene expression level)) that fit the graph from the paper I found
#x_values <- c(0.1, 0.5, 1, 5, 10)

# values of explanatory variable (log(mean gene expression level)) that fit the graph from the paper I found
x_values <- seq(from = 0.1, to = 10, by =  0.1)

length_x_values <- length(x_values)

# number of points to sample for each value of shape parameter at each x
sample_size <- 5

# x vector
x <- NULL
for (i in 1:length_x_values)
   x <- c(x, rep(x_values[i], sample_size))
(length_x <- length(x))


beta_0 <- 0
beta_1 <- 1

# expectations:  log(mu) = beta_0 + beta_1 * x
(mu_values <- exp(beta_0 + beta_1 * x_values))
(mu <- exp(beta_0 + beta_1 * x))


# shape parameters (dispersion parameters); the number of successes to acheive

#vector of length 10
#n <- c(0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000)

# values of shape parameter that fit the values in the paper I found
#n <- c(0.05, 0.1, 0.5, 1, 5, 10)
n <- c(0.25, 0.5, 1, 5, 10, 20)
(length_n <- length(n))

# probability of success matrix.
# Each row corresponds to an x.
# Each column corresponds to a value of n.
# p  <-  n / (n + mu)
p <- matrix(rep(0, length_x * length_n), nrow = length_x)
for (i in 1:length_x){
   for(j in 1:length_n){
      p[i, j] <- n[j] / (n[j] + mu[i])
   }
}
dim(p)


# run n_simulations simulations
n_simulations <- 100


# AIC matrices.
# Each column corresponds to a value of n.
# Each row corresponds to the number of simulation.
AIC_NB_glm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
AIC_VarStab_lm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
AIC_Poisson_glm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
AIC_log_lm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
AIC_sqrt_lm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
AIC_lm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)

# Mean Squared Error matrices.
# Each column corresponds to a value of n.
# Each row corresponds to the number of simulation.
MSE_NB_glm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
MSE_VarStab_lm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
MSE_Poisson_glm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
MSE_log_lm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
MSE_sqrt_lm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
MSE_lm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)

# Negative Log likelihood (NLL) matrices.
# Each column corresponds to a value of n.
# Each row corresponds to the number of simulation.
NLL_NB_glm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
NLL_VarStab_lm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
NLL_Poisson_glm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
NLL_log_lm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
NLL_sqrt_lm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)
NLL_lm <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)



# estimated theta (n, dispersion parameter), estimated by the NB model, to be used by VarStab 
NB_estimated_theta <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)


# mu estimated by the NB model 
#NB_estimated_mu <- matrix(rep(0, n_simulations * length_n), nrow = n_simulations)



for (num_of_simulation in 1:n_simulations){
  
   # mu estimated by the NB model 
   NB_estimated_mu <- matrix(rep(0, length_x * length_n), nrow = length_n)
  
   
   # for each shape parameter, sample sample_size=5 points at each of the x values.
   # So we have length_n=6 samples of size 100*5=500 each.
   # That is a total of 3000 points.
   # In the samples matrix we have  length_x = sample_size * length_x_values  rows,
   # and length_n columns.  
   NB_sample <- matrix(rep(0, length_n * length_x), nrow = length_x)
   NB_sample_expectations <- matrix(rep(0, length_n * length_x), nrow = length_x)
   for (j in 1:length_n){
      for(x_value in 1:length_x_values){
        NB_sample[(1 + (x_value - 1)*sample_size) : (x_value * sample_size), j] <- rnbinom(sample_size, mu=mu_values[x_value], size=n[j])
        NB_sample_expectations[(1 + (x_value - 1)*sample_size) : (x_value * sample_size), j] <- mu_values[x_value]
      }
   } 
   
   # test sample
   NB_sample_new <- matrix(rep(0, length_n * length_x), nrow = length_x)
   NB_sample_new_expectations <- matrix(rep(0, length_n * length_x), nrow = length_x)
   for (j in 1:length_n){
     for(x_value in 1:length_x_values){
        NB_sample_new[(1 + (x_value - 1)*sample_size) : (x_value * sample_size), j] <- rnbinom(sample_size, mu=mu_values[x_value], size=n[j])
        NB_sample_new_expectations[(1 + (x_value - 1)*sample_size) : (x_value * sample_size), j] <- mu_values[x_value]
     }
   }
  
   for(j in 1:length_n){
       
     
      # NB glm model
      NB_glm_model <- glm.nb(NB_sample[, j] ~ x, control=glm.control(maxit=30), init.theta = n[j])
      
      # Keren's version, which also doesn't work for theta = 0.05:
      #data_frame <- as.data.frame(x = cbind(NB_sample[, j], x))
      #colnames(data_frame) <- c("NB_sample", "x")
      #formula_for_model <- as.formula("NB_sample ~ x")
      #NB_glm_model <- find_model.nb.data(curr_data = data_frame, curr_model_names = formula_for_model)
     
      # keep the estimated theta (n, the dispersion) for use by VarStab
      #NB_estimated_theta[num_of_simulation, j] <- (summary(NB_glm_model))$theta
      NB_estimated_theta[num_of_simulation, j] <- (summary(NB_glm_model))$theta
     
     
      # keep the estimated mu
      NB_estimated_mu[j, ] <- predict(NB_glm_model, type = "response")
      
      # AIC of NB glm model
      AIC_NB_glm[num_of_simulation, j] <- AIC(NB_glm_model)
     
      # mean squared loss for NB glm model on original scale
      # MSE_NB_glm[num_of_simulation, j] <- mean( (predict.glm(NB_glm_model, type = "response") - NB_sample_new[, j])^2 )
      MSE_NB_glm[num_of_simulation, j] <- mean( (predict.glm(NB_glm_model, type = "response") - NB_sample_new_expectations[, j])^2 )
      
      # Negative log likelihood of the test sample for NB model on original scale
      NLL_NB_glm[num_of_simulation, j] <- sum(-dnbinom(x = NB_sample_new[, j],
                                                       size = NB_estimated_theta[num_of_simulation, j],
                                                       mu = NB_estimated_mu[j, ], log = TRUE))
      #SRNLL_NB_glm[num_of_simulation, j] <- sum(-dnbinom(x = NB_sample_new_expectations[, j],
      #                                                 size = NB_estimated_theta[num_of_simulation, j],
      #                                                 mu = NB_estimated_mu[j, ], log = TRUE))

   
      # AIC of known variance stabilizing transformation for NB and linear regression
      #AIC_VarStab_lm[num_of_simulation, j] <- AIC(lm((sqrt(n[j])*asinh(sqrt(NB_sample[, j] / n[j]))) ~ x))
     
      # lm model after variance stabilizing transformation for NB (using estimated theta from NB model) 
      #VarStab_lm_model <- lm((sqrt(NB_estimated_theta[num_of_simulation, j])*asinh(sqrt(NB_sample[, j] / NB_estimated_theta[num_of_simulation, j]))) ~ x)
      
      # lm model after variance stabilizing transformation for NB (using true parameter theta) 
      VarStab_lm_model <- lm((sqrt(n[j])*asinh(sqrt(NB_sample[, j] / n[j]))) ~ x)
      
      # AIC of variance stabilizing transformation for NB (using estimated theta from NB model)
      # and linear regression
      AIC_VarStab_lm[num_of_simulation, j] <- AIC(VarStab_lm_model)
      
      # mean squared loss for NB glm model on original scale
      #MSE_VarStab_lm[num_of_simulation, j] <- mean( (NB_estimated_theta[num_of_simulation, j]*(sinh(predict.lm(VarStab_lm_model, type = "response")/sqrt(NB_estimated_theta[num_of_simulation, j]))^2 - NB_sample_new[, j])^2 ))
      #MSE_VarStab_lm[num_of_simulation, j] <- mean( (n[j]*(sinh(predict.lm(VarStab_lm_model, type = "response")/sqrt(n[j]))^2 - NB_sample_new[, j])^2 ))
      MSE_VarStab_lm[num_of_simulation, j] <- mean( (n[j]*(sinh(predict.lm(VarStab_lm_model, type = "response")/sqrt(n[j]))^2 - NB_sample_new_expectations[, j])^2 ))
      
      
      # Negative log likelihood of the test sample for varstab + lm  model ON VARSTAB SCALE (perhaps should be on original scale)
      #NLL_VarStab_lm[num_of_simulation, j] <- sum(-dnorm(x = (sqrt(NB_estimated_theta[num_of_simulation, j])*asinh(sqrt(NB_sample_new[, j] / NB_estimated_theta[num_of_simulation, j]))),
      #                                                   mean = predict.lm(VarStab_lm_model, type = "response"),
      #                                                   sd = predict.lm(VarStab_lm_model, type = "response", se.fit = TRUE)$se.fit,
      #                                             log = TRUE))
      NLL_VarStab_lm[num_of_simulation, j] <- sum(-dnorm(x = (sqrt(n[j])*asinh(sqrt(NB_sample_new[, j] / n[j]))),
                                                         mean = predict.lm(VarStab_lm_model, type = "response"),
                                                         sd = predict.lm(VarStab_lm_model, type = "response", se.fit = TRUE)$se.fit,
                                                   log = TRUE))
      #SR NLL_VarStab_lm[num_of_simulation, j] <- sum(-dnorm(x = (sqrt(n[j])*asinh(sqrt(NB_sample_new_expectations[, j] / n[j]))),
      #                                                   mean = predict.lm(VarStab_lm_model, type = "response"),
      #                                                   sd = predict.lm(VarStab_lm_model, type = "response", se.fit = TRUE)$se.fit,
      #                                                  log = TRUE))
      
      
      
      # Poisson glm model
      Poisson_glm_model <- glm(NB_sample[, j] ~ x, family="poisson")
      
      # AIC of Poisson glm model
      AIC_Poisson_glm[num_of_simulation, j] <- AIC(Poisson_glm_model)
      
      # mean squared loss for NB glm model on original scale
      # MSE_Poisson_glm[num_of_simulation, j] <- mean( (predict.glm(Poisson_glm_model, type = "response") - NB_sample_new[, j])^2 )
      MSE_Poisson_glm[num_of_simulation, j] <- mean( (predict.glm(Poisson_glm_model, type = "response") - NB_sample_new_expectations[, j])^2 )
      
      # Negative log likelihood of the test sample for Poisson glm  model original scale
      NLL_Poisson_glm[num_of_simulation, j] <- sum(-dpois(x = NB_sample_new[, j],
                                                        lambda = predict.glm(Poisson_glm_model, type = "response", se.fit = TRUE)$se.fit,
                                                    log = TRUE))
     ###SR why se.fit?????                                                
      NLL_Poisson_glm[num_of_simulation, j] <- sum(-dpois(x = NB_sample_new[, j],
                                                        lambda = predict.glm(Poisson_glm_model, type = "response", se.fit = TRUE)$fit,
                                                    log = TRUE))
      
      #SRNLL_Poisson_glm[num_of_simulation, j] <- sum(-dpois(x = NB_sample_new_expectations[, j],
       #                                                   lambda = predict.glm(Poisson_glm_model, type = "response", se.fit = TRUE)$se.fit,
       #                                                   log = TRUE))
      
      
      
      
      # lm model after log (x+1) transformation
      log_lm_model <- lm(log(NB_sample[, j] + 1) ~ x)
      
      # AIC of log transformation and linear regression
      AIC_log_lm[num_of_simulation, j] <- AIC(log_lm_model)
       
      # mean squared loss for NB glm model on original scale
      #MSE_log_lm[num_of_simulation, j] <- mean( (exp(predict.lm(log_lm_model, type = "response")) - (NB_sample_new[, j] + 1))^2 )
      MSE_log_lm[num_of_simulation, j] <- mean( (exp(predict.lm(log_lm_model, type = "response")) - (NB_sample_new_expectations[, j] + 1))^2 )
      
      
      # Negative log likelihood of the test sample for log + lm  model ON LOG SCALE (perhaps should be on original scale)
      NLL_log_lm[num_of_simulation, j] <- sum(-dnorm(x = (log(NB_sample_new[, j] + 1)),
                                                         mean = predict.lm(log_lm_model, type = "response"),
                                                         sd = predict.lm(log_lm_model, type = "response", se.fit = TRUE)$se.fit,
                                                   log = TRUE))
      #SRNLL_log_lm[num_of_simulation, j] <- sum(-dnorm(x = (log(NB_sample_new_expectations[, j] + 1)),
      #                                               mean = predict.lm(log_lm_model, type = "response"),
      #                                               sd = predict.lm(log_lm_model, type = "response", se.fit = TRUE)$se.fit,
      #                                              log = TRUE))
      
      
      
      # lm model after sqrt transformation
      sqrt_lm_model <- lm(sqrt(NB_sample[, j]) ~ x)
      
      # AIC of sqrt transformation and linear regression
      AIC_sqrt_lm[num_of_simulation, j] <- AIC(sqrt_lm_model)
       
      # mean squared loss for NB glm model on original scale
      #MSE_sqrt_lm[num_of_simulation, j] <- mean( ((predict.lm(sqrt_lm_model, type = "response"))^2 - NB_sample_new[, j])^2 )
      MSE_sqrt_lm[num_of_simulation, j] <- mean( ((predict.lm(sqrt_lm_model, type = "response"))^2 - NB_sample_new_expectations[, j])^2 )
      
      # Negative log likelihood of the test sample for sqrt + lm  model ON SQRT SCALE (perhaps should be on original scale)
      NLL_sqrt_lm[num_of_simulation, j] <- sum(-dnorm(x = (sqrt(NB_sample_new[, j])),
                                                         mean = predict.lm(sqrt_lm_model, type = "response"),
                                                         sd = predict.lm(sqrt_lm_model, type = "response", se.fit = TRUE)$se.fit,
                                                log = TRUE))
      
      
      #SRNLL_sqrt_lm[num_of_simulation, j] <- sum(-dnorm(x = (sqrt(NB_sample_new_expectations[, j])),
      #                                                mean = predict.lm(sqrt_lm_model, type = "response"),
      #                                                sd = predict.lm(sqrt_lm_model, type = "response", se.fit = TRUE)$se.fit,
      #                                                log = TRUE))
      
      
      
      
      # lm model
      lm_model <- lm(NB_sample[, j] ~ x)
      
      # AIC of lm
      AIC_lm[num_of_simulation, j] <- AIC(lm_model)
      
      # mean squared loss for NB glm model on original scale
      #MSE_lm[num_of_simulation, j] <- mean( (predict.lm(lm_model, type = "response") - NB_sample_new[, j])^2 )
      MSE_lm[num_of_simulation, j] <- mean( (predict.lm(lm_model, type = "response") - NB_sample_new_expectations[, j])^2 )
      
      # Negative log likelihood of the test sample for lm  model
      NLL_lm[num_of_simulation, j] <- sum(-dnorm(x = NB_sample_new[, j]),
                                                      mean = predict.lm(lm_model, type = "response"),
                                                      sd = predict.lm(lm_model, type = "response", se.fit = TRUE)$se.fit,
                                                log = TRUE)
      
      #SRNLL_lm[num_of_simulation, j] <- sum(-dnorm(x = NB_sample_new_expectations[, j]),
      #                                    mean = predict.lm(lm_model, type = "response"),
      #                                    sd = predict.lm(lm_model, type = "response", se.fit = TRUE)$se.fit,
      #                                    log = TRUE)
      
   }
   
}





## boxplots of AICs for each value of the shape parameter n, with log scale on y axis
## (boxplots AIC for models for NB.pdf)

# global minimum and maximum of AIC
(min_AIC <- min(AIC_NB_glm, AIC_VarStab_lm,  AIC_Poisson_glm, AIC_log_lm, AIC_sqrt_lm, AIC_lm))
(max_AIC <- max(AIC_NB_glm, AIC_VarStab_lm,  AIC_Poisson_glm, AIC_log_lm, AIC_sqrt_lm, AIC_lm))

par(mfcol=c(3, 1))

x_lim_min <- 0
x_lim_max <- 15

x_tic_marks_location <- c(1:6, 9:14)
x_tic_marks_values <- c("NB \nglm", "Var Stab. \nlm", "Pois. \nlm.", "log \nlm", "sqrt \nlm", "lm", "NB \nglm", "Var Stab. \nlm", "Pois. \nlm.", "log \nlm", "sqrt \nlm", "lm")

#y_lim_min <- 0.5
#y_lim_max <- 1+log(80000000)

#y_tic_marks_location <- c(log(2000+1), log(5000+1), log(10000+1), log(20000+1), log(50000+1), log(100000+1), log(200000+1), log(500000+1), log(1000000+1), log(2000000+1), log(5000000+1), log(10000000+1), log(20000000+1), log(40000000+1), log(80000000+1))
#y_tic_marks_values <- c(2000 ,5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 40000000, 80000000)

y_lim_min <- 0.5
y_lim_max <- 1+log(40000000)

y_tic_marks_location <- c(log(100+1), log(1000+1), log(10000+1), log(100000+1), log(1000000+1), log(10000000+1), log(40000000+1))
y_tic_marks_values <- c(100, 1000, 10000, 100000, 1000000, 10000000, 40000000)

plot(x=NA, y=NA, xlab="shape parameter = 0.25                                           shape parameter = 0.5",
     ylab='AIC',
     main="AIC boxplots for the different models",
     cex.main=0.85, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))

axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.7)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.7)

bplot(add=T, x=log(cbind(AIC_NB_glm[, 1], AIC_VarStab_lm[, 1],  AIC_Poisson_glm[, 1], AIC_log_lm[, 1], AIC_sqrt_lm[, 1], AIC_lm[, 1],
                         AIC_NB_glm[, 2], AIC_VarStab_lm[, 2],  AIC_Poisson_glm[, 2], AIC_log_lm[, 2], AIC_sqrt_lm[, 2], AIC_lm[, 2])),
      xlim=c(0,15),
      pos=x_tic_marks_location,
      border=c("black", "blue", "purple", "red", "orange", "green", "black", "blue", "purple", "red", "orange", "green", "pink", "cyan"),
      ann=F, xaxt='n', yaxt='n')

plot(x=NA, y=NA, xlab="shape parameter = 1                                            shape parameter = 5",
     ylab='AIC',
     cex.main=0.85, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))

axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.7)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.7)

bplot(add=T, x=log(cbind(AIC_NB_glm[, 3], AIC_VarStab_lm[, 3],  AIC_Poisson_glm[, 3], AIC_log_lm[, 3], AIC_sqrt_lm[, 3], AIC_lm[, 3],
                         AIC_NB_glm[, 4], AIC_VarStab_lm[, 4],  AIC_Poisson_glm[, 4], AIC_log_lm[, 4], AIC_sqrt_lm[, 4], AIC_lm[, 4])),
      xlim=c(0,15),
      pos=x_tic_marks_location,
      border=c("black", "blue", "purple", "red", "orange", "green", "black", "blue", "purple", "red", "orange", "green", "pink", "cyan"),
      ann=F, xaxt='n', yaxt='n')

plot(x=NA, y=NA, xlab="shape parameter = 10                                    shape parameter = 20",
     ylab='AIC',
     cex.main=0.85, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))

axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.7)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.7)

bplot(add=T, x=log(cbind(AIC_NB_glm[, 5], AIC_VarStab_lm[, 5],  AIC_Poisson_glm[, 5], AIC_log_lm[, 5], AIC_sqrt_lm[, 5], AIC_lm[, 5],
                         AIC_NB_glm[, 6], AIC_VarStab_lm[, 6],  AIC_Poisson_glm[, 6], AIC_log_lm[, 6], AIC_sqrt_lm[, 6], AIC_lm[, 6])),
      xlim=c(0,15),
      pos=x_tic_marks_location,
      border=c("black", "blue", "purple", "red", "orange", "green", "black", "blue", "purple", "red", "orange", "green", "pink", "cyan"),
      ann=F, xaxt='n', yaxt='n')




## boxplots of MSEs for each value of the shape parameter n, with log scale on y axis
## (boxplot MSE for models for NB.pdf)

# Turn outliers into NA's
for(i in 1:n_simulations){
   for (j in 1:length_n){
      if(MSE_VarStab_lm[i, j] > 10^10) MSE_VarStab_lm[i, j] <- NA
      if(NLL_NB_glm[i, j] > 10^6) NLL_NB_glm[i, j] <- NA
   }
}
  
  

# global minimum and maximum of MSE
(min_MSE <- min(MSE_NB_glm, MSE_VarStab_lm,  MSE_Poisson_glm, MSE_log_lm, MSE_sqrt_lm, MSE_lm, na.rm = TRUE))
(max_MSE <- max(MSE_NB_glm, MSE_VarStab_lm,  MSE_Poisson_glm, MSE_log_lm, MSE_sqrt_lm, MSE_lm, na.rm = TRUE))

par(mfcol=c(3, 1))

x_lim_min <- 0
x_lim_max <- 15

x_tic_marks_location <- c(1:6, 9:14)
x_tic_marks_values <- c("NB \nglm", "Var Stab. \nlm", "Pois. \nlm.", "log \nlm", "sqrt \nlm", "lm", "NB \nglm", "Var Stab. \nlm", "Pois. \nlm.", "log \nlm", "sqrt \nlm", "lm")

#y_lim_min <- 0.5
#y_lim_max <- 1+log(500000)

#y_tic_marks_location <- c(log(10+1), log(100+1), log(1000+1), log(10000+1), log(100000+1), log(500000+1))
#y_tic_marks_values <- c(10, 100, 1000, 10000, 100000, 500000)

y_lim_min <- 0.5
y_lim_max <- 1+log(max_MSE)

y_tic_marks_location <- c(log(10^6+1), log(10^7+1), log(10^8+1), log(10^9+1), log(10^10+1), log(10^11+1))
y_tic_marks_values <- c(10^6, 10^7, 10^8, 10^9, 10^10, 10^11)


plot(x=NA, y=NA, xlab="shape parameter = 0.25                                                  shape parameter = 0.5",
     ylab='MSE',
     main="MSE boxplots for the different models",
     cex.main=0.85, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))

axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.7)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.7)

bplot(add=T, x=log(cbind(MSE_NB_glm[, 1], MSE_VarStab_lm[, 1],  MSE_Poisson_glm[, 1], MSE_log_lm[, 1], MSE_sqrt_lm[, 1], MSE_lm[, 1],
                         MSE_NB_glm[, 2], MSE_VarStab_lm[, 2],  MSE_Poisson_glm[, 2], MSE_log_lm[, 2], MSE_sqrt_lm[, 2], MSE_lm[, 2])),
      xlim=c(0,15),
      pos=x_tic_marks_location,
      border=c("black", "blue", "purple", "red", "orange", "green", "black", "blue", "purple", "red", "orange", "green", "pink", "cyan"),
      ann=F, xaxt='n', yaxt='n')

plot(x=NA, y=NA, xlab="shape parameter = 1                                                   shape parameter = 5",
     ylab='MSE',
     cex.main=0.85, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))

axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.7)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.7)

bplot(add=T, x=log(cbind(MSE_NB_glm[, 3], MSE_VarStab_lm[, 3],  MSE_Poisson_glm[, 3], MSE_log_lm[, 3], MSE_sqrt_lm[, 3], MSE_lm[, 3],
                         MSE_NB_glm[, 4], MSE_VarStab_lm[, 4],  MSE_Poisson_glm[, 4], MSE_log_lm[, 4], MSE_sqrt_lm[, 4], MSE_lm[, 4])),
      xlim=c(0,15),
      pos=x_tic_marks_location,
      border=c("black", "blue", "purple", "red", "orange", "green", "black", "blue", "purple", "red", "orange", "green", "pink", "cyan"),
      ann=F, xaxt='n', yaxt='n')

plot(x=NA, y=NA, xlab="shape parameter = 10                                                    shape parameter = 20",
     ylab='MSE',
     cex.main=0.85, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))

axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.7)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.7)

bplot(add=T, x=log(cbind(MSE_NB_glm[, 5], MSE_VarStab_lm[, 5],  MSE_Poisson_glm[, 5], MSE_log_lm[, 5], MSE_sqrt_lm[, 5], MSE_lm[, 5],
                         MSE_NB_glm[, 6], MSE_VarStab_lm[, 6],  MSE_Poisson_glm[, 6], MSE_log_lm[, 6], MSE_sqrt_lm[, 6], MSE_lm[, 6])),
      xlim=c(0,15),
      pos=x_tic_marks_location,
      border=c("black", "blue", "purple", "red", "orange", "green", "black", "blue", "purple", "red", "orange", "green", "pink", "cyan"),
      ann=F, xaxt='n', yaxt='n')




## boxplots of NLLs for each value of the shape parameter n, with log scale on y axis
## (boxplot NLL for models for NB.pdf)

# global minimum and maximum of NLL
(min_NLL <- min(NLL_NB_glm, NLL_VarStab_lm,  NLL_Poisson_glm, NLL_log_lm, NLL_sqrt_lm, NLL_lm, na.rm = TRUE))
(max_NLL <- max(NLL_NB_glm, NLL_VarStab_lm,  NLL_Poisson_glm, NLL_log_lm, NLL_sqrt_lm, NLL_lm, na.rm = TRUE))

par(mfcol=c(3, 1))

x_lim_min <- 0
x_lim_max <- 15

x_tic_marks_location <- c(1:6, 9:14)
x_tic_marks_values <- c("NB \nglm", "Var Stab. \nlm", "Pois. \nlm.", "log \nlm", "sqrt \nlm", "lm", "NB \nglm", "Var. Stab. \nlm", "Pois. \nlm.", "log \nlm", "sqrt \nlm", "lm")

#y_lim_min <- 0.5
#y_lim_max <- 1+log(500000)

#y_tic_marks_location <- c(log(10+1), log(100+1), log(1000+1), log(10000+1), log(100000+1), log(500000+1))
#y_tic_marks_values <- c(10, 100, 1000, 10000, 100000, 500000)

y_lim_min <- 0.5
y_lim_max <- 1+log(max_NLL)

y_tic_marks_location <- c(log(10^3+1), log(10^4+1), log(10^5+1), log(10^6+1), log(10^7+1), log(10^8+1), log(10^9+1))
y_tic_marks_values <- c(10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9)


plot(x=NA, y=NA, xlab="shape parameter = 0.05                                           shape parameter = 0.1",
     ylab='NLL',
     main="NLL boxplots for the different models",
     cex.main=0.85, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))

axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.7)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.7)

bplot(add=T, x=log(cbind(NLL_NB_glm[, 1], NLL_VarStab_lm[, 1],  NLL_Poisson_glm[, 1], NLL_log_lm[, 1], NLL_sqrt_lm[, 1], NLL_lm[, 1],
                         NLL_NB_glm[, 2], NLL_VarStab_lm[, 2],  NLL_Poisson_glm[, 2], NLL_log_lm[, 2], NLL_sqrt_lm[, 2], NLL_lm[, 2])),
      xlim=c(0,15),
      pos=x_tic_marks_location,
      border=c("black", "blue", "purple", "red", "orange", "green", "black", "blue", "purple", "red", "orange", "green", "pink", "cyan"),
      ann=F, xaxt='n', yaxt='n')

plot(x=NA, y=NA, xlab="shape parameter = 0.25                                           shape parameter = 0.5",
     ylab='NLL',
     cex.main=0.85, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))

axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.7)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.7)

bplot(add=T, x=log(cbind(NLL_NB_glm[, 3], NLL_VarStab_lm[, 3],  NLL_Poisson_glm[, 3], NLL_log_lm[, 3], NLL_sqrt_lm[, 3], NLL_lm[, 3],
                         NLL_NB_glm[, 4], NLL_VarStab_lm[, 4],  NLL_Poisson_glm[, 4], NLL_log_lm[, 4], NLL_sqrt_lm[, 4], NLL_lm[, 4])),
      xlim=c(0,15),
      pos=x_tic_marks_location,
      border=c("black", "blue", "purple", "red", "orange", "green", "black", "blue", "purple", "red", "orange", "green", "pink", "cyan"),
      ann=F, xaxt='n', yaxt='n')

plot(x=NA, y=NA, xlab="shape parameter = 1                                    shape parameter = 5",
     ylab='NLL',
     cex.main=0.85, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))

axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.7)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.7)

bplot(add=T, x=log(cbind(NLL_NB_glm[, 5], NLL_VarStab_lm[, 5],  NLL_Poisson_glm[, 5], NLL_log_lm[, 5], NLL_sqrt_lm[, 5], NLL_lm[, 5],
                         NLL_NB_glm[, 6], NLL_VarStab_lm[, 6],  NLL_Poisson_glm[, 6], NLL_log_lm[, 6], NLL_sqrt_lm[, 6], NLL_lm[, 6])),
      xlim=c(0,15),
      pos=x_tic_marks_location,
      border=c("black", "blue", "purple", "red", "orange", "green", "black", "blue", "purple", "red", "orange", "green", "pink", "cyan"),
      ann=F, xaxt='n', yaxt='n')



## Now just for the interesting transformations: Var Stab, log and sqrt, and without log scale on y axis
## boxplots of NLLs for each value of the shape parameter n, with log scale on y axis
## (boxplot NLL for models for NB.pdf)

# global minimum and maximum of NLL
(min_NLL <- min(NLL_VarStab_lm, NLL_log_lm, NLL_sqrt_lm, na.rm = TRUE))
(max_NLL <- max(NLL_VarStab_lm, NLL_log_lm, NLL_sqrt_lm, na.rm = TRUE))

par(mfcol=c(3, 1))

x_lim_min <- 0
x_lim_max <- 9

x_tic_marks_location <- c(1:3, 6:8)
x_tic_marks_values <- c("Var Stab. \nlm", "log \nlm", "sqrt \nlm", "Var. Stab. \nlm", "log \nlm", "sqrt \nlm")


y_lim_min <- min_NLL
y_lim_max <- max_NLL

y_tic_marks_location <- c(10^4, 5*10^4, 10^5, 1.5*10^5, 2*10^5)
y_tic_marks_values <- c(10^4, 5*10^4, 10^5, 1.5*10^5, 2*10^5)


plot(x=NA, y=NA, xlab="shape parameter = 0.25                                           shape parameter = 0.05",
     ylab='NLL',
     main="NLL boxplots for the different models",
     cex.main=0.85, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))

axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.7)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.7)

bplot(add=T, x=(cbind(NLL_VarStab_lm[, 1], NLL_log_lm[, 1], NLL_sqrt_lm[, 1],
                         NLL_VarStab_lm[, 2], NLL_log_lm[, 2], NLL_sqrt_lm[, 2])),
      xlim=c(0, 9),
      pos=x_tic_marks_location,
      border=c("blue", "red", "orange", "blue", "red", "orange"),
      ann=F, xaxt='n', yaxt='n')

plot(x=NA, y=NA, xlab="shape parameter = 1                                           shape parameter = 5",
     ylab='NLL',
     cex.main=0.85, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))

axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.7)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.7)

bplot(add=T, x=(cbind(NLL_VarStab_lm[, 3], NLL_log_lm[, 3], NLL_sqrt_lm[, 3],
                         NLL_VarStab_lm[, 4], NLL_log_lm[, 4], NLL_sqrt_lm[, 4])),
      xlim=c(0, 9),
      pos=x_tic_marks_location,
      border=c("blue", "red", "orange", "blue", "red", "orange"),
      ann=F, xaxt='n', yaxt='n')

plot(x=NA, y=NA, xlab="shape parameter = 10                                    shape parameter = 20",
     ylab='NLL',
     cex.main=0.85, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))

axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.7)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.7)

bplot(add=T, x=(cbind(NLL_VarStab_lm[, 5], NLL_log_lm[, 5], NLL_sqrt_lm[, 5],
                         NLL_VarStab_lm[, 6], NLL_log_lm[, 6], NLL_sqrt_lm[, 6])),
      xlim=c(0, 9),
      pos=x_tic_marks_location,
      border=c("blue", "red", "orange", "blue", "red", "orange"),
      ann=F, xaxt='n', yaxt='n')
