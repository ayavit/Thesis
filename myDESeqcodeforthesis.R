



library("SummarizedExperiment")
library("MASS")
#library ("RCpp")
library ("Rcpp")
library("genefilter")

setwd("Z:/Students/Aya/DESEq2")
source("AllClasses.r")
source("AllGenerics.R")
source("core.r")
source("expanded.r")
source("helper.R")  
source("lfcShrink.R" )
source("Methods.R" )
#source("load.r")
source("parallel.R")
source("plots.R")
source("RcppExports.R")
source("results.R")
source("wrappers.R")
source("rlog.R")
source("vst.R")
source("fitNbinomGLMs.R")


# library for boxplots
library(fields)


# Package for binomial CI
library(binom)

# package for plotting with CIs
library(ggplot2)

# Package for Cochran Armitage test and Somers' d with CI
library(DescTools)

# Package for bootstrap percentile confidence intervals
library(boot)

# Package for calculating Somer's d
library(ryouready)

## Package for computation of odds ratio, CI for odds ratio and Fisher's exact test
library(questionr)

# Package for converting a contingency table to a long format table
library(splitstackshape)

# Package for z.test
library(BSDA)


#num_of_simuls <- 100
num_of_simuls <- 1


source("core.r")






############ Simulation of DESeq2, corrected DESeq2 and our method on data without differential expression, in order to
############ find calibration




diff=FALSE
orig.calibs=corr.calibs=sqrt1.calibs=NULL
orig.p.threshold=corr.p.threshold=sqrt1.p.threshold=NULL
orig.powers=corr.powers=sqrt1.powers=NULL
for (simul in 1:num_of_simuls){
   
   
   ###### corrected DESeq2
   
   
   n_simulations <- 1
   n = 1000
   N <- n_simulations * n
   data_lengths <- rep(n, n_simulations)
   
   
   mus = rep(round(rgamma((2*N), 1/4,0.0003)+1),each=4)
   diff_no_yes_preliminary <- rep(FALSE, 2*N)
   data.NB_all_preliminary = matrix(rnegbin(4*N*2, mu=mus,theta=548),nrow=(N*2),byrow=T)
   data.NB_all <- NULL
   diff_no_yes <- NULL
   iii <- 0
   i <- 1
   while((iii < N) & (i <= N*2)){
      if (sum(data.NB_all_preliminary[i,])>=100){
         iii <- iii + 1
         data.NB_all <- rbind(data.NB_all, data.NB_all_preliminary[i, ])
         diff_no_yes <- c(diff_no_yes, diff_no_yes_preliminary[i])
      }
      i <- i + 1
   }  
   data_lengths <- rep(NA, n_simulations)
   
   
   if (diff){
      # for ~0.5 of the genes, create differential expression: for the treated, multiply mu by 1.5 or 0.5.
      mus = rep(round(rgamma((2*N), 1/4,0.0003)+1),each=4)
      diff_no_yes_preliminary <- as.logical(rbinom(2*N, 1, p = 0.5))
      diff_const_preliminary <- 0.3 * (2 * rbinom(2*N, 1, p = 0.5) - 1) # SR
      for(i in 1:(2*N)){
         if(diff_no_yes_preliminary[i]){
            mus[4 * (i - 1) + 3]  <-  mus[4 * (i - 1) + 3] * (1 + diff_const_preliminary[i])
            mus[4 * (i - 1) + 4]  <-  mus[4 * (i - 1) + 4] * (1 + diff_const_preliminary[i])
         }
      }   
      
      data.NB_all_preliminary = matrix(rnegbin(4*N*2, mu=mus,theta=548),nrow=(N*2),byrow=T)
      data.NB_all <- NULL
      diff_no_yes <- NULL
      iii <- 0
      i <- 1
      while((iii < N) & (i <= N*2)){
         if (sum(data.NB_all_preliminary[i,])>=100){
            iii <- iii + 1
            data.NB_all <- rbind(data.NB_all, data.NB_all_preliminary[i, ])
            diff_no_yes <- c(diff_no_yes, diff_no_yes_preliminary[i])
         }
         i <- i + 1
      }  
      data_lengths <- rep(NA, n_simulations)
   }
   #write.csv(x = data.NB_all, file = "D:/Aya/Thesis/simulated_data_without_diff_exp_1_3_22.csv")
   #write.csv(x = diff_no_yes, file = "D:/Aya/Thesis/diff_no_yes_1_2_22.csv")
   
   
   
   
   for (ver in c("orig","corr")){
      sourceCpp(paste("DESeq2_",ver,".cpp",sep=""))
      
      ######## Corrected DESeq2 on simulated data
      mu_Anders_Huber <- NULL
      alpha_Anders_Huber <- NULL
      p_values_Anders_Huber <- NULL
      p_adj_values_Anders_Huber <- NULL
      num_of_tests_Anders_Huber <- rep(NA, n_simulations)
      k <- 1
      
      data.NB <- data.NB_all[(1+n*(k-1)):(n*k), ]
      passilaCountTable = data.NB
      dim(passilaCountTable)
      row.names(passilaCountTable) <- c(1:dim(data.NB)[1])
      # Get a description of the sample types
      passilaDesign = data.frame(
         row.names = colnames(passilaCountTable),
         condition = c("untreated", "untreated", "treated", "treated"),
         libType = c("paired-end", "paired-end", "paired-end", "paired-end")
      )
      #passilaDesign
      coldata <- passilaDesign
      colnames(coldata) <- c("condition","type")
      coldata$condition <- factor(coldata$condition)
      coldata$type <- factor(coldata$type)
      #coldata
      # Use only paired-end
      coldata <- coldata[coldata$type=="paired-end", ]
      #coldata
      pairedSamples = passilaDesign$libType == "paired-end"
      cts = passilaCountTable[ , pairedSamples]
      condition = passilaDesign$condition[pairedSamples]
      #head(cts)
      # construct a DESeqDataSet
      dds <- DESeqDataSetFromMatrix(countData = data.NB,
                                    colData = coldata,
                                    design = ~ condition)
      #dds
      
      # Pre-filtering: While it is not necessary to pre-filter low count genes before running the DESeq2 functions, there are two reasons which
      # make pre-filtering useful: by removing rows in which there are very few reads, we reduce the memory size of the dds data object, and we
      # increase the speed of the transformation and testing functions within DESeq2. Here we perform a minimal pre-filtering to keep only rows
      # that have at least 100 reads total. More strict filtering to increase power is automatically applied via independent filtering
      # on the mean of normalized counts within the results function.
      keep <- rowSums(counts(dds)) >= 100
      dds <- dds[keep,]
      
      # Define the reference level
      dds$condition <- relevel(dds$condition, ref = "untreated")
      
      # Differential expression analysis
      sizeFactors(dds)=rep(1, dim(passilaCountTable)[2])
      dds <- DESeq(dds)
      res <- results(dds)
      #res
      
      # The estimated mu's
      mu_Anders_Huber <- c(mu_Anders_Huber, mcols(dds,use.names=TRUE)$baseMean)
      #length(mcols(dds,use.names=TRUE)$baseMean)
      
      
      # Anders and Huber estimate the dispersion alpha for each gene, assuming that the variance is  v = s * mu + alpha * s^2 * mu^2, where mu is the expected normalized count data and s is the size factor. Then they fit a curve through the estimate, using a parametric fit: A gamma-family GLM, where two coefficients alpha_0 and alpha_1 are found to parametrize the fit as  alpha = alpha_0 + alpha_1 / mu. Then, for each gene, if the per-gene estimate lies below the regression line, they assume that this might be sampling variance and shift the estimate upwards to the value predicted by the regression line. If, however, the per-gene estimate lies above the line, they keep it as is. So they are prudent in estimating the dispersion as the bigger value.
      
      # The estimated dispersions(alpha's)
      alpha_Anders_Huber <- c(alpha_Anders_Huber, dispersions(dds))
      #length(alpha_Anders_Huber)
      #summary(alpha_Anders_Huber)
      
      # Plot dispersion estimates and fitted values
      plotDispEsts(dds)
      abline(a = 1/548, b = 0, col = "green", untf = TRUE)
      axis(2, at = 1/548, label = "1/548")
      dispersions(dds)
      
      # Final dispersions - the dispersions according to Saharon
      mcols(dds)$dispGeneEst
      
      
      cat (ver, mean(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[!diff_no_yes]<0.05),"\n")
      cat (ver, "calib:", quantile(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[!diff_no_yes],0.05),"\n")
      
      
      if (ver=="orig"){ 
         DESeq2_final_dispersions <- mcols(dds)$dispGeneEst
         orig.calibs = c(orig.calibs, quantile(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[!diff_no_yes],0.05))
         #orig.powers = c(orig.powers, mean(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[diff_no_yes]<0.001))
      }
            
      if (ver=="corr"){
         DESeq2_corrected_final_dispersions <- mcols(dds)$dispGeneEst
         corr.calibs = c(corr.calibs, quantile(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[!diff_no_yes],0.05))
         #corr.powers = c(corr.powers, mean(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[diff_no_yes]<0.001))
      }
      
      
      if (ver=="orig"){
         
         penalized.lik = function(means,  x, theta, cases){
            #  iii = iii+1
            #  cat (means, theta,iii,"\n")
            means = rep(means,each=2)
            w = 1/(1/means + 1/theta)
            return (-sum(log(dnbinom(x,mu=means,size=theta)))+0.5*log(sum(w[cases])*sum(w[-cases])))
         }
         
         infun_vec1 = function(row, scales){
            x = c(0,0,1,1)
            #  mod = glm(row[-5] ~ x,family=negative.binomial(row[5]), offset=log(scales))
            par = c(mean(row[1:2]), mean(row[3:4]))
            res=optim(par, penalized.lik, control = list(reltol=1e-5),x=row[1:4],theta=row[5],cases=3:4)
            res$value  
         }
         
         
         
         model_all_sqrt_tran1 = function(par, scales=NULL)
         {
            
            disp = 1/(par[1] + par[2] * sqrt(pmax(1, apply(data.NB, 1, sum))))
            data32 = cbind(data.NB,disp)
            if (min(disp) < 0.0001){
               Inf
            } else {
               
               # The result of our MLE
               #    res = sum(apply (data32, 1, infun_vec, scales))
               
               # The result of MLE with their penalized likelihood
               res1 = sum(apply (data32, 1, infun_vec1, scales))	
               #    cat(par, res, res1,"\n")
               res1
            }
         }
         
         par1 = c(1/548, 0)
         res1 = optim(par1, model_all_sqrt_tran1, control = list(reltol=1e-5),scales = NULL)
         cat ('moments:', res1$par,'\n')
         
         disp_model_sqrt1 <- 1/(res1$par[1] + res1$par[2] * sqrt(pmax(1, apply(data.NB, 1, sum))))
         final_dispersion_model_sqrt1 <- 1/disp_model_sqrt1
         x = c(0,0,1,1)
         p_values_sqrt1 = NULL
         for (ii in 1: dim(data.NB)[1]){
            model0 <- glm(data.NB[ii, ] ~ 1, family=negative.binomial(disp_model_sqrt1[ii]))    
            model1 <- glm(data.NB[ii, ] ~ x, family=negative.binomial(disp_model_sqrt1[ii]))              
            p_values_sqrt1 <- c(p_values_sqrt1, 1 - pchisq(2 * (logLik(model1) - logLik(model0)), df = 1))
         }
         
         #####
         ## SR
         #####
         
         penalized.loglik = function(theta, x, means, cases){
            w = 1/(1/means + 1/theta)
            return (sum(log(dnbinom(x,mu=means,size=theta)))-0.5*log(sum(w[cases])*sum(w[-cases])))
         }
         
         raw_disps = NULL
         dispersion_model_sqrt1 = numeric ( dim(data.NB)[1])
         for (ii in 1: dim(data.NB)[1]){
            row = data.NB[ii, ]
            mus = c(rep(mean(row[1:2]),2),rep(mean(row[3:4]),2))
            res = optimize(penalized.loglik, interval=c(1,10^6), x=data.NB[ii,], means=mus, cases=c(3,4), maximum=TRUE)
            dispersion_model_sqrt1[ii] = 1/res$maximum
         }
         
         ####
         ####
         
         
         sqrt1.calibs = c(sqrt1.calibs, quantile(p_values_sqrt1[!diff_no_yes],0.05))
         #sqrt1.powers = c(sqrt1.powers, mean(p_values_sqrt1[diff_no_yes]<0.001))
      }   
   }
   #  if (simul %%10 == 0){
   cat (simul,"\n")
   print (summary(orig.calibs))
   print (summary(corr.calibs))
   print (summary(sqrt1.calibs))
   
   #print (summary(orig.powers))
   #print (summary(corr.powers))
   #print (summary(sqrt1.powers))
   
   #   }
}






diff=TRUE
orig.calibs=corr.calibs=sqrt1.calibs=NULL
orig.powers=corr.powers=sqrt1.powers=NULL
orig.powers.calib=corr.powers.calib=NULL
for (simul in 1:num_of_simuls){

   ############ Simulation of corrected DESeq2 on data with differential expression


   n_simulations <- 1
   n = 1000
   N <- n_simulations * n
   data_lengths <- rep(n, n_simulations)


   mus = rep(round(rgamma((2*N), 1/4,0.0003)+1),each=4)
   diff_no_yes_preliminary <- rep(FALSE, 2*N)
   data.NB_all_preliminary = matrix(rnegbin(4*N*2, mu=mus,theta=548),nrow=(N*2),byrow=T)
   data.NB_all <- NULL
   diff_no_yes <- NULL
   iii <- 0
   i <- 1
   while((iii < N) & (i <= N*2)){
      if (sum(data.NB_all_preliminary[i,])>=100){
         iii <- iii + 1
         data.NB_all <- rbind(data.NB_all, data.NB_all_preliminary[i, ])
         diff_no_yes <- c(diff_no_yes, diff_no_yes_preliminary[i])
      }
      i <- i + 1
   }  
   data_lengths <- rep(NA, n_simulations)


   if (diff){
      # for ~0.5 of the genes, create differential expression: for the treated, multiply mu by 1.5 or 0.5.
      mus = rep(round(rgamma((2*N), 1/4,0.0003)+1),each=4)
      diff_no_yes_preliminary <- as.logical(rbinom(2*N, 1, p = 0.5))
      diff_const_preliminary <- 0.3 * (2 * rbinom(2*N, 1, p = 0.5) - 1) # SR
      for(i in 1:(2*N)){
         if(diff_no_yes_preliminary[i]){
            mus[4 * (i - 1) + 3]  <-  mus[4 * (i - 1) + 3] * (1 + diff_const_preliminary[i])
            mus[4 * (i - 1) + 4]  <-  mus[4 * (i - 1) + 4] * (1 + diff_const_preliminary[i])
         }
      }   

      data.NB_all_preliminary = matrix(rnegbin(4*N*2, mu=mus,theta=548),nrow=(N*2),byrow=T)
      data.NB_all <- NULL
      diff_no_yes <- NULL
      iii <- 0
      i <- 1
      while((iii < N) & (i <= N*2)){
         if (sum(data.NB_all_preliminary[i,])>=100){
            iii <- iii + 1
            data.NB_all <- rbind(data.NB_all, data.NB_all_preliminary[i, ])
            diff_no_yes <- c(diff_no_yes, diff_no_yes_preliminary[i])
         }
         i <- i + 1
      }  
      data_lengths <- rep(NA, n_simulations)
   }
   #write.csv(x = data.NB_all, file = "D:/Aya/Thesis/simulated_data_without_diff_exp_1_3_22.csv")
   #write.csv(x = diff_no_yes, file = "D:/Aya/Thesis/diff_no_yes_1_2_22.csv")




   for (ver in c("orig","corr")){
      sourceCpp(paste("DESeq2_",ver,".cpp",sep=""))

      ######## Corrected DESeq2 on simulated data
      mu_Anders_Huber <- NULL
      alpha_Anders_Huber <- NULL
      p_values_Anders_Huber <- NULL
      p_adj_values_Anders_Huber <- NULL
      num_of_tests_Anders_Huber <- rep(NA, n_simulations)
      k <- 1

      data.NB <- data.NB_all[(1+n*(k-1)):(n*k), ]
      passilaCountTable = data.NB
      dim(passilaCountTable)
      row.names(passilaCountTable) <- c(1:dim(data.NB)[1])
      # Get a description of the sample types
      passilaDesign = data.frame(
         row.names = colnames(passilaCountTable),
         condition = c("untreated", "untreated", "treated", "treated"),
         libType = c("paired-end", "paired-end", "paired-end", "paired-end")
      )
      #passilaDesign
      coldata <- passilaDesign
      colnames(coldata) <- c("condition","type")
      coldata$condition <- factor(coldata$condition)
      coldata$type <- factor(coldata$type)
      #coldata
      # Use only paired-end
      coldata <- coldata[coldata$type=="paired-end", ]
      #coldata
      pairedSamples = passilaDesign$libType == "paired-end"
      cts = passilaCountTable[ , pairedSamples]
      condition = passilaDesign$condition[pairedSamples]
      #head(cts)
      # construct a DESeqDataSet
      dds <- DESeqDataSetFromMatrix(countData = data.NB,
                              colData = coldata,
                              design = ~ condition)
      #dds

      # Pre-filtering: While it is not necessary to pre-filter low count genes before running the DESeq2 functions, there are two reasons which
      # make pre-filtering useful: by removing rows in which there are very few reads, we reduce the memory size of the dds data object, and we
      # increase the speed of the transformation and testing functions within DESeq2. Here we perform a minimal pre-filtering to keep only rows
      # that have at least 100 reads total. More strict filtering to increase power is automatically applied via independent filtering
      # on the mean of normalized counts within the results function.
      keep <- rowSums(counts(dds)) >= 100
      dds <- dds[keep,]

      # Define the reference level
      dds$condition <- relevel(dds$condition, ref = "untreated")

      # Differential expression analysis
      sizeFactors(dds)=rep(1, dim(passilaCountTable)[2])
      dds <- DESeq(dds)
      res <- results(dds)
      #res

      # The estimated mu's
      mu_Anders_Huber <- c(mu_Anders_Huber, mcols(dds,use.names=TRUE)$baseMean)
      #length(mcols(dds,use.names=TRUE)$baseMean)


      # Anders and Huber estimate the dispersion alpha for each gene, assuming that the variance is  v = s * mu + alpha * s^2 * mu^2, where mu is the expected normalized count data and s is the size factor. Then they fit a curve through the estimate, using a parametric fit: A gamma-family GLM, where two coefficients alpha_0 and alpha_1 are found to parametrize the fit as  alpha = alpha_0 + alpha_1 / mu. Then, for each gene, if the per-gene estimate lies below the regression line, they assume that this might be sampling variance and shift the estimate upwards to the value predicted by the regression line. If, however, the per-gene estimate lies above the line, they keep it as is. So they are prudent in estimating the dispersion as the bigger value.

      # The estimated dispersions(alpha's)
      alpha_Anders_Huber <- c(alpha_Anders_Huber, dispersions(dds))
      #length(alpha_Anders_Huber)
      #summary(alpha_Anders_Huber)

      # Plot dispersion estimates and fitted values
      plotDispEsts(dds)
      abline(a = 1/548, b = 0, col = "green", untf = TRUE)
      axis(2, at = 1/548, label = "1/548")
      dispersions(dds)
      
      # Final dispersions - the dispersions according to Saharon
      mcols(dds)$dispGeneEst
      
      
      cat (ver, mean(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[!diff_no_yes]<0.05),"\n")
      cat (ver, "calib:", quantile(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[!diff_no_yes],0.05),"\n")

      
      if (ver=="orig"){ 
         DESeq2_final_dispersions <- mcols(dds)$dispGeneEst
         orig.calibs = c(orig.calibs, quantile(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[!diff_no_yes],0.05))
         orig.powers = c(orig.powers, mean(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[diff_no_yes]<0.001))
         orig.powers.calib = c(orig.powers, mean(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[diff_no_yes]<0.000832))
      }
      if (ver=="corr"){
         DESeq2_corrected_final_dispersions <- mcols(dds)$dispGeneEst
         corr.calibs = c(corr.calibs, quantile(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[!diff_no_yes],0.05))
         corr.powers = c(corr.powers, mean(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[diff_no_yes]<0.001))
         corr.powers.calib = c(corr.powers, mean(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[diff_no_yes]<0.000738))
      }


      if (ver=="orig"){

         penalized.lik = function(means,  x, theta, cases){
            #  iii = iii+1
            #  cat (means, theta,iii,"\n")
            means = rep(means,each=2)
            w = 1/(1/means + 1/theta)
            return (-sum(log(dnbinom(x,mu=means,size=theta)))+0.5*log(sum(w[cases])*sum(w[-cases])))
         }
    
          infun_vec1 = function(row, scales){
             x = c(0,0,1,1)
             #  mod = glm(row[-5] ~ x,family=negative.binomial(row[5]), offset=log(scales))
             par = c(mean(row[1:2]), mean(row[3:4]))
             res=optim(par, penalized.lik, control = list(reltol=1e-5),x=row[1:4],theta=row[5],cases=3:4)
             res$value  
          }
    
    
    
          model_all_sqrt_tran1 = function(par, scales=NULL)
          {
      
             disp = 1/(par[1] + par[2] * sqrt(pmax(1, apply(data.NB, 1, sum))))
             data32 = cbind(data.NB,disp)
             if (min(disp) < 0.0001){
               Inf
             } else {
        
                # The result of our MLE
                #    res = sum(apply (data32, 1, infun_vec, scales))
        
                # The result of MLE with their penalized likelihood
                res1 = sum(apply (data32, 1, infun_vec1, scales))	
                #    cat(par, res, res1,"\n")
                res1
             }
          }
    
          par1 = c(1/548, 0)
          res1 = optim(par1, model_all_sqrt_tran1, control = list(reltol=1e-5),scales = NULL)
          cat ('moments:', res1$par,'\n')
    
          disp_model_sqrt1 <- 1/(res1$par[1] + res1$par[2] * sqrt(pmax(1, apply(data.NB, 1, sum))))
          final_dispersion_model_sqrt1 <- 1/disp_model_sqrt1
          x = c(0,0,1,1)
          p_values_sqrt1 = NULL
          for (ii in 1: dim(data.NB)[1]){
             model0 <- glm(data.NB[ii, ] ~ 1, family=negative.binomial(disp_model_sqrt1[ii]))    
             model1 <- glm(data.NB[ii, ] ~ x, family=negative.binomial(disp_model_sqrt1[ii]))              
             p_values_sqrt1 <- c(p_values_sqrt1, 1 - pchisq(2 * (logLik(model1) - logLik(model0)), df = 1))
          }
    
          #####
          ## SR
          #####
   
          penalized.loglik = function(theta, x, means, cases){
             w = 1/(1/means + 1/theta)
             return (sum(log(dnbinom(x,mu=means,size=theta)))-0.5*log(sum(w[cases])*sum(w[-cases])))
          }

          raw_disps = NULL
          dispersion_model_sqrt1 = numeric ( dim(data.NB)[1])
          for (ii in 1: dim(data.NB)[1]){
             row = data.NB[ii, ]
             mus = c(rep(mean(row[1:2]),2),rep(mean(row[3:4]),2))
             res = optimize(penalized.loglik, interval=c(1,10^6), x=data.NB[ii,], means=mus, cases=c(3,4), maximum=TRUE)
             dispersion_model_sqrt1[ii] = 1/res$maximum
          }
          
          ####
          ####

        
          sqrt1.calibs = c(sqrt1.calibs, quantile(p_values_sqrt1[!diff_no_yes],0.05))
          sqrt1.powers = c(sqrt1.powers, mean(p_values_sqrt1[diff_no_yes]<0.001))

       }   
   }
#  if (simul %%10 == 0){
      cat (simul,"\n")
      print (summary(orig.calibs))
      print (summary(corr.calibs))
      print (summary(sqrt1.calibs))

      print (summary(orig.powers))
      print (summary(corr.powers))
      print (summary(sqrt1.powers))
      
      print (summary(orig.powers.calib))
      print (summary(corr.powers.calib))
      

#   }
}


plot(DESeq2_final_dispersions ~ dispersion_model_sqrt1, ylab = "DESeq2", xlab = "Improved method")
plot(DESeq2_corrected_final_dispersions ~ dispersion_model_sqrt1, ylab = "corrected DESeq2", xlab = "Improved method")


#plot(final_dispersion_model_sqrt1 ~ rowMeans(data.NB))
plot(final_dispersion_model_sqrt1 ~ rowMeans(data.NB), xlim = c(20, 20000), ylim = c(10^(-6), 1), log = "xy",
     xlab = "Mean counts", ylab = "Improved method")
abline(a = 1/548, b = 0, col = "green", untf = TRUE)
axis(2, at = 1/548, label = "1/548")






diff=TRUE
orig.calibs=corr.calibs=sqrt1.calibs=NULL
orig.discoveries=corr.discoveries=sqrt1.discoveries=NULL
orig.discoveries.calib=corr.discoveries.calib=NULL

num_of_simuls <- 1

for (simul in 1:num_of_simuls){
   
   ############ Real data with differential expression
   
   
   #### Read the data
   #data1 <- read.csv(file = "z:/students/Aya/DESEq2/real_data_original_4_3_20.csv")
   data1 <- read.csv(file = "D:/Aya/Thesis/real_data_original_4_3_20.csv")
   head(data1)
   data1 <- data1[1:14590, ]
   dim(data1)
   
   passilaCountTable = data1[ , -1]
   passilaCountTable[1:100, ]
   dim(passilaCountTable)
   row.names(passilaCountTable) <- data1[, 1]
   passilaCountTable[1:100, ]
   
   # Get a description of the sample types
   passilaDesign = data.frame(
      row.names = colnames(passilaCountTable),
      condition = c("untreated", "untreated", "untreated", "untreated", "treated", "treated", "treated"),
      libType = c("single-end", "single-end", "paired-end", "paired-end", "single-end", "paired-end", "paired-end")
   )
   passilaDesign
   
   
   coldata <- passilaDesign
   colnames(coldata) <- c("condition","type")
   coldata$condition <- factor(coldata$condition)
   coldata$type <- factor(coldata$type)
   coldata
   
   
   # Use only paired-end
   coldata <- coldata[coldata$type=="paired-end", ]
   coldata
   pairedSamples = passilaDesign$libType == "paired-end"
   cts = passilaCountTable[ , pairedSamples]
   condition = passilaDesign$condition[pairedSamples]
   head(cts)
   
   
   
   # construct a DESeqDataSet
   dds <- DESeqDataSetFromMatrix(countData = cts,
                                 colData = coldata,
                                 design = ~ condition)
   dds
   
   
   # Pre-filtering: While it is not necessary to pre-filter low count genes before running the DESeq2 functions, there are two reasons which
   # make pre-filtering useful: by removing rows in which there are very few reads, we reduce the memory size of the dds data object, and we
   # increase the speed of the transformation and testing functions within DESeq2. Here we perform a minimal pre-filtering to keep only rows
   # that have at least 100 reads total. More strict filtering to increase power is automatically applied via independent filtering
   # on the mean of normalized counts within the results function.
   keep <- rowSums(counts(dds)) >= 100
   dds <- dds[keep,]
   
   # Do the same for the matrix that we are going to use
   row_at_least_100 <- rep(FALSE, dim(cts)[1])
   for(i in 1:dim(cts)[1])
      if (sum(cts[i, ]) >= 100)
         row_at_least_100[i] <- TRUE
   data2 <- cts[row_at_least_100, ]
   
   dim(data2)
   # 7141 4
   
   
   
   
   n_simulations <- 1
   #n = 1000
   n <- dim(data2)[1]
   N <- n_simulations * n
   data_lengths <- rep(n, n_simulations)
   
   data.NB_all <- data2
   
   
   

   for (ver in c("orig","corr")){
      sourceCpp(paste("DESeq2_",ver,".cpp",sep=""))
      
      ######## Corrected DESeq2 on simulated data
      mu_Anders_Huber <- NULL
      alpha_Anders_Huber <- NULL
      p_values_Anders_Huber <- NULL
      p_adj_values_Anders_Huber <- NULL
      num_of_tests_Anders_Huber <- rep(NA, n_simulations)
      k <- 1
      
      data.NB <- data.NB_all[(1+n*(k-1)):(n*k), ]
      passilaCountTable = data.NB
      dim(passilaCountTable)
      row.names(passilaCountTable) <- c(1:dim(data.NB)[1])
      # Get a description of the sample types
      passilaDesign = data.frame(
         row.names = colnames(passilaCountTable),
         condition = c("untreated", "untreated", "treated", "treated"),
         libType = c("paired-end", "paired-end", "paired-end", "paired-end")
      )
      #passilaDesign
      coldata <- passilaDesign
      colnames(coldata) <- c("condition","type")
      coldata$condition <- factor(coldata$condition)
      coldata$type <- factor(coldata$type)
      #coldata
      # Use only paired-end
      coldata <- coldata[coldata$type=="paired-end", ]
      #coldata
      pairedSamples = passilaDesign$libType == "paired-end"
      cts = passilaCountTable[ , pairedSamples]
      condition = passilaDesign$condition[pairedSamples]
      #head(cts)
      # construct a DESeqDataSet
      dds <- DESeqDataSetFromMatrix(countData = data.NB,
                                    colData = coldata,
                                    design = ~ condition)
      #dds
      
      # Pre-filtering: While it is not necessary to pre-filter low count genes before running the DESeq2 functions, there are two reasons which
      # make pre-filtering useful: by removing rows in which there are very few reads, we reduce the memory size of the dds data object, and we
      # increase the speed of the transformation and testing functions within DESeq2. Here we perform a minimal pre-filtering to keep only rows
      # that have at least 100 reads total. More strict filtering to increase power is automatically applied via independent filtering
      # on the mean of normalized counts within the results function.
      keep <- rowSums(counts(dds)) >= 100
      dds <- dds[keep,]
      
      # Define the reference level
      dds$condition <- relevel(dds$condition, ref = "untreated")
      
      # Differential expression analysis
      sizeFactors(dds)=rep(1, dim(passilaCountTable)[2])
      dds <- DESeq(dds)
      res <- results(dds)
      #res
      
      # The estimated mu's
      mu_Anders_Huber <- c(mu_Anders_Huber, mcols(dds,use.names=TRUE)$baseMean)
      #length(mcols(dds,use.names=TRUE)$baseMean)
      
      
      # Anders and Huber estimate the dispersion alpha for each gene, assuming that the variance is  v = s * mu + alpha * s^2 * mu^2, where mu is the expected normalized count data and s is the size factor. Then they fit a curve through the estimate, using a parametric fit: A gamma-family GLM, where two coefficients alpha_0 and alpha_1 are found to parametrize the fit as  alpha = alpha_0 + alpha_1 / mu. Then, for each gene, if the per-gene estimate lies below the regression line, they assume that this might be sampling variance and shift the estimate upwards to the value predicted by the regression line. If, however, the per-gene estimate lies above the line, they keep it as is. So they are prudent in estimating the dispersion as the bigger value.
      
      # The estimated dispersions(alpha's)
      alpha_Anders_Huber <- c(alpha_Anders_Huber, dispersions(dds))
      #length(alpha_Anders_Huber)
      #summary(alpha_Anders_Huber)
      
      # Plot dispersion estimates and fitted values
      plotDispEsts(dds)
      #abline(a = 1/548, b = 0, col = "green", untf = TRUE)
      #axis(2, at = 1/548, label = "1/548")
      dispersions(dds)
      
      # Final dispersions - the dispersions according to Saharon
      mcols(dds)$dispGeneEst
      
      
      #cat (ver, mean(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[!diff_no_yes]<0.05),"\n")
      #cat (ver, "calib:", quantile(mcols(dds)$WaldPvalue_condition_treated_vs_untreated[!diff_no_yes],0.05),"\n")
      
      
      if (ver=="orig"){ 
         DESeq2_final_dispersions <- mcols(dds)$dispGeneEst
         orig.discoveries = c(orig.discoveries, mean(mcols(dds)$WaldPvalue_condition_treated_vs_untreated<(0.05/N)))
         orig.discoveries.calib = c(orig.discoveries, mean(mcols(dds)$WaldPvalue_condition_treated_vs_untreated<(0.0416/N)))
      }
      if (ver=="corr"){
         DESeq2_corrected_final_dispersions <- mcols(dds)$dispGeneEst
         corr.discoveries = c(corr.discoveries, mean(mcols(dds)$WaldPvalue_condition_treated_vs_untreated<(0.05/N)))
         corr.discoveries.calib = c(corr.discoveries, mean(mcols(dds)$WaldPvalue_condition_treated_vs_untreated<(0.0369/N)))
      }
      
      
      if (ver=="orig"){
         
         penalized.lik = function(means,  x, theta, cases){
            #  iii = iii+1
            #  cat (means, theta,iii,"\n")
            means = rep(means,each=2)
            w = 1/(1/means + 1/theta)
            return (-sum(log(dnbinom(x,mu=means,size=theta)))+0.5*log(sum(w[cases])*sum(w[-cases])))
         }
         
         infun_vec1 = function(row, scales){
            x = c(0,0,1,1)
            #  mod = glm(row[-5] ~ x,family=negative.binomial(row[5]), offset=log(scales))
            par = c(mean(row[1:2]), mean(row[3:4]))
            res=optim(par, penalized.lik, control = list(reltol=1e-5),x=row[1:4],theta=row[5],cases=3:4)
            res$value  
         }
         
         
         
         model_all_sqrt_tran1 = function(par, scales=NULL)
         {
            
            disp = 1/(par[1] + par[2] * sqrt(pmax(1, apply(data.NB, 1, sum))))
            data32 = cbind(data.NB,disp)
            if (min(disp) < 0.0001){
               Inf
            } else {
               
               # The result of our MLE
               #    res = sum(apply (data32, 1, infun_vec, scales))
               
               # The result of MLE with their penalized likelihood
               res1 = sum(apply (data32, 1, infun_vec1, scales))	
               #    cat(par, res, res1,"\n")
               res1
            }
         }
         
         par1 = c(1/548, 0)
         res1 = optim(par1, model_all_sqrt_tran1, control = list(reltol=1e-5),scales = NULL)
         cat ('moments:', res1$par,'\n')
         
         disp_model_sqrt1 <- 1/(res1$par[1] + res1$par[2] * sqrt(pmax(1, apply(data.NB, 1, sum))))
         final_dispersion_model_sqrt1 <- 1/disp_model_sqrt1
         x = c(0,0,1,1)
         p_values_sqrt1 = NULL
         for (ii in 1: dim(data.NB)[1]){
            row = as.numeric(data.NB[ii,])
            model0 <- glm(row ~ 1, family=negative.binomial(disp_model_sqrt1[ii]))    
            model1 <- glm(row ~ x, family=negative.binomial(disp_model_sqrt1[ii]))              
            p_values_sqrt1 <- c(p_values_sqrt1, 1 - pchisq(2 * (logLik(model1) - logLik(model0)), df = 1))
         }
         
         #####
         ## SR
         #####
         
         penalized.loglik = function(theta, x, means, cases){
            w = 1/(1/means + 1/theta)
            return (sum(log(dnbinom(x,mu=means,size=theta)))-0.5*log(sum(w[cases])*sum(w[-cases])))
         }
         
         raw_disps = NULL
         dispersion_model_sqrt1 = numeric ( dim(data.NB)[1])
         for (ii in 1: dim(data.NB)[1]){
            row = as.numeric(data.NB[ii, ])
            mus = c(rep(mean(row[1:2]),2),rep(mean(row[3:4]),2))
            res = optimize(penalized.loglik, interval=c(1,10^6), x=row, means=mus, cases=c(3,4), maximum=TRUE)
            dispersion_model_sqrt1[ii] = 1/res$maximum
         }
         
         ####
         ####
         
         
         #sqrt1.calibs = c(sqrt1.calibs, quantile(p_values_sqrt1[!diff_no_yes],0.05))
         sqrt1.discoveries = c(sqrt1.discoveries, mean(p_values_sqrt1<(0.05/N)))
         
      }   
   }
   #  if (simul %%10 == 0){
   cat (simul,"\n")
   #print (summary(orig.calibs))
   #print (summary(corr.calibs))
   #print (summary(sqrt1.calibs))
   
   print (summary(orig.discoveries))
   print (summary(corr.discoveries))
   print (summary(sqrt1.discoveries))
   
   print (summary(orig.discoveries.calib))
   print (summary(corr.discoveries.calib))
   
   
   #   }
}
