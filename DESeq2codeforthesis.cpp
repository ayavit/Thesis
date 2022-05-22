
//********************************************************************************
//*** this is an extracted function from the file DESeq2.cpp, it is not standalone
//*** code!
//*** our changes are highlighted
//********************************************************************************
// Obtain the MLE or MAP dispersion estimate using line search.
// fitting occurs on the scale of log(alpha)
//
// [[Rcpp::export]]
List fitDisp(SEXP ySEXP, SEXP xSEXP, SEXP mu_hatSEXP, SEXP log_alphaSEXP,
 SEXP log_alpha_prior_meanSEXP, SEXP log_alpha_prior_sigmasqSEXP,
 SEXP min_log_alphaSEXP, SEXP kappa_0SEXP, SEXP tolSEXP, SEXP maxitSEXP,
 SEXP usePriorSEXP, SEXP weightsSEXP, SEXP useWeightsSEXP,
 SEXP weightThresholdSEXP, SEXP useCRSEXP) {
  NumericMatrix y(ySEXP);
  arma::mat x = as<arma::mat>(xSEXP);
  int y_n = y.nrow();
  NumericVector log_alpha(clone(log_alphaSEXP));
  NumericMatrix mu_hat(mu_hatSEXP);
  NumericVector log_alpha_prior_mean(log_alpha_prior_meanSEXP);
  double log_alpha_prior_sigmasq = as<double>(log_alpha_prior_sigmasqSEXP);
  double min_log_alpha = as<double>(min_log_alphaSEXP);
  double kappa_0 = as<double>(kappa_0SEXP);
  int maxit = as<int>(maxitSEXP);
  double epsilon = 1.0e-4;
  double a, a_propose, kappa, lp, lpnew, dlp, theta_kappa, theta_hat_kappa,
 change;
  // record log posterior values
  NumericVector initial_lp(y_n);
  NumericVector initial_dlp(y_n);
  NumericVector last_lp(y_n);
  NumericVector last_dlp(y_n);
  NumericVector last_d2lp(y_n);
  NumericVector last_change(y_n);
  IntegerVector iter(y_n);
  IntegerVector iter_accept(y_n);
  double tol = as<double>(tolSEXP);
  bool usePrior = as<bool>(usePriorSEXP);
  // observation weights
  NumericMatrix weights(weightsSEXP);
  bool useWeights = as<bool>(useWeightsSEXP);
  double weightThreshold = as<double>(weightThresholdSEXP);
  bool useCR = as<bool>(useCRSEXP);


// outer loop
  for (int i = 0; i < y_n; i++) {
    if (i % 100 == 0) checkUserInterrupt();
    NumericMatrix::Row yrow = y(i,_);
    NumericMatrix::Row mu_hat_row = mu_hat(i,_);
    // maximize the log likelihood over the variable a, the log of alpha, the
 dispersion parameter.
    // in order to express the optimization in a typical manner, 
    // for calculating theta(kappa) we multiple the log likelihood by -1 and seek
 a minimum
    a = log_alpha(i);
    // we use a line search based on the Armijo rule.
    // define a function theta(kappa) = f(a + kappa * d), where d is the search
 direction.
    // in this case the search direction is taken by the first derivative of the
 log likelihood
    
    
    lp = log_posterior(a, yrow, mu_hat_row, x, log_alpha_prior_mean(i),
 log_alpha_prior_sigmasq, usePrior, weights.row(i), useWeights, weightThreshold,
 useCR);
    dlp = dlog_posterior(a, yrow, mu_hat_row, x, log_alpha_prior_mean(i),
 log_alpha_prior_sigmasq, usePrior, weights.row(i), useWeights, weightThreshold,
 useCR);
    kappa = kappa_0;
    initial_lp(i) = lp;
    initial_dlp(i) = dlp;
    change = -1.0;
    last_change(i) = -1.0;
    

//#### 
//## The following is our new heuristic
//####
    bool NewHeuristic=true; // Make it false to get DESeq2 standard optimization
    if( NewHeuristic){
        double besta = a, bestlp=lp;
        
        double center = -10;
        for (double scale =1; scale >= 0.001; scale/=10){
            for (int ii =0;ii<=20;ii++){
                double anow = center+(ii-10)*scale;
                double lpnow = log_posterior(anow, yrow, mu_hat_row, x,
 log_alpha_prior_mean(i), log_alpha_prior_sigmasq, usePrior, weights.row(i),
 useWeights, weightThreshold, useCR);
    
                if (lpnow >bestlp){
                    bestlp = lpnow;
                    besta=anow;
                }
            }
            center=besta;
        }
        log_alpha(i) = besta;
// dummy values
        last_lp(i) = bestlp+1;
        last change indicates the change for the final iteration
        last_change(i) = 0;
        iter(i)=2;        
// This line makes sure we skip over the DESeq2 standard solution, which
 is the rest of the outer loop
        continue;
    }
//#### 
//## end new heuristic
//####

    
//#### 
//Back to DESeq2's standard code
//####
	
    for (int t = 0; t < maxit; t++) {

      // iter counts the number of steps taken out of  maxit;
      iter(i)++;
//      if (fabs(kappa*dlp)<0.1){ //SRSR
  //     kappa = 0.1/fabs(dlp);
    //   if (i<10) Rcout<<"kappa changed to: "<<kappa<<"\n";
      //}  //SRSR
      a_propose = a + kappa * dlp;
//      if (i<10) Rcout<<i<<" "<<iter(i)<<" "<<a<<" "<<a_propose-a<<"
 "<<kappa<<"\n"; //SR
      // note: lgamma is unstable for values around 1e17, where there
 is a switch in lgamma.c
      // we limit log alpha from going lower than -30
      if (a_propose < -30.0) {
	      kappa = (-30.0 - a)/dlp;
      }
      // note: we limit log alpha from going higher than 10
      if (a_propose > 10.0) {
	      kappa = (10.0 - a)/dlp;
      }
      theta_kappa = -1.0 * log_posterior(a + kappa*dlp, yrow, mu_hat_row,
 x, log_alpha_prior_mean(i), log_alpha_prior_sigmasq, usePrior,
 weights.row(i), useWeights, weightThreshold, useCR);
      theta_hat_kappa = -1.0 * lp - kappa * epsilon * R_pow_di(dlp, 2);
//      if (i<10) Rcout<<kappa<<" "<<theta_kappa<<"
 "<<theta_hat_kappa-theta_kappa<<"\n"; //SR
      // if this inequality is true, we have satisfied the Armijo rule and 
      // accept the step size kappa, otherwise we halve kappa
      if (theta_kappa <= theta_hat_kappa) {
      	// iter_accept counts the number of accepted proposals;
      	iter_accept(i)++;
      	a = a + kappa * dlp;
      	lpnew = log_posterior(a, yrow, mu_hat_row, x, log_alpha_prior_mean(i),
 log_alpha_prior_sigmasq, usePrior, weights.row(i), useWeights,
 weightThreshold, useCR);
      	// look for change in log likelihood
      	change = lpnew - lp;
//        if (i<10) Rcout<<change<<" "<<tol<<"\n"; //SR

      	if (change < tol) {
      	  lp = lpnew;
      	  break;
      	}
      	// if log(alpha) is going to -infinity
      	// break the loop
      	if (a < min_log_alpha) {
      	  break;
      	}
      	lp = lpnew;
      	dlp = dlog_posterior(a, yrow, mu_hat_row, x, log_alpha_prior_mean(i),
 log_alpha_prior_sigmasq, usePrior, weights.row(i), useWeights,
 weightThreshold, useCR);
      	// instead of resetting kappa to kappa_0 
      	// multiple kappa by 1.1
      	kappa = fmin(kappa * 1.1, kappa_0);
      	// every 5 accepts, halve kappa
      	// to prevent slow convergence
      	// due to overshooting
      	if (iter_accept(i) % 5 == 0) {
      	  kappa = kappa / 2.0;
      	}
      } else {
	kappa = kappa / 2.0;
      }
    }
    last_lp(i) = lp;
    last_dlp(i) = dlp;
    last_d2lp(i) = d2log_posterior(a, yrow, mu_hat_row, x,
 log_alpha_prior_mean(i), log_alpha_prior_sigmasq, usePrior, weights.row(i),
 useWeights, weightThreshold, useCR);
    log_alpha(i) = a;
    // last change indicates the change for the final iteration
    last_change(i) = change;
  }

  return List::create(Named("log_alpha",log_alpha),
		      Named("iter",iter),
		      Named("iter_accept",iter_accept),
		      Named("last_change",last_change),
		      Named("initial_lp",initial_lp),
		      Named("initial_dlp",initial_dlp),
		      Named("last_lp",last_lp),
		      Named("last_dlp",last_dlp),
		      Named("last_d2lp",last_d2lp));
}
