MD_OLS <-
  function(type_inference, level_significance, null_value, X, Y, V, W)	{

    n_obs<-NROW(X)
    p_dim<-NCOL(X)

    if (missing(V)){V <- NA
    }else{
      n_obs_1<-NROW(V)
      p_dim_1<-NCOL(V)
    }

    if (missing(W)){W <- NA}
    if (missing(level_significance)){level_significance<-0.05}

    exam_value<-c(type_inference %in% c(0:7))


    if(p_dim >= n_obs){

      print("The programm can't be completed since p > n")

    }else if(NROW(Y)!=n_obs){

      print("The programm can't be completed since row dimension of Y isn't consistent to that of X")

    }else if(exam_value=="FALSE"){

      print("Please input suitable range of type_inference")

    }else if(level_significance>1 || level_significance<0){

      print("Please input suitable range of level_significance")

    }else{

      tau <- p_dim / n_obs

      tolerance <- 1e-5

      upper_quantile = stats::qnorm(1 - level_significance / 2)



      XtX <- t(X)%*%X
      VtV <- t(V)%*%V

      inv_XtX <- solve(XtX)
      inv_VtV <- solve(VtV)

      trace_inv_XtX <- sum(diag(inv_XtX))
      trace_inv_VtV <- sum(diag(inv_VtV))

      trace_inv_XtX_square  <- sum(diag(inv_XtX%*%inv_XtX))
      trace_inv_VtV_square  <- sum(diag(inv_VtV%*%inv_VtV))

      beta_hat  = solve(XtX)%*%t(X)%*%Y
      gamma_hat = solve(VtV) %*% t(V)%*%W

      residual_1 <- Y - X%*%beta_hat
      residual_2 <- W - V %*%gamma_hat

      sigma_square_hat <- (sum((residual_1)^2)) / (n_obs - p_dim)
      delta_square_hat <- (sum((residual_2)^2)) / (n_obs - p_dim)

      nu_4_hat <- (1 - p_dim / n_obs) ^ (-4) * (mean((residual_1)^ 4) - 3 * sigma_square_hat^2 * tau * (1 - tau)^2 * (2 - tau) )


      if(type_inference == 0) # CI for quadratic functional
      {

        eta_hat_proposed <- t(beta_hat)%*%(XtX / n_obs) %*% beta_hat  - sigma_square_hat * tau
        # First formula in page 15

        var_proposed_test  <- n_obs ^ (-1) * (mean(Y^4) - nu_4_hat- 2*sigma_square_hat * eta_hat_proposed- eta_hat_proposed ^ 2 + 2 * sigma_square_hat ^ 2 * tau / (1 - tau))
        # Denominator of formula (3.4)


        if(var_proposed_test<tolerance){

          print("var_proposed_test is lower than tolerable level we allow")

          var_proposed_test<-NA
        }


        CI_proposed_method <- matrix(eta_hat_proposed,nrow=2,ncol=1) +c((var_proposed_test)^(1 / 2))* seq(-upper_quantile, upper_quantile,length=2)
        # Interval

        test_stat_proposed  <- (var_proposed_test )^(-1/2)*(eta_hat_proposed - null_value)
        # Formula (3.4)

        Result<- matrix(c(eta_hat_proposed,var_proposed_test,CI_proposed_method,test_stat_proposed),nrow=5,ncol=1)
        rownames(Result) <- c("eta_hat_proposed","var_proposed_test","CI_proposed_method_lower","CI_proposed_method_upper","test_stat_proposed")

      }
      if(type_inference == 1) #signal detection
      {
        var_proposed_test  <- 2 * sigma_square_hat^2 * trace_inv_XtX_square + 2 * sigma_square_hat^2 / (n_obs - p_dim) * trace_inv_XtX^2
        # bottom of page 12

        var_conventional_test  <- 4 * sigma_square_hat* t(beta_hat) %*% solve(XtX) %*% beta_hat
        # Theorem 4
        if(var_proposed_test<tolerance){
          print("var_proposed_test is lower than tolerable level we allow")
          var_proposed_test<-NA
        }
        if(var_conventional_test <tolerance){
          print("var_conventional_test is lower than tolerable level we allow")
          var_conventional_test<-NA
        }

        test_stat_proposed  <- (var_proposed_test )^(-1/2) * ((sum(beta_hat^2)) - trace_inv_XtX * sigma_square_hat - null_value)
        # Theorem 4

        test_stat_conventional  <- (var_conventional_test )^(-1/2) * (sum(beta_hat^2) -null_value)

        p_value_proposed  <- 1 - stats::pnorm(test_stat_proposed)
        p_value_conventional  <- 1 - stats::pnorm(test_stat_conventional)

        Result<- matrix(c(var_proposed_test,var_conventional_test,test_stat_proposed,test_stat_conventional,p_value_proposed,p_value_conventional, NA, NA),nrow=8,ncol=1)
        rownames(Result) <- c("var_proposed_test","var_conventional_test","test_stat_proposed","test_stat_conventional","p_value_proposed","p_value_conventional", "CI_proposed_method_lower","CI_proposed_method_upper")

      }
      if (type_inference == 2) # FVE
      {
        eta_hat_proposed <-(t(beta_hat) %*% (XtX / n_obs) %*% beta_hat - sigma_square_hat) * tau
        rho_hat_proposed <- eta_hat_proposed /(eta_hat_proposed + sigma_square_hat)
        # First and second formula in page 15

        var_proposed_test  <- n_obs ^(-1) *(eta_hat_proposed + sigma_square_hat)^(-4)*    (2 * sigma_square_hat^4 * p_dim / (n_obs - p_dim)-(2 + 4 * p_dim / (p_dim - n_obs))*sigma_square_hat^3 * eta_hat_proposed     + sigma_square_hat^2 * ( mean(Y^4) - nu_4_hat + eta_hat_proposed^2 * (4 * tau - 2) / (1-tau) )+ eta_hat_proposed ^ 2 * nu_4_hat)
        # Denominator of formula in theorem 5

        var_conventional_test  <- n_obs ^ (-1) *(t(beta_hat) %*% (XtX / n_obs) %*% beta_hat + sigma_square_hat)^(-4)*    ( (-2)*sigma_square_hat^3*(t(beta_hat) %*%( XtX / n_obs) %*% beta_hat)+ sigma_square_hat^2 * ( mean(Y^4) - mean( residual_1^4)))-2 * (t(beta_hat) %*% (XtX / n_obs) %*% beta_hat^2 ) + (t(beta_hat) %*% (XtX / n_obs) %*% beta_hat)^2 * mean(residual_1^4)
        # Denominator of formula (3.3) in page 47

        if(var_proposed_test<tolerance){
          print("var_proposed_test is lower than tolerable level we allow")
          var_proposed_test<-NA
        }
        if(var_conventional_test <tolerance){
          print("var_conventional_test is lower than tolerable level we allow")
          var_conventional_test<-NA
        }

        test_stat_proposed  <- (var_proposed_test )^(-1/2)*(rho_hat_proposed - null_value)
        # Formula in theorem 5

        test_stat_conventional  <- (var_conventional_test )^(-1/2) *(t(beta_hat) %*% (t(X)%*% X / n_obs) %*% beta_hat ) /( (t(beta_hat) %*%  (t(X)%*% X / n_obs) %*% beta_hat + sigma_square_hat)- null_value)
        # Formula (3.3) in page 47

        p_value_proposed  <- 1 - stats::pnorm(test_stat_proposed)

        p_value_conventional  <- 1 - stats::pnorm(test_stat_conventional)



        Result<- matrix(c(eta_hat_proposed,rho_hat_proposed,var_proposed_test,var_conventional_test,test_stat_proposed,test_stat_conventional,p_value_proposed,p_value_conventional, NA, NA),nrow=10,ncol=1)
        rownames(Result) <- c("eta_hat_proposed","rho_hat_proposed","var_proposed_test","var_conventional_test","test_stat_proposed","test_stat_conventional","p_value_proposed","p_value_conventional", "CI_proposed_method_lower","CI_proposed_method_upper")

      }
      if (type_inference == 3) # error variance
      {
        var_proposed_test  <- (nu_4_hat +sigma_square_hat^2 * (3 * tau - 1) / (1 - tau) ) / n_obs
        # Proposition 1

        if(var_proposed_test<tolerance){
          print("var_proposed_test is lower than tolerable level we allow")
          var_proposed_test<-NA
        }

        test_stat_proposed  <- (var_proposed_test )^(-1/2) * (sigma_square_hat - null_value)
        # Proposition 1

        p_value_proposed  <- 2 * stats::pnorm(- abs(test_stat_proposed ))

        Result<- matrix(c(var_proposed_test,test_stat_proposed,p_value_proposed, NA, NA),nrow=5,ncol=1)
        rownames(Result) <- c("var_proposed_test","test_stat_proposed","p_value_proposed", "CI_proposed_method_lower","CI_proposed_method_upper")

      }
      if (type_inference == 4) # equality
      {

        if(p_dim!=p_dim_1){

          print("The programm can't be completed since column dimension of V isn't consistent to that of X")

        }else{

          trace_inv_XtX_inv_VtV <- sum(diag(inv_XtX%*%solve(VtV)))

          var_proposed_test  <-    2 * sigma_square_hat^2 * (-trace_inv_XtX_square+ 1 / (n_obs - p_dim) * trace_inv_XtX^2)  + 2 *  delta_square_hat^2 * (- trace_inv_VtV_square + 1 / (n_obs_1 - p_dim) * trace_inv_VtV^2)  - 4 * sigma_square_hat * delta_square_hat * trace_inv_XtX_inv_VtV + 4 * sigma_square_hat * (t(beta_hat - gamma_hat) %*% solve(XtX) %*% (beta_hat - gamma_hat))+ 4 * delta_square_hat * (t(beta_hat - gamma_hat)%*% solve(VtV) %*% (beta_hat - gamma_hat))
          # Fifth formula in page 48

          if(var_proposed_test<tolerance){
            print("var_proposed_test is lower than tolerable level we allow")
            var_proposed_test<-NA
          }

          test_stat_proposed  <- (var_proposed_test )^(-1/2) *
            (sum((beta_hat - gamma_hat)^2) - trace_inv_XtX  * sigma_square_hat -  trace_inv_VtV * delta_square_hat)
          # Theorem S.2

          p_value_proposed  <- 1 - stats::pnorm(test_stat_proposed)

          Result<- matrix(c(var_proposed_test,test_stat_proposed ,p_value_proposed, NA, NA),nrow=5,ncol=1)
          rownames(Result) <- c("var_proposed_test","test_stat_proposed ","p_value_proposed", "CI_proposed_method_lower","CI_proposed_method_upper")

        }
      }





      if (type_inference == 5) # co-heritability
      {

        if(p_dim!=p_dim_1){

          print("The programm can't be completed since column dimension of V isn't consistent to that of X")

        }else{

          trace_inv_XtX_inv_VtV <- sum(diag(inv_XtX%*%solve(VtV)))

          norm_beta_estimated <-  (max(sum(beta_hat^2) - trace_inv_XtX * sigma_square_hat, tolerance))^(1/2)
          # Proposition S.2

          norm_gamma_estimated <- (max(sum(gamma_hat^2) - trace_inv_VtV * delta_square_hat, tolerance))^(1/2)
          # Proposition S.2

          co_herit_hat_proposed <- t(beta_hat) %*% gamma_hat/ (norm_beta_estimated * norm_gamma_estimated)
          # Proposition S.2

          var_proposed_test  <- -sigma_square_hat * delta_square_hat * trace_inv_XtX_inv_VtV/(norm_beta_estimated * norm_gamma_estimated)^2 + 1 / (norm_beta_estimated * norm_gamma_estimated)^2 *    delta_square_hat * t(beta_hat - gamma_hat %*% t(gamma_hat)%*% beta_hat / norm_gamma_estimated^2 )%*% solve(VtV) %*% (beta_hat - gamma_hat %*% t(gamma_hat) %*%beta_hat  / norm_gamma_estimated^2 )  +1/(norm_beta_estimated*norm_gamma_estimated)^2*sigma_square_hat*t(gamma_hat - beta_hat %*% t(beta_hat)%*%gamma_hat/norm_beta_estimated^2 )%*% solve(XtX) %*% (gamma_hat - beta_hat %*% t(beta_hat)%*% gamma_hat/ norm_beta_estimated^2 )+ (t(gamma_hat) %*% beta_hat)^2/(2 * norm_beta_estimated * norm_gamma_estimated^3)^2 * 2 * delta_square_hat^2 * ( - trace_inv_VtV_square + 1/(n_obs_1 - p_dim) * trace_inv_VtV^2 )  + (t(gamma_hat) %*% beta_hat)^2/ (2 * norm_beta_estimated^3 * norm_gamma_estimated)^2 * 2 * sigma_square_hat^2 *    ( - trace_inv_XtX_square + 1 / (n_obs - p_dim) *  trace_inv_XtX^2 )
          # Proposition S.2

          var_conventional_test  <-  1/(sum(beta_hat^2)) * (sum(gamma_hat^2)) * delta_square_hat * t(beta_hat - gamma_hat %*% t(gamma_hat)%*%beta_hat / (sum(beta_hat^2)) ) %*% solve(VtV) %*% (beta_hat - gamma_hat%*% t(gamma_hat)%*%beta_hat/(sum(gamma_hat^2)) )+ 1 / (sum(beta_hat^2))  * (sum(gamma_hat^2))*sigma_square_hat * t(gamma_hat - beta_hat %*% t(beta_hat) %*%gamma_hat / (sum(gamma_hat^2))) %*% solve(XtX) %*% (gamma_hat - beta_hat %*% t(beta_hat) %*% gamma_hat /(sum(beta_hat^2)))

          if(var_proposed_test<tolerance){
            print("var_proposed_test is lower than tolerable level we allow")
            var_proposed_test<-NA
          }
          if(var_conventional_test <tolerance){
            print("var_conventional_test is lower than tolerable level we allow")
            var_conventional_test<-NA
          }

          test_stat_proposed  <- (var_proposed_test )^(-1/2)* (co_herit_hat_proposed - null_value)
          # Proposition S.2

          co_herit_hat_conventional <- t(beta_hat) %*% gamma_hat  / (sqrt(sum(beta_hat^2)) * sqrt(sum(gamma_hat^2)))
          # Proposition S.2


          test_stat_conventional <-(var_conventional_test )^(-1/2)* (co_herit_hat_conventional - null_value)

          p_value_proposed <- 2 * stats::pnorm(-abs(test_stat_proposed ))

          p_value_conventional  <- 2 * stats::pnorm(-abs(test_stat_conventional ))

          Result<- matrix(c(norm_beta_estimated,norm_gamma_estimated,co_herit_hat_proposed,co_herit_hat_conventional,var_proposed_test,var_conventional_test,test_stat_proposed,test_stat_conventional,p_value_proposed,p_value_conventional, NA, NA),nrow=12,ncol=1)
          rownames(Result) <- c("norm_beta_estimated","norm_gamma_estimated","co_herit_hat_proposed","co_herit_hat_conventional","var_proposed_test","var_conventional_test","test_stat_proposed","test_stat_conventional","p_value_proposed","p_value_conventional", "CI_proposed_method_lower","CI_proposed_method_upper")

        }
      }
      if (type_inference == 6) # confidence ball
      {

        sigma_star_hat_square <- 2 * sigma_square_hat^2 * trace_inv_XtX_square+ 2 * sigma_square_hat^2 / (n_obs - p_dim)* trace_inv_XtX ^ 2
        # Formula (3.1)

        ID_correctly_covered_by_confidence_ball_two_sided  <-    abs(sum((beta_hat - null_value)^2) - trace_inv_XtX * sigma_square_hat) <= upper_quantile * sigma_star_hat_square^(1/2);
        # Formula (3.2)

        ID_correctly_covered_by_confidence_ball_one_sided  =    sum((beta_hat - null_value)^2) <= trace_inv_XtX * sigma_square_hat+ stats::qnorm(1 - level_significance) * sigma_star_hat_square^(1/2)
        # Formula (3.2)


        Result<- matrix(c(sigma_star_hat_square,ID_correctly_covered_by_confidence_ball_two_sided,ID_correctly_covered_by_confidence_ball_one_sided, NA, NA),nrow=5,ncol=1)
        rownames(Result) <- c("sigma_star_hat_square","ID_correctly_covered_by_confidence_ball_two_sided","ID_correctly_covered_by_confidence_ball_one_sided", "CI_proposed_method_lower","CI_proposed_method_upper")

      }
      if (type_inference ==7){

        var_proposed_test<- 4 * sigma_square_hat * (t(beta_hat) %*% solve(t(X)%*%X) %*% beta_hat) - 2 * sigma_square_hat^2 * trace_inv_XtX_square + 2 * sigma_square_hat^2 * trace_inv_XtX^2 * (n_obs-p_dim)^(-1)

        if(var_proposed_test<tolerance)
        {

          print("var_proposed_test is lower than tolerable level we allow")

          var_proposed_test<-NA

        }



        test_stat_proposed  <- (var_proposed_test )^(-1/2)* (sum(beta_hat^2) - trace_inv_XtX * sigma_square_hat - null_value)
        # Theorem S.2

        p_value_proposed  <- 1 - stats::pnorm(test_stat_proposed)

        Result<- matrix(c(var_proposed_test,test_stat_proposed,p_value_proposed, NA, NA),nrow=5,ncol=1)
        rownames(Result) <- c("var_proposed_test","test_stat_proposed","p_value_proposed", "CI_proposed_method_lower","CI_proposed_method_upper")

      }
      Result
    }
  }
