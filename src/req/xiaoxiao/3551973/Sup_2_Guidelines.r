## Loading required packages. May need to be installed from CRAN.
library(nlrwr)
library(boot)
## Function to carry out the analytical process described in General Guidelines
## Input: x - vector, explanatory variable
##        y - vector, response variable
##        CI_boot - whether confidence intervals should be computed using bootstrapping
##                  for model averaging. Default is TRUE.
##        diagno - residual plots as visual diagnostics to check the assumptions of 
##                 normality and homoscedasticity for LR and NLR. Default is TRUE. 
##        output_plot - whether a scatter plot with fitted relationship is desired.
##                      Default is FALSE.
## Output: method - method used in analysis
##         a, b - estimated parameters
##         a_confint, b_confint - 95% confidence intervals for a & b
##                                (optional if method is "Model Averaging")
power_analysis = function (x, y, CI_boot = TRUE, diagno = TRUE, output_plot = FALSE){
  ## Step 1: Likelihood analysis
  model_lr = lm(log(y) ~ log(x))
  a_lr = exp(coef(summary(model_lr))[1, 1])
  b_lr = coef(summary(model_lr))[2, 1]
  sd_lr = sd(log(y) - (log(a_lr) + b_lr * log(x)))
  
  model_nlr = nls(y ~ a1 * x ^ a2, start = list(a1 = a_lr, a2 = b_lr),
    control = nls.control(maxiter = 2000, warnOnly = TRUE))
  a_nlr = coef(summary(model_nlr))[1, 1]
  b_nlr = coef(summary(model_nlr))[2, 1]
  sd_nlr = sd(y - a_nlr * x ^ b_nlr)

  l_logn = sum(log(dlnorm(y, log(a_lr * x ^ b_lr), sd_lr)))
  l_norm = sum(log(dnorm(y, a_nlr * x ^ b_nlr, sd_nlr)))

  n = length(x)
  k = 3
  AICc_logn = 2 * k - 2 * l_logn + 2 * k * (k + 1) / (n - k - 1)
  AICc_norm = 2 * k - 2 * l_norm + 2 * k * (k + 1) / (n - k - 1)
  delta_AICc = AICc_norm - AICc_logn
  writeLines(paste("AICc_logn: ", AICc_logn, "\nAICc_norm: ", AICc_norm))
  w_logn = exp(-(AICc_logn - min(AICc_logn, AICc_norm)) / 2)
  w_norm = exp(-(AICc_norm - min(AICc_logn, AICc_norm)) / 2)
  weight_logn = w_logn / (w_logn + w_norm)
  weight_norm = w_norm / (w_logn + w_norm)

  ## Step 2a: Analysis with NLR
  if (delta_AICc < - 2){
    writeLines("The assumption of additive normal error is better supported. \nProceed with NLR.")
    method = "NLR"
    a = a_nlr
    b = b_nlr
    a_confint = confint2(model_nlr)[1, ]
    b_confint = confint2(model_nlr)[2, ]
  }
  ## Step 2b: Analysis with LR
  else if (delta_AICc > 2){
    writeLines("The assumption of multiplicative log-normal error is better supported. \nProceed with LR.")
    method = "LR"
    a = a_lr
    b = b_lr
    a_confint = confint(model_lr)[1, ]
    b_confint = confint(model_lr)[2, ]
  }
  ## Step 2c: Analysis with model averaging
  else {
    writeLines("The two error distributions have similar support. \nProceed with model averaging.")
    method = "Model Averaging"
    a = a_lr * weight_logn + a_nlr * weight_norm
    b = b_lr * weight_logn + b_nlr * weight_norm
    if (!CI_boot){
      a_confint = NA
      b_confint = NA
    }
    else {
      boot.est=function(dat, indices) {
        dat.sub=dat[indices, ]
        names(dat.sub) = c("x", "y")
        model.lr = lm(log(y) ~ log(x), dat = dat.sub)
        a.lr = exp(coef(summary(model.lr))[1, 1])
        b.lr = coef(summary(model.lr))[2, 1]
        sd.lr = sd(log(dat.sub$y) - (log(a.lr) + b.lr * log(dat.sub$x)))
        a.lr.CI = confint(model.lr)[1, ]
        b.lr.CI = confint(model.lr)[2, ]
        model.nlr = nls(y ~ a1 * x ^ a2, start = list(a1 = a.lr, a2 = b.lr), dat = dat.sub,
         control = nls.control(maxiter = 2000, warnOnly = TRUE))
        a.nlr = coef(summary(model.nlr))[1, 1]
        b.nlr = coef(summary(model.nlr))[2, 1]
        sd.nlr = sd(dat.sub$y - a.nlr * dat.sub$x ^ b.nlr)
        a.nlr.CI = confint2(model.nlr)[1, ]
        b.nlr.CI = confint2(model.nlr)[2, ]

        l.logn = sum(log(dlnorm(dat.sub$y, log(a_lr * dat.sub$x ^ b_lr), sd_lr)))
        l.norm = sum(log(dnorm(dat.sub$y, a_nlr * dat.sub$x ^ b_nlr, sd_nlr)))
        AICc.logn = 2 * k - 2 * l.logn + 2 * k * (k + 1) / (n - k - 1)
        AICc.norm = 2 * k - 2 * l.norm + 2 * k * (k + 1) / (n - k - 1)
        AICc.min = min(AICc.logn, AICc.norm)
        weight.logn = exp(-(AICc.logn - AICc.min)/2)
        weight.norm = exp(-(AICc.norm - AICc.min)/2)
        logn.w = weight.logn / (weight.logn + weight.norm)
        norm.w = weight.norm / (weight.logn + weight.norm)

        a.boot = a.lr * logn.w + a.nlr * norm.w
        b.boot = b.lr * logn.w + b.nlr * norm.w
        return(c(a.boot, b.boot))
      }
      dat.boot=boot(data = as.data.frame(cbind(x, y)), statistic = boot.est, R = 1000)
      a_confint = boot.ci(dat.boot, index = 1, type = "perc")$perc[4:5]
      b_confint = boot.ci(dat.boot, index = 2, type = "perc")$perc[4:5]
    }
  }
  writeLines(paste("a: ", a, "\nb: ", b))

  ## Step 3: residual plots as visual diagnostics
  if (diagno){
    lr_res = resid(model_lr)
    lr_pred = predict(model_lr)
    nlr_res = resid(model_nlr)
    nlr_pred = predict(model_nlr)
    par(mfrow = c(2, 2), oma = c(0, 4, 2, 0))
    hist(lr_res, freq = FALSE, main = "", xlab = "Residuals")
    curve(dnorm(x, mean(lr_res), sd(lr_res)), add = TRUE)
    plot(lr_pred, lr_res, xlab = "Predicted y", ylab = "Residuals")
    abline(h = 0)
    hist(nlr_res, freq = FALSE, main = "", xlab = "Residuals")
    curve(dnorm(x, mean(nlr_res), sd(nlr_res)), add = TRUE)
    plot(nlr_pred, nlr_res, xlab = "Predicted y", ylab = "Residuals")
    abline(h = 0)
    mtext("      Normality                                            Homoscedasticity", 
      cex = 1.5, side = 3, outer = TRUE)
    mtext("NLR                                        LR", cex = 1.5, side = 2, outer = TRUE)
  }
  if (output_plot){                                                                                                                   
    par(mfrow = c(1, 2), oma = c(0, 4, 2, 0))
    plot(x, y, log = "xy", xlab = "x", ylab = "y", pch = 20, main = "Logarithmic Scale")
    curve(a * x ^ b, add = TRUE, lty = "dashed")
    plot(x, y, xlab = "x", ylab = "y", pch = 20, main = "Arithmetic Scale")
    curve(a * x ^ b, add = TRUE, lty = "dashed")
    title(main = paste("Fitting Power-Law Data with ", method), outer = TRUE)
  }
  
  return (list(method = method, a = a, b = b, a_confint = a_confint, b_confint = b_confint))
}
