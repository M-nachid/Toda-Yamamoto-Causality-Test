# Toda-Yamamoto-Causality-Test

This is an R function to perform the Toda-Yamamoto causality test (Toda & Yamamoto, 1995), a test of the null hypothesis than one time series does not "Granger-cause" another one. A time series X is considered a Granger cause of another time series Y if past values of X and Y predicts Y significantly better than past values of Y alone (Granger, 1969).

The concept and method of Granger causality was originally developed for stationary time series. The Toda-Yamamoto approach overcomes some issues that can arise when testing non-stationary series, a condition leading to the risk of spurious causality (He, Maekawa, 2001). Toda and Yamamoto propose that their method:

is applicable whether the VAR’s may be stationary (around a deterministic trend), integrated of an arbitrary order, or cointegrated of an arbitrary order. Consequently, one can test linear or nonlinear restrictions on the coefficients by estimating a levels VAR and applying the Wald criterion, paying little attention to the integration and cointegration properties of the time series data in hand. (Toda & Yamamoto, 1995, pp. 245-246).

The procedure is as follows (it is also clearly explained on the Econometrics Beat blog by prof. Dave Giles):

(...) we can estimate levels VAR’s and test general restrictions on the parameter matrices even if the processes may be integrated or cointegrated of an arbitrary order; we can apply the usual lag selection procedure (...) to a possibly integrated or cointegrated VAR (as far as the order of integration of the process does not exceed the true lag length of the model). Having chosen a lag length k, we then estimate a (k + equation)th-order VAR where equation is the maximal order of integration that we suspect might occur in the process. The coefficient matrices of the last equation lagged vectors in the model are ignored (since these are regarded as zeros), and we can test linear or nonlinear restrictions on the first k coefficient matrices using the standard asymptotic theory. (Toda & Yamamoto, 1995, p. 227).

An R code to run a bivariate version of the Toda Yamamoto test was published by Christoph Pfeiffer on his website. More recently, Josephine Lukito shared an R function to run a Toda-Yamamoto method for two-way Granger Causality in the ira_3media public repository for replicating results from her paper (Lukito, 2020).

The code I am posting here is a blended and modified version of the code shared by Pfeiffer and Lukito. It creates an R loop to apply the test to a multivariate VAR (it is therefore appliable both to bivariate and multivariate VARs), using the approach and functions suggested by Pfeiffer. The function can be applied to any VAR model and makes it easier and faster to run the analysis.

In particular, the function includes an authomated way to detect the maximum order of integration equation that must be added to the lag length k of the original VAR model to estimate the (k + equation)th-order VAR required by the Toda & Yamamoto procedure (Toda & Yamamoto, 1995). More specifically, the order of integration is tested by using the function ndiffs from the Rob J. Hyndman's forecast package. The function ndiffs uses a unit root test to determine the number of differences required for time series x to be made stationary. The default test is Kwiatkowski–Phillips–Schmidt–Shin (KPSS), but it can also be used the Augmented Dickey-Fuller test (test="adf") or the Phillips-Perron test (test="pp").

Besides forecast, the function requires the library vars, to fit a VAR model augmented of equation lags, and aod to perform the Wald Test ignoring the equation lags as required by the Toda-Yamamoto procedure.

The main argument of the function is a VAR model fitted with the R library vars.

The Toda-Yamamoto causality test is applied to the multivariate case following the equation for the direct Granger procedure with more than two variables, as reported in Kirchgässner & Wolters (2007, p. 114):

Let z_1, ..., z_m be additional variables. According to the definition of Granger causality, the estimation equation (3.21)

equation

can be extended to

equation

if we test for simple Granger causal relations, with equation, being the coefficients of the additional variables. It does not matter whether the additional variables are endogenous or exogenous since only lagged values are considered. After determining the numbers of lags equation, ..., (3.23) can be estimated using OLS. As in the bivariate case, it can be checked via an F test whether the coefficients of the lagged values of x are jointly significantly different from zero. By interchanging x and y in (3.23), it can be tested whether there exists a simple Granger causal relation from y to x and/or feedback.

Accordingly to the above equation, when there are equation additional time series, the test is performed by controlling for the effect of the equation series. Moreover, using the Toda-Yamamoto procedure, the coefficients of the lagged values of X are tested excluding the additional lags equation and using the Wald Chi-Squared Test instead of the F test.

# Example from the Christoph Pfeiffer's blog
list.packages  <-  c("fUnitRoots", "urca", "vars", "aod", "zoo", "tseries")
new.packages <- list.packages[!(list.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(fUnitRoots)
library(urca)
library(vars)
library(aod)
library(zoo)
library(tseries)

#Load data
cof <- read.csv("http://christophpfeiffer.org/wp-content/uploads/2012/11/coffee_data.csv", header=T,sep=";")
cof <- cof[complete.cases(cof),]
names(cof)[1] <- "Date"

#Adjust Date format
cof["Date"]<- paste(sub("M","-",cof$Date),"-01",sep="")

#Visualize
plot(as.Date(cof$Date),cof$Arabica,type="l",col="black",lwd=2)
lines(as.Date(cof$Date),cof$Robusta,col="blue",lty=2,lwd=1)
legend("topleft",c("Arabica","Robusta"),col=c("black","blue"),lty=c(1,2),lwd=c(2,1),bty="n")

#Possible structural break in 1970s. Therefore only values from 1976:01 onwards are regarded
cof1<-cof[193:615,]

#Visualize
plot(as.Date(cof1$Date),cof1$Arabica,type="l",col="black",lwd=2,ylim=range(cof1$Robusta))
lines(as.Date(cof1$Date),cof1$Robusta,col="blue",lty=2,lwd=1)
legend("topright",c("Arabica","Robusta"),col=c("black","blue"),lty=c(1,2),lwd=c(2,1),bty="n")

#Test for unit roots
adf.test(cof$Arabica)
adf.test(cof$Robusta)
kpss.test(cof$Arabica)
kpss.test(cof$Arabica)

adf.test(diff(cof$Arabica,1))
adf.test(diff(cof$Robusta,1))
kpss.test(diff(cof$Arabica,1))
kpss.test(diff(cof$Robusta,1))

# Since first order differencing eliminates the unit root, the maximum order of integration
# is concluded to be I(1).

#Set up VAR-Model
#select lag order // either 2 or 6
VARselect(cof1[,2:3],lag=20,type="both")

#VAR Model, lag=2
V.2<-VAR(cof1[,2:3],p=2,type="both")
serial.test(V.2)

#VAR-Model, lag=6
V.6<-VAR(cof1[,2:3],p=6,type="both") 
serial.test(V.6) #Stability analysis 1/roots(V.6)[[1]] # "&gt;1"
1/roots(V.6)[[2]] # "&gt;1"

#Alternative stability analyis
plot(stability(V.6)) ## looks fine

# Model with p=6 is less likely to be serially correlated. Thus model with p=6 is selected.

# Wald-test for the first 6 lags
# The test can be directly done with the VAR model, however using the correct
# variables is a little more tricky

#VAR-Model, lag=7 (additional lag, though not tested)
V.7<-VAR(cof1[,2:3],p=7,type="both")
V.7$varresult
summary(V.7)

#Wald-test (H0: Robusta does not Granger-cause Arabica)
wald.test(b=coef(V.7$varresult[[1]]), Sigma=vcov(V.7$varresult[[1]]), Terms=c(2,4,6,8,10,12))
# Could not be rejected (X2=8.6; p=0.2)

#Wald.test (H0: Arabica does not Granger-cause Robusta)
wald.test(b=coef(V.7$varresult[[2]]), Sigma=vcov(V.7$varresult[[2]]), Terms= c(1,3,5,7,9,11))
# Could be rejected at 10% (X2=12.3; p=0.056)

# It seems that Arabica Granger-causes Robusta prices, but not the other way around.

#################################################################################################################################
# TODA-YAMAMOTO CAUSALITY TEST ----
# Function to implement the Toda and Yamamoto version of Granger causality test.
# Toda, H. Y., & Yamamoto, T. (1995). Statistical inference in vector autoregressions with possibly integrated processes. 
# Journal of econometrics, 66(1-2), 225-250.
#################################################################################################################################

toda.yamamoto <- function(var.model, test = c("kpss","adf","pp")) {
  
  require(forecast)
  require(vars)
  require(aod)
  
  ty.df <- data.frame(var.model$y)
  ty.varnames <- colnames(ty.df)
  
  # estimates the number of first differences d_max required to make a given time series stationary. 
  d.max <- max(sapply(ty.df, function(x) forecast::ndiffs(x, test = test)))
  
  # k + d_max according to Toda & Yamamoto (1995)
  ty.lags <- var.model$p + d.max
  ty.augmented_var <- vars::VAR(ty.df, ty.lags, type=var.model$type)
  
  ty.results <- data.frame(cause = character(0), 
                           effect = character(0), 
                           chisq = numeric(0), 
                           pvalue = numeric(0))
  
  
  for (k in 1:length(ty.varnames)) {
    for (j in 1:length(ty.varnames)) {
      if (k != j){
        
        # coefficients to test, ignoring the d_max lags
        ty.coefres <- head(grep(ty.varnames[j], 
                                setdiff(colnames(ty.augmented_var$datamat), 
                                        colnames(ty.augmented_var$y))))
        
        wald.res <- wald.test(b=coef(ty.augmented_var$varresult[[k]]), 
                              Sigma=vcov(ty.augmented_var$varresult[[k]]),
                              Terms = ty.coefres) 
        
        ty.results <- rbind(ty.results, data.frame(
          cause = ty.varnames[j], 
          effect = ty.varnames[k], 
          chisq = as.numeric(wald.res$result$chi2[1]),
          pvalue = as.numeric(wald.res$result$chi2[3]))
        )
      }
    }
  }
  return(ty.results)
}

# -------------------------------------
# Using the toda.yamamoto function ####
toda.yamamoto(V.6)

For any comments or observation, you can drop me a message

References

Granger, C. W. (1969). Investigating causal relations by econometric models and cross-spectral methods. Econometrica: journal of the Econometric Society, 424-438.

He, Z., & Maekawa, K. (2001). On spurious Granger causality. Economics Letters, 73(3), 307-313.

Lukito, J. (2020). Coordinating a multi-platform disinformation campaign: Internet Research Agency activity on three US social media platforms, 2015 to 2017. Political Communication, 37(2), 238-255.

Kirchgässner, G., & Wolters, J. (2007). Introduction to Modern Time Series Analysis. Springer.

Toda, H. Y., & Yamamoto, T. (1995). Statistical inference in vector autoregressions with possibly integrated processes. Journal of econometrics, 66(1-2), 225-250.
