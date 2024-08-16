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
