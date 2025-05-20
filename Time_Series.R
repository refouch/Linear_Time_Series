##################################
### LINEAR TIME SERIES PROJECT ###
##################################
# Rémi Fouchérand / Galeran Subileau


################
### I. SETUP ###
################

#Downloading packages (code borrowed online)
list.of.packages <- c("zoo", "tseries", "stargazer", "fUnitRoots", "dplyr",
                      "aTSA", "xtable", "forecast", "ellipse", "graphics", "knitr","rstudioapi")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

invisible(lapply(list.of.packages, library, character.only = TRUE))

#Set working directory to current
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Reset environment
rm(list=ls())

#######################################
### II. IMPORTING AND CLEANING DATA ###
#######################################

#Importing the series
series <- read.csv("monthly_values.csv", sep=";")

#Cleaning unwanted data and renaming the columns
series <- series[-c(1:3), ]
series <- series[, !sapply(series, function(col) all(col == 'A'))]
colnames(series) <- c('Mois','Valeur')
series$Valeur <- as.numeric(series$Valeur) #Converting the variables to the right type

#Inverting the time values
series <- series[dim(series)[1]:1,]

# Creating the Time index variable
series$Time <- 1:nrow(series)

#Creating dates vector for clean plotting.
dates <- as.yearmon(seq(from=1990+0/12,to=2025+1/12,by=1/12)) 

#Transforming into a zoo object
production <- zoo(series$Valeur, order.by = dates)

#Dropping the last 2 values for prediction + creating differentiated series
production_train <- production[1:(length(production)-2)]
dproduction <- diff(production_train,1)

########################
### PLOTTING BOTH SERIES

plot(production_train)

data_plot <- cbind(production_train,dproduction)
colnames(data_plot) <- c('Raw Series','Diff. Series')
plot(data_plot, 
     xlab = 'Time',
     main = "Series before and after differentiation",
     col = 'blue')

#########################
### III. STATIONARITY ###
#########################

#Regressing to check for a trend
rg <- lm(Valeur ~ Time, data = series)
stargazer::stargazer(rg, type = 'text')
#stargazer::stargazer(rg,type = 'latex')

# -> The coefficient of Time is not significant, there is no trend. But we remark a significant constant
# Doing the ADF test with only a constant.
adf <- adfTest(production_train, lag = 0, type = 'c')
adf

# P value is very small, we can reject the null hypothesis -> The series is likely stationary
# BUT: we need to check for autocorrelation of the residuals,

#We define functions to determine the number of lags we have to consider to get rid of autocorrelation

#Qtest
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients)) 
#Every p-value is very small, thus all residuals are autocorrelated !
#We then define a function to run the adf test with lags until residuals aren't correlated anymore.

#adfTest_valid
adfTest_valid <- function(series,kmax,type){ #ADF tests until no more autocorrelated residuals
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k, " lags: residuals OK? "))
    adf <- adfTest(series,lags=k,type=type)
    pvals <- Qtests(adf@test$lm$residuals,kmax,fitdf=length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T) == 0) {
      noautocorr <- 1; cat("OK \n")}
    else cat("No \n")
    k <- k + 1
  }
  return(adf)
}


adfTest_valid(production_train, 24, "c") 
# After 11 lags we got rid of autocorrelation -> can indicate weak result.
# We get a p-value of 0.01. The series thus could be stationary but 11 lags is a lot.

#We decide to investigate further by doing a PP test
pp.test(as.vector(production_train), output = TRUE)
#We get a p-value of 0.01, indicateing the series is also likely stationary
# We still perform a KPSS test for good measure
kpss.test(x=as.vector(production_train))
# We get a p-value of 0.0775 -> also indicates stationarity, but quite close to 0.05 treshold.

# We stil decide to differenciate for good measure
# Creating lags and first difference in order to check for trend with a linear regression
series$lag1 <- lag(series$Valeur)
series$dif <- series$Valeur - series$lag1

#We run a regression to test for the presence of a time trend
reg_diff <- lm(dif ~ Time, data=series)
stargazer(reg_diff, type="text")
#stargazer(reg_diff, type="latex")
# We again find no trend, and this time no constant

#There is no time trend and no significant constant so we run an adf test of type "nc"
dadf <- adfTest(dproduction, lag=0, type="nc") 
dadf

#Testing again for autocorrelation
Qtests(dadf@test$lm$residuals, 24, fitdf = length(dadf@test$lm$coefficients))
# All residuals are again correlated, the previous result is not valid.

#Same as before, as there is autocorrelation we run the test to determine the number of lags we have to consider
adfTest_valid(dproduction,24,"nc")
# After 4 lags we get a small p-value at less than 0.01. The series is likely stationary.

#Just to be sure we also run PP + KPSS tests.
pp.test(x=as.vector(dproduction), output=TRUE) 
kpss.test(x=as.vector(dproduction))#KPSS
# In both cases we reject (or not) H0 -> our differentiated series is stationary.


#######################
### IV. ARMA MODELS ###
#######################

#We Begin by plotting ACF/PACF to determine the order of ARIMA we will consider
acf(dproduction,24)
pacf(dproduction,24)


#Judging by the graph we then choose
qmax <- 2
pmax <- 5

# We will now compute every possible model and check for validity + well adjusted.
# We then define 3 useful functions in order to test every combination of possible ARMA models.

p_value <- #This function returns the p-values of the estimated ARMA -> to determine if it is well adjusted
  function(estim) {
    coef <- estim$coef 
    se <- sqrt(diag(estim$var.coef)) 
    t <- coef / se 
    pval <- (1 - pnorm(abs(t)))*2 
    return(rbind(coef, se, pval))
  }

model_choice <- #This function estimates an arima and checks the fit and validity of the model with p-value + Autocorrelation of residuals (Qtest)
  function(p, q, data = production_train, k = 24) {
    estim <-
      try(arima(data, c(p, 1, q), optim.control = list(maxit = 20000)))
    if (class(estim) == "try-error")
      return(c(
        "p" = p,
        "q" = q,
        "arsignif" = NA,
        "masignif" = NA,
        "resnocorr" = NA,
        "ok" = NA
      ))
    arsignif <- if (p == 0) 
      NA
    else
      p_value(estim)[3, p] <= 0.05 
    masignif <- if (q == 0) 
      NA
    else
      p_value(estim)[3, p + q] <= 0.05 
    resnocorr <-
      sum(Qtests(estim$residuals, 24, length(estim$coef) - 1)[, 2] <= 0.05, na.rm =
            T) == 0 
    checks <- c(arsignif, masignif, resnocorr) 
    ok <-
      as.numeric(sum(checks, na.rm = T) == (3 - sum(is.na(checks)))) 
    return(
      c(
        "p" = p,
        "q" = q,
        "arsignif" = arsignif,
        "masignif" = masignif,
        "resnocorr" = resnocorr,
        "ok" = ok
      )
    )
  }

arma_model_choice <- #This function runs the previous one with all p<pmax & q<qmax
  function(pmax, qmax,data = production_train) {
    pqs <- expand.grid(0:pmax, 0:qmax) 
    t(apply(matrix(1:dim(pqs)[1]), 1, function(row) { 
      p <- pqs[row, 1]
      q <- pqs[row, 2]
      cat(paste0("Computing ARMA(", p, ",", q, ") \n"))
      model_choice(p, q, data) 
    }))
  }

# We calculate all possible ARMA models
armamodels <- arma_model_choice(pmax,qmax,production_train) #estime tous les arima (patienter...)

#We then select only models that are valid and well adjusted
selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),]
selec
### We get 3 well adjusted and valid models:
# ARIMA(5,1,0)
# ARIMA(4,1,2)
# ARIMA(4,1,1)

############
## AIC / BIC
## Now that we have 3 models, we need to determine which one we will choose using the AIC and BIC.

#We create a loop to calculate the AIC/BIC for each possibility in the grid (qmax,pmax)=(2,5)
mat <- matrix(NA,nrow=pmax+1,ncol=qmax+1) #empty matrix to fill
rownames(mat) <- paste0("p=",0:pmax) #renames lines
colnames(mat) <- paste0("q=",0:qmax) #renames columns
AICs <- mat #AIC matrix not filled 
BICs <- mat #BIC matrix not filled 
pqs <- expand.grid(0:pmax,0:qmax) #all possible combinations of p and q
for (row in 1:dim(pqs)[1]){ #loop for each (p,q)
  p <- pqs[row,1] #gets p
  q <- pqs[row,2] #gets q
  estim <- try(arima(production_train,c(p,1,q),include.mean = F)) #tries ARIMA estimation
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic #assigns the AIC
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim) #assigns the BIC
}

#Print the AICs
AICs

AICs==min(AICs)

#xtable(AICs) #We get the latex table

#Print the AICs
BICs
BICs==min(BICs)

#xtable(BICs) #We get the latex table

# We get that ARIMA(5,1,0) minimizes AIC and ARIMA(1,1,1) minimizes the BIC
arima510 <- arima(production_train,c(5,1,0))
arima111 <- arima(production_train,c(1,1,1))

#ARIMA(5,1,0) is valid ! (We show tables for the report)
arima510
qtest_arima510 <- Qtests(arima510$residuals,24,fitdf= length(arima510$coef)-1)
xtable(qtest_arima510)

#ARIMA(1,1,1) is NOT valid ! (showing tables for the report)
arima111
qtest_arima111 <- Qtests(arima111$residuals,24,fitdf= length(arima111$coef)-1)
xtable(qtest_arima111)

#We take a look at the residuals of our arima, to see if they are gaussian (for prediction)
checkresiduals(arima510)

#We plot the inverse of our roots (to check the existence of a solution)
model <- Arima(production_train,c(5,1,0))
plot(model)
#All the roots of our model are > 1, there exists a causal solution.
#ARIMA(5,1,0) seems right for prediction.

#######################
### V. MODELIZATION ###
#######################

#We plot our model against the actuel series 
series$lag1 <- lag(series$Valeur)
series$lag2 <- lag(series$lag1)
series$lag3 <- lag(series$lag2)
series$lag4 <- lag(series$lag3)
series$lag5 <- lag(series$lag4)
series$lag6 <- lag(series$lag5)

phi1 <- arima510$coef[1]
phi2 <- arima510$coef[2]
phi3 <- arima510$coef[3]
phi4 <- arima510$coef[4]
phi5 <- arima510$coef[5]

series$pred <- (1+phi1)*series$lag1 + (phi2-phi1)*series$lag2 +
  (phi3-phi2)*series$lag3 + (phi4-phi3)*series$lag4 +
  (phi5-phi4)*series$lag5 - phi5*series$lag6

predict <- zoo(series$pred,order.by=dates)
predict <- predict[1:(length(predict)-2)]

plot(production_train, col = "black", xlab = "Time", ylab = "Values")
lines(predict, col = "red")
legend("topleft", legend = c(expression(X[t]), "ARIMA(5,1,0)"), col = c("black", "red"), lty = 1)

######################
### VI. PREDICTION ###
######################

#We supposed that our residuals are gaussian
#We want to estimate their standard error
residuals <- production - predict
residuals <- na.omit(residuals)
sigma <- sd(residuals)
sigma

#The standard error of our residuals is sigma=10.26202

#We seek to predict the two next values of our series thanks to the ARIMA model
arima510 <- arima(production_train,c(5,1,0))
model_pred <- predict(arima510, n.ahead=2) #prediction
pred <- zoo(model_pred$pred, order.by=as.yearmon(c(2025+0/12,2025+1/12))) #transforming our predictions into a zoo object
link <- rbind(production_train[length(production_train)],pred[1]) 

pred
# We get the value for bot predictions:
# T+1 = 77.21
# T+2 = 79.51

#We create a function that plots the prediction, starting the series at a given date
plot_pred <- function(start){
  plot(production, col = "black", ylab = "Values", xlab="Time",
       xlim = c(start, 2025+3/12))
  U <- model_pred$pred + 1.96*model_pred$se
  L <- model_pred$pred - 1.96*model_pred$se
  Upper1 <- model_pred$pred[1] + 1.96*sigma
  Lower1 <- model_pred$pred[1] - 1.96*sigma
  Upper2 <- model_pred$pred[2] + 1.96*sigma*sqrt((1+(1+phi1)^2))
  Lower2 <- model_pred$pred[2] - 1.96*sigma*sqrt((1+(1+phi1)^2))
  Upper <- c(Upper1,Upper2)
  Upper <- zoo(Upper, order.by=as.yearmon(c(2025+0/12,2025+1/12)))
  Lower <- c(Lower1,Lower2)
  xx <- c(time(Upper), rev(time(Upper)))
  yy <- c(Lower,rev(Upper))
  polygon(xx,yy,border="8", col=gray(0.6,alpha=0.2))
  lines(pred, type = "p", col = "darkgreen")
  lines(pred, type = "l", col = "darkgreen")
  lines(link, type = "l", col = "darkgreen")
  legend("topleft", legend=c(expression(X[t]), "Prediction"), col=c("black", "darkgreen"), lty=1)
}

#We plot our prediction, starting the series in January 2022
plot_pred(2022)

