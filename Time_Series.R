##############################
### LINEAR TIME SERIES PROJECT
##############################


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
dates <- as.yearmon(seq(from=1990+0/12,to=2023+1/12,by=1/12)) 

#Transforming into a zoo object
production <- zoo(series$Valeur, order.by = dates)

#Dropping the last 4 values for prediction + creating differentiated series
production_train <- production[1:(length(production)-4)]
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

# -> The coefficient of Time is not significant, there is no trend.
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
# After 5 lags we got rid of autocorrelation.
# We get a p-value of 0.0376. The series thus could be stationary but the result is quite close to 0.05.

#We decide to investigate further by doing a PP test
pp.test(as.vector(production_train), output = TRUE)
#We get a p-value of 0.486, we are not able to reject HO and thus the series is celarly NOT stationary.
# Thus the need to differentiate our series.

#Creating lags and first difference in order to check for trend with a linear regression
series$lag1 <- lag(series$Valeur)
series$dif <- series$Valeur - series$lag1

#We run a regression to test for the presence of a time trend
reg_diff <- lm(dif ~ Time, data=series)
stargazer(reg_diff, type="text")
# We again find no trend, which is good.

#There is no time trend and no significant constant so we run an adf test of type "nc"
dadf <- adfTest(dproduction, lag=0, type="nc") 
dadf

#We reject the non stationarity but we didn't consider lags so the test isn't valid

#Testing again for autocorrelation
Qtests(dadf@test$lm$residuals, 24, fitdf = length(dadf@test$lm$coefficients)) 

#Same as before, as there is autocorrelation we run the test to determine the number of lags we have to consider
adfTest_valid(dproduction,24,"nc")
# After 4 lags we get a small p-value at less than 0.01. OUr series is then likely stationary?

#Just to be sure we can aso run PP + KPSS tests to further strenghten our result.
pp.test(x=as.vector(dproduction), output=TRUE) 
kpss.test(x=as.vector(dproduction)) #KPSS
# In both cases we reject (or not) H0 leading to the conclusion that our differentiated series is stationary.


#######################
### IV. ARMA MODELS ###
#######################

#We Begin by plotting ACF/PACF to determine the order of ARIMA we will consider
par(mfrow=c(1,2))
acf(dproduction,24);pacf(dproduction,24)
dev.off()

#Judging by the graph we then choose
qmax <- 4
pmax <- 5

# Computing an ARIMA(5,1,4) + Box test to check autocorrelation of residuals
arima <- arima(dproduction,c(5,1,4))
Qtests(arima$residuals, 24, 5)
# All p-values are higher than 0.05, which mean ARMA(5,4) would be a valid model

#We now check is the model is well adjusted:
arima
# We find that multiple coefficient are not significant. The model is thus not well adjusted, wo do not keep it

# We will now generalize this reasoning to every possible model.
# We then define 3 useful functions in order to test every combination of possible ARMA models.
p_value <- #This function returns the p-values of the estimated ARMA
  function(estim) {
    coef <- estim$coef 
    se <- sqrt(diag(estim$var.coef)) 
    t <- coef / se 
    pval <- (1 - pnorm(abs(t)))*2 
    return(rbind(coef, se, pval))
  }

model_choice <- #This function estimates an arima and checks the fit and validity of the model with p-value
  function(p, q, data = dproduction, k = 24) {
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
  function(pmax, qmax,data = dproduction) {
    pqs <- expand.grid(0:pmax, 0:qmax) 
    t(apply(matrix(1:dim(pqs)[1]), 1, function(row) { 
      p <- pqs[row, 1]
      q <- pqs[row, 2]
      cat(paste0("Computing ARMA(", p, ",", q, ") \n"))
      model_choice(p, q, data) 
    }))
  }

# We calculate all possible ARMA models
armamodels <- arma_model_choice(pmax,qmax,dproduction) #estime tous les arima (patienter...)

#We then select only models that are valid and well adjusted
selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),]
selec
### We get 4 well adjusted and valid models:
# ARMA(5,1)
# ARMA(1,2)
# ARMA(4,2)
# ARMA(0,3)

############
## AIC / BIC
## Now that we have 4 models, we need to determine which one we will choose using the AIC and BIC criterion.

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
  estim <- try(arima(dproduction,c(p,1,q),include.mean = F)) #tries ARIMA estimation
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic #assigns the AIC
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim) #assigns the BIC
}

#Print the AICs
AICs
AICs==min(AICs)

xtable(AICs) #We get the latex table

#Print the AICs
BICs
BICs==min(BICs)

xtable(BICs) #We get the latex table

# We get that ARIMA(5,1,1) minimizes AIC and ARIMA(1,1,2) minimizes the BIC
# Both models are also valid and well adjusted. We decide to keep them both for the moment.
arima511 <- arima(dproduction,c(5,1,1))
arima112 <- arima(dproduction,c(1,1,2))


# TODO: distinguish the best model by calculating R2


#######################
### V. MODELIZATION ###
#######################

# TODO: Attention a bien check sur quel séries on calcule les ARMA
# Typiquement peut être qu'on devrait le faire sur la série pas différenciée en fait...
