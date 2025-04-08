##############################
### LINEAR TIME SERIES PROJECT
##############################


################
### I. SETUP ###
################

#Downloading packages (code borrowed online)
list.of.packages <- c("readr", "zoo", "tseries", "stargazer", "fUnitRoots", "dplyr",
                      "aTSA", "xtable", "forecast", "ellipse", "graphics", "knitr", 
                      "rmarkdown", "markdown","rstudioapi")

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
series <- read_csv2('valeurs_mensuelles_jus.csv')

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
ddproduction <- diff(dproduction,1)

########################
### PLOTTING BOTH SERIES
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

#We define two functions to determine the number of lags we have to consider to get rid of autocorrelation

#Qtest
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

#adfTest
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

Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients)) 
#Tout est autocorrélé donc il va falloir tout faire les petits lags.
# NOTE A MOI MEME -> check le box test et H0 pour interpréter

#Donc on corrige pour 24 lag
adfTest_valid(production_train, 24, "c") 
# P VALUE 0.3937 -> On rejette PAS H0 du test ADF, donc c'est pas stationnaire c'est parfait,on va devoir différencier.

#We will thus differenciate our serie to obtain a better rejected adf test

#Creating lags and first difference
series$lag1 <- lag(series$Valeur)
series$dif <- series$Valeur - series$lag1

#Again, we run a regression to test for the presence of a time trend
reg_diff <- lm(dif ~ Time, data=series)
stargazer(reg_diff, type="text")
# We find no time trned, this is perfection.

#There is no time trend and significant constant so we run an adf test of type "nc"
dadf <- adfTest(dproduction, lag=0, type="nc") 
dadf

#We reject the non stationarity but we didn't consider lags so the test isn't valid
#Again, all the residuals are correlated if we use an adf test with 0 lags

Qtests(dadf@test$lm$residuals, 24, fitdf = length(dadf@test$lm$coefficients)) 
#Oh no, there is autocorrelation
#We run the same test as before to determine the number of lags we have to consider

adfTest_valid(dproduction,24,"nc")

#Do the pp test to be a good fayot
pp.test(x=as.vector(dproduction), output=TRUE) #Phillips-Perron test
kpss.test(x=as.vector(dproduction)) #KPSS

######################
### IV ARMA MODELS ###
######################

#We center the series just to test because we did it in the TD
y <- dproduction - mean(dproduction)
par(mfrow=c(1,2))
acf(y,24);pacf(y,24)
dev.off()

#Judging by the graph we then choose
qmax <- 1
pmax <- 5

arima <- arima(y,c(1,0,3))
Box.test(arima$residuals, lag=6, type="Ljung-Box", fitdf=5) #

Qtests(arima$residuals, 24, 5)

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
      try(arima(data, c(p, 0, q), optim.control = list(maxit = 20000)))
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

armamodels <- arma_model_choice(pmax,qmax,dproduction) #estime tous les arima (patienter...)

print(armamodels)
selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),] #modeles bien ajustes et valides
selec
### On a ? modeles bien ajustes et valides

### PAS DE MODELES BIEN AJUSTÉS ET VALIDES A LA FOIS -> PROBLEMATIQUE
## Peut-être que la Série n'est en fin de compte pas stationnaire ? On essaie avec la série différenciée deux fois juste pour être sûr...

par(mfrow=c(1,2))
acf(ddproduction,24);pacf(ddproduction,24)
dev.off()

plot(ddproduction)

pmax <- 7
qmax <- 2

armamodels <- arma_model_choice(pmax,qmax, ddproduction) #estime tous les arima (patienter...)

print(armamodels)
selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),] #modeles bien ajustes et valides
selec
