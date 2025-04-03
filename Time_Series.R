##############################
### LINEAR TIME SERIES PROJECT
##############################


################
### I. SETUP ###
################

#Downloading packages (code borrowed online)
list.of.packages <- c("readr", "zoo", "tseries", "stargazer", "fUnitRoots", "dplyr",
                      "aTSA", "xtable", "forecast", "ellipse", "graphics", "knitr", 
                      "rmarkdown", "markdown")

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
series <- read_csv2('valeurs_mensuelles.csv')

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
qmax <- 4
pmax <- 1

