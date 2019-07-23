
                  ##### GIACOMINI AND WHITE (2006) TEST #######
                  
                  
#Before simulating MC I have to set the coefficients of the simulated regressions to predefined values.
#All this values are given/calculated by given numbers in paper by Paye.

#Set the coefficients to prespecified values..
const2=0.27
coef_ar1_x=0.60
unconexp_x=const2/(1-coef_ar1_x)

#Set the autoregresseive coefficients to predetermined values.
const1=-0.61
coef_ar1_y=0.52
coef_ar2_y=0.22
unconexp_y=(const1+unconexp_x)/(1-coef_ar1_y-coef_ar2_y)

# To generate the error terms one need to specify the variances/std.dev of error terms.
var_et=0.092
var_vt=0.160
cov_et_vt=0.026
corr_et_vt=cov_et_vt/sqrt(var_et*var_vt)

# Generate the distribution according to standard normal distribution. I set n=2000
set.seed(123456)
zt1<-rnorm(2000,0,1)
zt2<-rnorm(2000,0,1)

et<-zt1*sqrt(var_et)
vt<-(zt1*corr_et_vt+zt2*sqrt(1-corr_et_vt))*sqrt(var_vt)


#For loop inside the prediction/estimation
set.seed(123456)
wndw = 20 #number of observations for the in-sample (estimation) analysis.
T = 120  #number of observations for the out-of-sample (forecast) analysis.
y <- matrix(NA, 2000,10000)  #Randomly generated dependent variable matrix.
y_sample <- matrix(NA,T+wndw,10000) #creating a matrix which incorporates both in-sample and out-of-sample analysis osbervations: T+wndw=150 observations.
x <-matrix(NA,2000,10000)
x_sample <-matrix(NA,T+wndw,10000)

for (ii in 1:10000) {
  # Generate first two datapoints	
  y[1,ii] <- unconexp_y
  y[2,ii] <- 1
  x[1,ii] <-unconexp_x
  x[2,ii] <-1
  #randomly generate the rest of the data for indep.var. starting from the third obs.
  for (tt in 3:2000) {
    x[tt,ii] = const2+coef_ar1_x*x[tt-1,ii]+rnorm(1)*sqrt(var_vt)
    y[tt,ii] = const1 + coef_ar1_y*y[tt-1,ii] + coef_ar2_y*y[tt-2,ii] + rnorm(1)*sqrt(var_et)
  }
  #Select last 150 observations for later in- and out-of-sample analysis and assign them to predefined y_sample matrix.
  y_sample[,ii] <- y[(2000 - wndw - T + 1) : 2000, ii]
  x_sample[,ii] <- x[(2000 - wndw - T + 1) : 2000, ii]
}


# Estimate model on rolling window 
prediction <- matrix(NA,T+wndw,10000) #create an empty matrix later to be filled by forecasting values. First 20 observations are not gonna be used, only T=120 observations are used.
losses<-matrix(NA,T+wndw,10000) #This command generates a matrix for the forecast errors where it is calculated by subtracting actual generated data from the forecasted values.
losses_sqr<-matrix(NA,T+wndw,10000) #This command generates a matrix for the forecast error squares where it is calculated by taking the square of 'losses' series.
prediction4 <- matrix(NA,T+wndw,10000)
losses4<-matrix(NA,T+wndw,10000)
losses4_sqr<-matrix(NA,T+wndw,10000)
for (ii in 1:10000) {
  for (kk in wndw:(T+wndw-1)) {
    # estimate with data from (kk-wndw+1) : kk.Model is defined as AR(2) model and estimated with linear regression.
    fit2<-lm(y_sample[(kk-wndw+1): kk,ii]~dplyr::lead(y_sample[(kk-wndw+1): kk,ii])+dplyr::lead(y_sample[(kk-wndw+1): kk,ii],2))
    #Forecasts are obtained by extracting the coefficients from the regression model of AR(2).
    prediction[kk+1,ii]<-fit2$coefficients[1]+fit2$coefficients[2]*y_sample[kk,ii]+fit2$coefficients[3]*y_sample[kk-1,ii]
    #Forecast errors are the difference between the actual and predicted(forecast) values.
    losses[kk+1,ii] <- y_sample[kk+1,ii] - prediction[kk+1,ii]
    losses_sqr[kk+1,ii] <- (y_sample[kk+1,ii] - prediction[kk+1,ii])^2
    
    fit<-lm(y_sample[(kk-wndw+1): kk,ii]~dplyr::lead(y_sample[(kk-wndw+1): kk,ii])+dplyr::lead(y_sample[(kk-wndw+1): kk,ii],2)+dplyr::lead(x_sample[(kk-wndw+1): kk,ii]))
    prediction4[kk+1,ii]<-fit2$coefficients[1]+fit2$coefficients[2]*y_sample[kk,ii]+fit$coefficients[3]*y_sample[kk-1,ii]+fit$coefficients[4]*x_sample[kk,ii]
    losses4[kk+1,ii] <- y_sample[kk+1,ii] - prediction4[kk+1,ii]
    losses4_sqr[kk+1,ii] <- (y_sample[kk+1,ii] - prediction4[kk+1,ii])^2
    
  }
}




gw_stats1<-matrix(NA,T+wndw,10000)
gw_stats1<-(losses_sqr-losses4_sqr)*100

count<-matrix(NA,1,10000)

for (ii in 1:10000) {
  gw_stats1_eq<-lm(gw_stats1[21:140,ii]~1)
  if (summary(gw_stats1_eq)$coefficients[1, 4]<0.10){
    count[,ii]=1
  }else{
    count[,ii]=0  
  }
}

#This code calculates the analysis percentage of rejection rates (at 10% sign. level)
#over 10000 simulated samples. The test is based on Giacomini and White (2006) test.
GW_test_stats1<-100*(rowSums(count)/10000)



#Here I am simulating from the model that sets the coefficient for the xt-1 variable (new variable) to zero.
#Therefore, I create a new variable called x.  
wndw = 20
T = 120
set.seed(123456)

#Set the coefficients to prespecified values
const2=0.27
coef_ar1_x=0.60
unconexp_x=const2/(1-coef_ar1_x)

#Set the autoregresseive coefficients to predetermined values.
const1=-0.61
coef_ar1_y=0.52
coef_ar2_y=0.22
unconexp_y=(const1+unconexp_x)/(1-coef_ar1_y-coef_ar2_y)

# To generate the error terms one need to specify the variances/std.dev of error terms.
var_et=0.092
var_vt=0.160
cov_et_vt=0.026
corr_et_vt=cov_et_vt/sqrt(var_et*var_vt)

# Generate the distribution according to standard normal distribution. I set n=2000
set.seed(123456)
zt1<-rnorm(2000,0,1)
zt2<-rnorm(2000,0,1)

et<-zt1*sqrt(var_et)
vt<-(zt1*corr_et_vt+zt2*sqrt(1-corr_et_vt))*sqrt(var_vt)


#For loop inside the prediction/estimation
set.seed(123456)
wndw = 20
T = 120
y2 <- matrix(NA, 2000,10000)
y2_sample <- matrix(NA,T+wndw,10000)
x <-matrix(NA,2000,10000)
x_sample <-matrix(NA,T+wndw,10000)
for (ii in 1:10000) {
  # Generate datapoints	
  y2[1,ii] <- unconexp_y
  y2[2,ii] <- 1
  x[1,ii] <-unconexp_x
  x[2,ii] <-1
  for (tt in 3:2000) {
    x[tt,ii] = const2+coef_ar1_x*x[tt-1,ii]+rnorm(1)*sqrt(var_vt)
    y2[tt,ii] = const1 + coef_ar1_y*y2[tt-1,ii] + coef_ar2_y*y2[tt-2,ii] +0.09*x[tt-1,ii]+ rnorm(1)*sqrt(var_et)
  }
  y2_sample[,ii] <- y2[(2000 - wndw - T + 1) : 2000, ii]
  x_sample[,ii] <- x[(2000 - wndw - T + 1) : 2000, ii]
}


# Estimate model on rolling window 
prediction2 <- matrix(NA,T+wndw,10000)
losses2<-matrix(NA,T+wndw,10000)
losses2_sqr<-matrix(NA,T+wndw,10000)
prediction4 <- matrix(NA,T+wndw,10000)
losses4<-matrix(NA,T+wndw,10000)
losses4_sqr<-matrix(NA,T+wndw,10000)
for (ii in 1:10000) {
  for (kk in wndw:(T+wndw-1)) {
    # estimate with data from (kk-wndw+1) : kk
    fit<-lm(y2_sample[(kk-wndw+1): kk,ii]~dplyr::lead(y2_sample[(kk-wndw+1): kk,ii])+dplyr::lead(y2_sample[(kk-wndw+1): kk,ii],2)+dplyr::lead(x_sample[(kk-wndw+1): kk,ii]))
    prediction2[kk+1,ii]<-fit$coefficients[1]+fit$coefficients[2]*y2_sample[kk,ii]+fit$coefficients[3]*y2_sample[kk-1,ii]+fit$coefficients[4]*x_sample[kk,ii]
    losses2[kk+1,ii] <- y2_sample[kk+1,ii] - prediction2[kk+1,ii]
    losses2_sqr[kk+1,ii] <- (y2_sample[kk+1,ii] - prediction2[kk+1,ii])^2
    fit2<-lm(y2_sample[(kk-wndw+1): kk,ii]~dplyr::lead(y2_sample[(kk-wndw+1): kk,ii])+dplyr::lead(y2_sample[(kk-wndw+1): kk,ii],2))
    prediction4[kk+1,ii]<-fit2$coefficients[1]+fit2$coefficients[2]*y2_sample[kk,ii]+fit$coefficients[3]*y2_sample[kk-1,ii]
    losses4[kk+1,ii] <- y2_sample[kk+1,ii] - prediction4[kk+1,ii]
    losses4_sqr[kk+1,ii] <- (y2_sample[kk+1,ii] - prediction4[kk+1,ii])^2
  }
}





gw_stats1<-matrix(NA,T+wndw,10000)
gw_stats1<-(losses4_sqr-losses2_sqr)*100

count<-matrix(NA,1,10000)

for (ii in 1:10000) {
  gw_stats1_eq<-lm(gw_stats1[21:140,ii]~1)
  if (summary(gw_stats1_eq)$coefficients[1, 4]<0.10){
    count[,ii]=1
  }else{
    count[,ii]=0  
  }
}

#This code calculates the analysis percentage of rejection rates (at 10% sign. level)
#over 10000 simulated samples. The test is based on Giacomini and White (2006) test.
GW_test_prctg_change<-100*(rowSums(count)/10000)





