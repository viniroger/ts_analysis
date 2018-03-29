#!/usr/bin/Rscript
# Script to analise time series and find correlations (two by two)
options(warn = 1) # If 'warn' is one, warnings are printed as they occur

# Define libraries
end_libs='~/Rpacks'
suppressPackageStartupMessages(require(timeDate,lib=end_libs))
suppressPackageStartupMessages(require(zoo,lib=end_libs))
suppressPackageStartupMessages(require(urca,lib=end_libs))
suppressPackageStartupMessages(require(randtests,lib=end_libs))
suppressPackageStartupMessages(require(PMCMR,lib=end_libs))
suppressPackageStartupMessages(require(forecast,lib=end_libs))

## FUNCTIONS ##

analysis = function(v_ts,var_name) {
	# Define output text file
	cat(paste('\n====================== Working with',var_name,'data ======================','\n\n'))
	# Define output plots file
	pdf(paste('Plots_',var_name,'.pdf',sep=''))
	
	# PLOT Decompose on randon/seasonal/trend/observed, multiplicative OR additive
	decompose = decompose(v_ts, type = 'additive')
	plot(decompose)
	
	# PLOT histogram
	hist(decompose$random)
	
	# Trend test: Wald-Wolfowitz, Cox-Stuart and Mann-Kendall (randtests package)
	print(runs.test(v_ts))
	print(cox.stuart.test(v_ts))
	print(rank.test(v_ts))
	
	# Seasonality tests: Kruskal-Wallis and Friedman (PMCMR package)
	month = format(as.Date(time(v_ts)), '%m') # create array with months
	year = format(as.Date(time(v_ts)), '%Y') # create array with years
	print(kruskal.test(as.matrix(v_ts) ~ as.integer(month), data))
	print(friedman.test(v_ts,as.factor(year),as.factor(month)))

	# Stationarity tests: ADF, PP and KPSS (urca package)
	print(summary(ur.df(y=v_ts, selectlags = c('AIC'))))
	print(summary(ur.pp(v_ts, type='Z-tau')))
	print(summary(ur.kpss(y=v_ts, type='tau')))

	# PLOT 1st difference
	plot(diff(v_ts), main='1st diff from time serie')
	
	# PLOT Box-Cox transformation (forecast package)
	lambda = BoxCox.lambda(v_ts)
	print(paste('Lambda =',lambda))
	plot(BoxCox(v_ts, lambda), main='Log transformation from time serie')
	
	# PLOT ACF/PACF graphics
	tsdisplay(v_ts, main='Time serie, ACF and PACF')
	tsdisplay(diff(v_ts,1), main='Time serie, ACF and PACF from 1st diff') # 1 is the difference order
	
	# Build a model
	model = auto.arima(v_ts) # automatically chosen parameters
	#model = Arima(v_ts, order=c(1,1,1), seasonal=c(0,1,1), method=c('ML')) # FILL IN THE VALUES ACCORDING TO YOUR ANALYZES
	print(model)
	
	# PLOT residuals
	tsdiag(model)
	
	# Make and PLOT a forecast
	forecast = forecast(model)
	plot(forecast)
	
	# Return data to CCF and ARIMA
	return(decompose$random)
	#q()
}

## MAIN ##

# Filenames
ts_names = c('date','Oil','CO2')
filename_x = paste('data_',ts_names[2],'.csv',sep='')
filename_y = paste('data_',ts_names[3],'.csv',sep='')

# Load data files (changing -99.99 to NA)
var_x = read.csv(filename_x, header = TRUE, as.is = TRUE, na.strings='-99.99')
var_y = read.csv(filename_y, header = TRUE, as.is = TRUE, na.strings='-99.99')

# Select/Merge data with common dates
data = merge(var_x,var_y, by='date')
colnames(data) = ts_names

# Start date and frequency
temp = strsplit(as.character(data$date[1]), split='-')
start_year_hist = as.numeric(temp[[1]][1])
start_month_hist = as.numeric(temp[[1]][2])
freq = 12
# Build ts + interpolate data missing (zoo/timeDate packages)
x_ts = ts(na.approx(data[,2], rule=2), frequency = freq, start = c(start_year_hist, start_month_hist))
y_ts = ts(na.approx(data[,3], rule=2), frequency = freq, start = c(start_year_hist, start_month_hist))

# Split data (reserve last 10% to calculate model accuracy)
end_hist = c(2014,12)
x_ts_obs = window(x_ts, end=end_hist)
y_ts_obs = window(y_ts, end=end_hist)
start_val = c(2015,1)
x_ts_val = window(x_ts, start=start_val)
y_ts_val = window(y_ts, start=start_val)

# Analysis for each time serie
vars_df = as.data.frame(matrix(ncol=2, nrow=length(x_ts_obs)))
vars_df[,1] = analysis(x_ts_obs,ts_names[2])
vars_df[,2] = analysis(y_ts_obs,ts_names[3])

cat(paste('\n====================== Working with both data ======================','\n\n'))

# PLOT cross-correlation and find max lag - ccf(x,y)
pdf(paste('Plots_',ts_names[2],'_',ts_names[3],'.pdf',sep=''))
vars_df = vars_df[complete.cases(vars_df), ] # remove rows with NAs
#ccf = ccf(data[,ts_names[2]],data[,ts_names[3]],50) # ccf with raw data - useless
ccf_obj = ccf(vars_df[,1],vars_df[,2], main = paste("CCF between",ts_names[2],"and",ts_names[3]))
print(ccf_obj)
cor = ccf_obj$acf[,,1]
lag = ccf_obj$lag[,,1]
res = data.frame(cor,lag)
lag_max = res[which.max(res$cor),]$lag

# Calulate lags (needs zoo package)
var_x_zoo = zoo(as.data.frame(x_ts_obs),as.Date(time(x_ts_obs)))
interval = c(-2:2) # Define interval between min and max lags
lags = lag(var_x_zoo,k=interval[1]:interval[length(interval)])
colnames(lags) = paste('Lag',interval[1]:interval[length(interval)],sep='')

# Choose optimal lag length for x based on AIC
# Restrict data so models use same fitting period
icol = 0
lowest_aicc = 99999
for (i in interval) {
	#print(paste('lag =',i))
	icol = icol + 1
	fit = auto.arima(as.data.frame(y_ts_obs), xreg=lags[,icol], d=0)
	#print(fit)
	if (fit$aicc <= lowest_aicc) {
		lowest_aicc = fit$aicc
		icol_best_k = icol
	}
}
#icol_best_k = 1 # Force column number with best k lag
print(paste('Best lag:',colnames(lags)[icol_best_k]))

# Fit a model to original serie, without and with regressors
print('Model without regressor')
fit = auto.arima(y_ts_obs)
print(fit)
print('Model with regressor')
xreg = lags[,icol_best_k]
fit_reg = auto.arima(y_ts_obs, xreg = xreg)
print(fit_reg)

# Box-Ljung test
BT = Box.test(fit$residuals, lag=30, type = 'Ljung-Box', fitdf=2)
print(BT)

# Forecast regressors
xreg_fc = forecast(x_ts_obs)

# Forecast with ARIMA model, without and with regressors, and PLOT
fc = forecast(fit)
plot(fc, main = 'Forecast without regressors')
fc_reg = forecast(fit_reg, xreg=xreg_fc$mean)
plot(fc_reg, main = 'Forecast with regressors')

# Calculate forecast accuracy measures
print(accuracy(fc$mean, y_ts_val))
print(accuracy(fc_reg$mean, y_ts_val))
