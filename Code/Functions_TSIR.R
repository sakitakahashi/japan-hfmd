#################################################
## Functions for main 1-serotype TSIR analysis ##
#################################################

Est_UR_1sero_linear <- function(data, which_serotype, which_proportion, start_year_infer, max_year) {
	
	## Subset the data to be only after start_year_infer and before max_year (inclusive)
	data <- data[data$Year >= start_year_infer & data$Year <= max_year, ]
	
	## Which serotype and which type of proportion?
	Ct <- data[,paste0("p_", which_serotype, "_", which_proportion)] * data$HFMD_per_sentinel * data$Sentinels
	
	## Get other pertinent data
	Bt <- data$Births_interpolated
	Nt <- data$Population_interpolated
	Nbar <- mean(Nt)
	week_of_year <- data$Week_of_year
	year <- data$Year
	weeks_per_year <- as.numeric(tapply(week_of_year, year, max))
	total_weeks <- length(week_of_year)
	total_weeks_minus_1 <- total_weeks-1
	min_year_minus_1 <- start_year_infer-1
	max_year_plus_1 <- max(year)+1
	max_week <- max(weeks_per_year[1:(length(weeks_per_year)-1)])
	max_week_plus_1 <- max_week+1
	
	## Get the under-reporting rate and deviances about Sbar
	cum_reg <- lm(cumsum(Bt) ~ cumsum(Ct))
	ur <- 1/cum_reg$coef[2]
	Dt <- residuals(cum_reg)
	
	## Correct time series for under-reporting
	It <- Ct/ur
	
	## More parameters for susceptible reconstruction
	which_old <- 1:total_weeks_minus_1
	which_new <- 2:total_weeks
	logIt_new <- log(It[which_new]); logIt_new[which(logIt_new==-Inf)] <- NA
	logIt_old <- log(It[which_old]); logIt_old[which(logIt_old==-Inf)] <- NA
	logNt_old <- log(Nt[which_old])
	Nt_old <- Nt[which_old]
	Dt_old <- Dt[which_old]
	seas_old <- week_of_year[which_old]
	year_old <- year[which_old]
	
	return(list(Ct=Ct, Bt=Bt, Nt=Nt, Nbar=Nbar, week_of_year=week_of_year, year=year, weeks_per_year=weeks_per_year, total_weeks=total_weeks, total_weeks_minus_1=total_weeks_minus_1, min_year_minus_1=min_year_minus_1, max_year_plus_1=max_year_plus_1, max_week=max_week, max_week_plus_1=max_week_plus_1, ur=ur, Dt=Dt, It=It, which_old=which_old, which_new=which_new, logIt_new=logIt_new, logIt_old=logIt_old, logNt_old=logNt_old, Dt_old=Dt_old, seas_old=seas_old, year_old=year_old, Nt_old=Nt_old))
	
}

Sim_TSIR_1sero <- function(startS, startI, alpha, beta, Bt, Nt, seas, deterministic=TRUE) {
	
	## Populate with initial conditions
	St <- It <- rep(NA, length(seas))
	St[1] <- startS
	It[1] <- startI
	
	for(i in 1:(length(seas)-1)) {
		
		## Infecteds
		muval <- beta[seas[i]] * (It[i]^alpha) * St[i] / Nt[i]
		if(deterministic) It[i+1] <- muval
		else It[i+1] <- rpois(1, max(muval,1))

		## Get over transients at start
		if(It[i+1]==0) It[i+1] <- 1
		
		## Susceptibles
		St[i+1] <- max(St[i] + Bt[i] - It[i+1], 0)
	
	}
	
	## Get back the number of susceptible and infected individuals
	return(list(St=St, It=It))
	
}

Run_TSIR_1sero <- function(data, start_year_infer, max_year, alpha, color=c("darkseagreen","firebrick"), which_plot_return=c("beta","fwdsim")) {
	
	##################
	## To plot beta ##
	##################
	
	## Loop over Sbar values
	Sbar <- seq(0.01, 0.3, by=0.001) * data$Nbar
	length_Sbar <- length(Sbar)
	loglik <- rep(NA,length_Sbar)
	
	for(i in 1:length_Sbar) {
	
		logSt_old_tmp <- log(Sbar[i] + data$Dt_old)
		logSt_old_tmp[which(logSt_old_tmp==-Inf)] <- NA
		glmfit_tmp <- lm(data$logIt_new ~ -1 + as.factor(data$seas_old) + offset(alpha * data$logIt_old + logSt_old_tmp - data$logNt_old))
		loglik[i] <- deviance(glmfit_tmp)
		rm(logSt_old_tmp)
		
	}
	
	## Get Sbar_mle
	Sbar_mle <- Sbar[which(loglik==min(loglik))]
	
	## Re-fit, now using Sbar_mle
	logSt_old <- log(Sbar_mle + data$Dt_old)
	St <- Sbar_mle + data$Dt
	glmfit <- lm(data$logIt_new ~ -1 + as.factor(data$seas_old) + offset(alpha * data$logIt_old + logSt_old - data$logNt_old))
	
	## Estimate of the beta_week
	beta_week <- exp(glmfit$coef[1:data$max_week])
	beta_week_CI <- exp(confint(glmfit))[1:data$max_week,]
	beta_to_plot <- data.frame(week=1:data$max_week, beta_mean=beta_week, beta_lb=beta_week_CI[,1], beta_ub=beta_week_CI[,2])
	
	## Plot the beta_week
	beta_to_plot %>%
		ggplot(aes(x=week, y=beta_mean)) +
			geom_errorbar(aes(ymin=beta_lb, ymax=beta_ub)) +
			theme_bw() +
			xlab("Week of year") +
			ylab(expression(beta[s])) +
			geom_line(color=color) +
			geom_point(color=color) +
			theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm")) -> p_beta
	
	#################
	## To fit TSIR ##
	#################
	
	## Fit the joint regression
	glmfit <- lm(data$logIt_new ~ -1 + as.factor(data$seas_old) + data$Dt_old + offset(alpha * data$logIt_old))
	
	## Estimate of sbar
	sbar <- (1/glmfit$coef["data$Dt_old"]) / data$Nbar
	St <- (1/glmfit$coef["data$Dt_old"]) + data$Dt
	
	## Estimate of the beta_week
	log_beta_week <- glmfit$coef[1:data$max_week]
	beta_week <- exp(log_beta_week)
	
	## Run a single deterministic simulation (beta = beta_week/sbar in this formulation)
	out_det <- Sim_TSIR_1sero(startS=St[1], startI=data$It[1], alpha=alpha, beta=beta_week/sbar, Bt=data$Bt, Nt=data$Nt, seas=data$week_of_year, deterministic=TRUE)
	
	## Down-sample
	out_det$It_samp <- round(out_det$It * data$ur)
	
	## Prepare to plot forward sim
	divisor <- tapply(data$week_of_year, data$year, max)
	divisor_rep <- as.numeric(rep(divisor, divisor))
	date_decimal <- data$year + ((data$week_of_year-1) / divisor_rep)
	
	## Get the data
	fwdsim_to_plot <- data.frame(time=date_decimal, obs=data$Ct, pred=-out_det$It_samp)
	y_range <- max(fwdsim_to_plot$obs, abs(fwdsim_to_plot$pred))
	
	## Change to long format
	gather(fwdsim_to_plot, type, value, -time) -> fwdsim_to_plot
	
	## Plot the forward sim
	fwdsim_to_plot %>%
		ggplot(aes(x=time, y=value, colour=type)) +
			theme_bw() +
			xlab("Time") +
			ylab("Cases") +
			geom_line(size=1) +
			scale_x_continuous(breaks=seq(start_year_infer,max_year+1,by=1), minor_breaks=seq(start_year_infer,max_year+1,by=1)) +
			scale_y_continuous(limits=c(-y_range,y_range), labels=abs) +
			scale_colour_manual(values=c("black",color), labels=c("Observed", "Predicted")) +
			theme(
				plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
				legend.title=element_blank(),
				legend.position=c(0.925,0.115),
				legend.background=element_rect(linetype="solid", colour="black")) -> p_fwdsim
	
	if(which_plot_return=="beta") return(p_beta)
	if(which_plot_return=="fwdsim") return(p_fwdsim)
	
}

#################################################
## Functions for main 2-serotype TSIR analysis ##
#################################################

Est_UR_2sero_linear <- function(data, which_proportion, min_year, start_year_infer, max_year, k_EVA71, k_CVA16) {
	
	## Subset the data to be only after min_year and before max_year (inclusive)
	data <- data[data$Year >= min_year & data$Year <= max_year, ]
	
	## Which type of proportion?
	Ct_EVA71 <- data[,paste0("p_", "EVA71", "_", which_proportion)] * data$HFMD_per_sentinel * data$Sentinels
	Ct_CVA16 <- data[,paste0("p_", "CVA16", "_", which_proportion)] * data$HFMD_per_sentinel * data$Sentinels
	
	## Get other pertinent data
	Bt <- data$Births_interpolated
	Nt <- data$Population_interpolated
	week_of_year <- data$Week_of_year
	year <- data$Year
	
	## Get the whole data set (i.e., including X years for fitting C-P)
	DAT_WHOLE <- data.frame(YEAR=year, WEEK=week_of_year, N.EVA71=Ct_EVA71, N.CVA16=Ct_CVA16, BIRTHS=Bt, POP=Nt)
	rownames(DAT_WHOLE) <- 1:nrow(DAT_WHOLE)
	length_whole <- nrow(DAT_WHOLE)
	
	## Make a subsetted data set for inference (i.e., taking out X years for fitting C-P)
	DAT_INFER <- DAT_WHOLE[which(DAT_WHOLE$YEAR>=start_year_infer),]
	rownames(DAT_INFER) <- 1:nrow(DAT_INFER)
	weeks_per_year <- as.numeric(tapply(DAT_INFER$WEEK, DAT_INFER$YEAR, max))
	total_weeks <- nrow(DAT_INFER)
	total_weeks_minus_1 <- total_weeks-1
	max_week <- max(weeks_per_year[1:(length(weeks_per_year)-1)])
	max_week_plus_1 <- max_week+1
	length_infer <- nrow(DAT_INFER)
	
	## Make vectors of k_EVA71 and k_CVA16 to loop over
	start_week_infer <- min(which(DAT_WHOLE$YEAR==start_year_infer))
	start_week_infer_minus_1 <- start_week_infer-1
	
	## Get the under-reporting rate
	Y.EVA71 <- Y.CVA16 <- cumsum(DAT_INFER$BIRTHS)
	
	X.EVA71 <- cumsum(DAT_INFER$N.EVA71) 
	X.CVA16 <- cumsum(DAT_INFER$N.CVA16)
	
	CP_EVA71 <- cumsum(DAT_INFER$N.EVA71 - DAT_WHOLE$N.EVA71[(start_week_infer-k_EVA71):(length_whole-k_EVA71)])	
	CP_CVA16 <- cumsum(DAT_INFER$N.CVA16 - DAT_WHOLE$N.CVA16[(start_week_infer-k_CVA16):(length_whole-k_CVA16)])
	
	cum.reg.EVA71 <- lm(Y.EVA71 ~ X.EVA71)
	rho1.start <- cum.reg.EVA71$coef[2]
	
	cum.reg.CVA16 <- lm(Y.CVA16 ~ X.CVA16)
	rho2.start <- cum.reg.CVA16$coef[2]
	
	## Iterate until estimates of under-reporting convergence 
	for(l in 1:100) {
		
		CP.EVA71 <- (1 * rho2.start * CP_CVA16)
		CP.CVA16 <- (1 * rho1.start * CP_EVA71)
		
		cum.reg.EVA71.update <- lm(Y.EVA71 ~ X.EVA71 + offset(CP.EVA71))
		rho1.update <- cum.reg.EVA71.update$coef[2]
		
		cum.reg.CVA16.update <- lm(Y.CVA16 ~ X.CVA16 + offset(CP.CVA16))
		rho2.update <- cum.reg.CVA16.update$coef[2]
		
		rho1.start <- rho1.update
		rho2.start <- rho2.update
		
		l <- l + 1
		
		rm(CP.EVA71, CP.CVA16, rho1.update, rho2.update)
		
	}
	
	## Get the deviances about Sbar
	DEV.EVA71 <- residuals(cum.reg.EVA71.update)
	DEV.CVA16 <- residuals(cum.reg.CVA16.update)
	
	## Correct time series for under-reporting
	I_EVA71 <- DAT_INFER$N.EVA71 * rho1.start
	I_CVA16 <- DAT_INFER$N.CVA16 * rho2.start
	
	## More parameters for susceptible reconstruction
	which_infer_old <- 1:(length_infer-1)
	which_infer_new <- 2:length_infer
	
	log.I.new.EVA71 <- log(I_EVA71)[which_infer_new]
	log.I.old.EVA71 <- log(I_EVA71)[which_infer_old]
	D.old.EVA71 <- DEV.EVA71[which_infer_old]
	
	log.I.new.CVA16 <- log(I_CVA16)[which_infer_new]
	log.I.old.CVA16 <- log(I_CVA16)[which_infer_old]
	D.old.CVA16 <- DEV.CVA16[which_infer_old]
	
	log.I.new.EVA71[which(log.I.new.EVA71==-Inf)] <- NA
	log.I.new.CVA16[which(log.I.new.CVA16==-Inf)] <- NA
	
	log.I.old.EVA71[which(log.I.old.EVA71==-Inf)] <- NA
	log.I.old.CVA16[which(log.I.old.CVA16==-Inf)] <- NA
	
	seas_old <- DAT_INFER$WEEK[which_infer_old]
	N.old <- DAT_INFER$POP[which_infer_old]
	log.N.old.EVA71 <- log.N.old.CVA16 <- log(DAT_INFER$POP)[which_infer_old]
	mean.N.old.EVA71 <- mean.N.old.CVA16 <- mean(DAT_INFER$POP[which_infer_old])
	Nbar <- mean(DAT_INFER$POP)
	
	return(list(DAT_WHOLE=DAT_WHOLE, DAT_INFER=DAT_INFER, rho1.start=rho1.start, rho2.start=rho2.start, DEV.EVA71=DEV.EVA71, DEV.CVA16=DEV.CVA16, D.old.EVA71=D.old.EVA71, D.old.CVA16=D.old.CVA16, log.I.new.EVA71=log.I.new.EVA71, log.I.old.EVA71=log.I.old.EVA71, log.I.new.CVA16=log.I.new.CVA16, log.I.old.CVA16=log.I.old.CVA16, seas_old=seas_old, log.N.old.EVA71=log.N.old.EVA71, log.N.old.CVA16=log.N.old.CVA16, mean.N.old.EVA71=mean.N.old.EVA71, mean.N.old.CVA16=mean.N.old.CVA16, start_week_infer=start_week_infer, start_week_infer_minus_1=start_week_infer_minus_1, length_whole=length_whole, length_infer=length_infer, I_EVA71=I_EVA71, I_CVA16=I_CVA16, which_infer_new=which_infer_new, N.old=N.old, Nbar=Nbar, length_old=length_infer-1))

}

Sim_TSIR_2sero <- function(startS_EVA71, startS_CVA16, startI_EVA71_with_offset, startI_CVA16_with_offset, alpha_EVA71, alpha_CVA16, beta_EVA71, beta_CVA16, delta_EVA71, delta_CVA16, k_EVA71, k_CVA16, Bt, Nt, seas, start_week_infer, start_week_infer_minus_1, length_infer, deterministic=TRUE) {
	
	## Populate with initial conditions
	Sval_EVA71 <- Sval_CVA16 <- rep(NA, length_infer)
	Ival_EVA71 <- Ival_CVA16 <- rep(NA, length_infer + start_week_infer_minus_1)
	Sval_EVA71[1] <- startS_EVA71
	Sval_CVA16[1] <- startS_CVA16
	Ival_EVA71[1:start_week_infer] <- startI_EVA71_with_offset[1:start_week_infer]
	Ival_CVA16[1:start_week_infer] <- startI_CVA16_with_offset[1:start_week_infer]
	
	for(i in 1:(length_infer-1)) {
		
		## Infecteds: EVA71
		muval.EVA71 <- (beta_EVA71[seas[i]] * (Ival_EVA71[start_week_infer_minus_1+i]^alpha_EVA71) * Sval_EVA71[i])/Nt[i]
		if(deterministic) Ival_EVA71[start_week_infer_minus_1+i+1] <- muval.EVA71
		else Ival_EVA71[start_week_infer_minus_1+i+1] <- rpois(1, muval.EVA71)
		
		## Infecteds: CVA16
		muval.CVA16 <- (beta_CVA16[seas[i]] * (Ival_CVA16[start_week_infer_minus_1+i]^alpha_CVA16) * Sval_CVA16[i])/Nt[i]
		if(deterministic) Ival_CVA16[start_week_infer_minus_1+i+1] <- muval.CVA16
		else Ival_CVA16[start_week_infer_minus_1+i+1] <- rpois(1, muval.CVA16)
		
		## Get over transients at start
		if (!is.na(Ival_EVA71[start_week_infer_minus_1+i+1]) && Ival_EVA71[start_week_infer_minus_1+i+1]==0) {Ival_EVA71[start_week_infer_minus_1+i+1] <- 1}
		if (!is.na(Ival_CVA16[start_week_infer_minus_1+i+1]) && Ival_CVA16[start_week_infer_minus_1+i+1]==0) {Ival_CVA16[start_week_infer_minus_1+i+1] <- 1}
		
		## Incorporate CP
		CP_EVA71 <- delta_CVA16 * (Ival_CVA16[start_week_infer_minus_1+i] - Ival_CVA16[start_week_infer_minus_1+i-k_CVA16])
		CP_CVA16 <- delta_EVA71 * (Ival_EVA71[start_week_infer_minus_1+i] - Ival_EVA71[start_week_infer_minus_1+i-k_EVA71])

		## Susceptibles: incorporating cross-protection
		Sval_EVA71[i+1] <- max(Sval_EVA71[i] + Bt[i] - Ival_EVA71[start_week_infer_minus_1+i] - CP_EVA71, 0)
		Sval_CVA16[i+1] <- max(Sval_CVA16[i] + Bt[i] - Ival_CVA16[start_week_infer_minus_1+i] - CP_CVA16, 0)
		
		rm(CP_EVA71, CP_CVA16)
		
	}
	
	## Get back the number of infected & susceptible individuals
	return(list(Sval_EVA71=Sval_EVA71,
				Sval_CVA16=Sval_CVA16,
				Ival_EVA71=Ival_EVA71[start_week_infer:length(Ival_EVA71)],
				Ival_CVA16=Ival_CVA16[start_week_infer:length(Ival_CVA16)]))
	
}

Sim_TSIR_2sero_reset <- function(startS_EVA71, startS_CVA16, startI_EVA71_with_offset, startI_CVA16_with_offset, alpha_EVA71, alpha_CVA16, beta_EVA71, beta_CVA16, delta_EVA71, delta_CVA16, k_EVA71, k_CVA16, Bt, Nt, seas, start_week_infer, start_week_infer_minus_1, length_infer, deterministic=TRUE, reset_week, reset_It_EVA71_with_lag, reset_It_CVA16_with_lag, reset_St_EVA71, reset_St_CVA16) {
	
	## Populate with initial conditions
	Sval_EVA71 <- Sval_CVA16 <- rep(NA, length_infer)
	Ival_EVA71 <- Ival_CVA16 <- rep(NA, length_infer + start_week_infer_minus_1)
	Sval_EVA71[1] <- startS_EVA71
	Sval_CVA16[1] <- startS_CVA16
	Ival_EVA71[1:start_week_infer] <- startI_EVA71_with_offset[1:start_week_infer]
	Ival_CVA16[1:start_week_infer] <- startI_CVA16_with_offset[1:start_week_infer]
	
	for (i in 1:(length_infer-1)) {
		
		## Infecteds: EVA71
		muval.EVA71 <- (beta_EVA71[seas[i]] * (Ival_EVA71[start_week_infer_minus_1+i]^alpha_EVA71) * Sval_EVA71[i])/Nt[i]
		if(deterministic) Ival_EVA71[start_week_infer_minus_1+i+1] <- muval.EVA71
		else Ival_EVA71[start_week_infer_minus_1+i+1] <- rpois(1, muval.EVA71)
		
		## Infecteds: CVA16
		muval.CVA16 <- (beta_CVA16[seas[i]] * (Ival_CVA16[start_week_infer_minus_1+i]^alpha_CVA16) * Sval_CVA16[i])/Nt[i]
		if(deterministic) Ival_CVA16[start_week_infer_minus_1+i+1] <- muval.CVA16
		else Ival_CVA16[start_week_infer_minus_1+i+1] <- rpois(1, muval.CVA16)
		
		## Get over transients at start
		if (!is.na(Ival_EVA71[start_week_infer_minus_1+i+1]) && Ival_EVA71[start_week_infer_minus_1+i+1]==0) {Ival_EVA71[start_week_infer_minus_1+i+1] <- 1}
		if (!is.na(Ival_CVA16[start_week_infer_minus_1+i+1]) && Ival_CVA16[start_week_infer_minus_1+i+1]==0) {Ival_CVA16[start_week_infer_minus_1+i+1] <- 1}
		
		## Incorporate CP
		CP_EVA71 <- delta_CVA16 * (Ival_CVA16[start_week_infer_minus_1+i] - Ival_CVA16[start_week_infer_minus_1+i-k_CVA16])
		CP_CVA16 <- delta_EVA71 * (Ival_EVA71[start_week_infer_minus_1+i] - Ival_EVA71[start_week_infer_minus_1+i-k_EVA71])
		
		## Susceptibles: incorporating cross-protection
		Sval_EVA71[i+1] <- max(Sval_EVA71[i] + Bt[i] - Ival_EVA71[start_week_infer_minus_1+i] - CP_EVA71, 0) 
		Sval_CVA16[i+1] <- max(Sval_CVA16[i] + Bt[i] - Ival_CVA16[start_week_infer_minus_1+i] - CP_CVA16, 0) 
		
		rm(CP_EVA71, CP_CVA16)
		
		if(i==reset_week) {
			
			Ival_EVA71[(start_week_infer_minus_1+i+1-k_CVA16):(start_week_infer_minus_1+i+1)] <- reset_It_EVA71_with_lag
			Ival_CVA16[(start_week_infer_minus_1+i+1-k_EVA71):(start_week_infer_minus_1+i+1)] <- reset_It_CVA16_with_lag
			Sval_EVA71[i+1] <- reset_St_EVA71
			Sval_CVA16[i+1] <- reset_St_CVA16
			
		}
		
	}
	
	## Get back the number of infected & susceptible individuals
	return(list(Sval_EVA71=Sval_EVA71,
				Sval_CVA16=Sval_CVA16,
				Ival_EVA71=Ival_EVA71[start_week_infer:length(Ival_EVA71)],
				Ival_CVA16=Ival_CVA16[start_week_infer:length(Ival_CVA16)]))
	
}

Run_TSIR_2sero <- function(data, data_train, data_test, start_year_infer, max_year, k_EVA71, k_CVA16, alpha) {
	
	##################
	## To plot beta ##
	##################
	
	## Loop over sbar values
	sbar.EVA71 <- seq(0.02, 0.5, by=0.001)
	sbar.CVA16 <- seq(0.02, 0.5, by=0.001)
	length_sbar <- length(sbar.EVA71)
	loglik.EVA71 <- loglik.CVA16 <- rep(NA,length_sbar)
	
	for(i in 1:length(sbar.EVA71)){
		
		log.S.old <- log(sbar.EVA71[i] * data$mean.N.old.EVA71 + data$D.old.EVA71)
		log.S.old[which(log.S.old==-Inf)] <- NA
		glmfit_tmp <- lm(data$log.I.new.EVA71 ~ -1 + as.factor(data$seas_old) + offset(alpha * data$log.I.old.EVA71 + log.S.old - data$log.N.old.EVA71))
		loglik.EVA71[i] <- deviance(glmfit_tmp)
		rm(log.S.old)
		
	}
	
	for(i in 1:length(sbar.CVA16)){
		
		log.S.old <- log(sbar.CVA16[i] * data$mean.N.old.CVA16 + data$D.old.CVA16)
		log.S.old[which(log.S.old==-Inf)] <- NA
		glmfit_tmp <- lm(data$log.I.new.CVA16 ~ -1 + as.factor(data$seas_old) + offset(alpha * data$log.I.old.CVA16 + log.S.old - data$log.N.old.CVA16))
		loglik.CVA16[i] <- deviance(glmfit_tmp)
		rm(log.S.old)
		
	}
	
	# Get sbar_mle for EVA71 and re-fit
	sbar_mle.EVA71 <- sbar.EVA71[which(loglik.EVA71==min(loglik.EVA71))] * data$mean.N.old.EVA71
	log.S.old.EVA71 <- log(sbar_mle.EVA71 + data$D.old.EVA71)
	glmfit.EVA71 <- lm(data$log.I.new.EVA71 ~ -1 + as.factor(data$seas_old) + offset(alpha * data$log.I.old.EVA71 + log.S.old.EVA71 - data$log.N.old.EVA71))
	alpha_EVA71 <- alpha
	beta_EVA71_plot <- exp(glmfit.EVA71$coef[1:53])
	beta_EVA71_plot_CI <- exp(confint(glmfit.EVA71))[1:53,]
	
	# Get sbar_mle for CVA16 and re-fit
	sbar_mle.CVA16 <- sbar.CVA16[which(loglik.CVA16==min(loglik.CVA16))] * data$mean.N.old.CVA16
	log.S.old.CVA16 <- log(sbar_mle.CVA16 + data$D.old.CVA16)
	glmfit.CVA16 <- lm(data$log.I.new.CVA16 ~ -1 + as.factor(data$seas_old) + offset(alpha * data$log.I.old.CVA16 + log.S.old.CVA16 - data$log.N.old.CVA16))
	alpha_CVA16 <- alpha
	beta_CVA16_plot <- exp(glmfit.CVA16$coef[1:53])
	beta_CVA16_plot_CI <- exp(confint(glmfit.CVA16))[1:53,]
	
	## Plot the beta_week for EVA71
	beta_EVA71 <- data.frame(week=1:53, beta_mean=beta_EVA71_plot, beta_lb=beta_EVA71_plot_CI[,1], beta_ub=beta_EVA71_plot_CI[,2])
	ggplot(beta_EVA71, aes(x=week, y=beta_mean)) +
		geom_errorbar(aes(ymin=beta_lb, ymax=beta_ub)) +
		theme_bw() +
		xlab("Week of year") +
		ylab(expression(beta[s])) +
		geom_line(colour="darkseagreen") +
		geom_point(colour="darkseagreen") +
		theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm")) -> p1
	
	## Plot the beta_week for CVA16
	beta_CVA16 <- data.frame(week=1:53, beta_mean=beta_CVA16_plot, beta_lb=beta_CVA16_plot_CI[,1], beta_ub=beta_CVA16_plot_CI[,2])
	ggplot(beta_CVA16, aes(x=week, y=beta_mean)) +
		geom_errorbar(aes(ymin=beta_lb, ymax=beta_ub)) +
		theme_bw() +
		xlab("Week of year") +
		ylab(expression(beta[s])) +
		geom_line(colour="firebrick") +
		geom_point(colour="firebrick") +
		theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm")) -> p4
	
	#################
	## To fit TSIR ##
	#################
	
	## Fit the joint regression for EVA71
	glmfit.EVA71 <- lm(data$log.I.new.EVA71 ~ -1 + as.factor(data$seas_old) + data$D.old.EVA71 + offset(alpha * data$log.I.old.EVA71))
	alpha_EVA71 <- alpha
	sbar_EVA71 <- (1/glmfit.EVA71$coef["data$D.old.EVA71"]) / data$Nbar
	S_EVA71 <- (1/glmfit.EVA71$coef["data$D.old.EVA71"]) + data$D.old.EVA71
	log_beta_EVA71 <- glmfit.EVA71$coef[1:53]
	beta_EVA71 <- exp(log_beta_EVA71)
	
	## Fit the joint regression for CVA16
	glmfit.CVA16 <- lm(data$log.I.new.CVA16 ~ -1 + as.factor(data$seas_old) + data$D.old.CVA16 + offset(alpha * data$log.I.old.CVA16))
	alpha_CVA16 <- alpha
	sbar_CVA16 <- (1/glmfit.CVA16$coef["data$D.old.CVA16"]) / data$Nbar
	S_CVA16 <- (1/glmfit.CVA16$coef["data$D.old.CVA16"]) + data$D.old.CVA16
	log_beta_CVA16 <- glmfit.CVA16$coef[1:53]
	beta_CVA16 <- exp(log_beta_CVA16)
	
	## Run a single deterministic simulation (beta = beta_week/sbar in this formulation)
	out_det <- Sim_TSIR_2sero(startS_EVA71=S_EVA71[1], startS_CVA16=S_CVA16[1], startI_EVA71_with_offset=data$DAT_WHOLE$N.EVA71[1:data$start_week_infer] * data$rho1.start, startI_CVA16_with_offset=data$DAT_WHOLE$N.CVA16[1:data$start_week_infer] * data$rho2.start, alpha_EVA71=alpha_EVA71, alpha_CVA16=alpha_CVA16, beta_EVA71=beta_EVA71/sbar_EVA71, beta_CVA16=beta_CVA16/sbar_CVA16, delta_EVA71=1, delta_CVA16=1, k_EVA71=k_EVA71, k_CVA16=k_CVA16, Bt=data$DAT_INFER$BIRTHS, Nt=data$DAT_INFER$POP, seas=data$DAT_INFER$WEEK, start_week_infer=data$start_week_infer, start_week_infer_minus_1=data$start_week_infer_minus_1, length_infer=data$length_infer, deterministic=TRUE)
	
	## Down-sample
	out_det$Ival_EVA71_samp <- round(out_det$Ival_EVA71 * 1/data$rho1.start)
	out_det$Ival_CVA16_samp <- round(out_det$Ival_CVA16 * 1/data$rho2.start)
	
	## Prepare to plot forward sims for EVA71
	divisor <- tapply(data$DAT_INFER$WEEK, data$DAT_INFER$YEAR, max)
	divisor_rep <- as.numeric(rep(divisor, divisor))
	date_decimal <- data$DAT_INFER$YEAR + ((data$DAT_INFER$WEEK-1) / divisor_rep)
	fwdsim_EVA71 <- data.frame(time=date_decimal, obs=data$DAT_INFER$N.EVA71, pred=out_det$Ival_EVA71_samp)
	y_range_EVA71 <- max(fwdsim_EVA71$obs, abs(fwdsim_EVA71$pred))
	fwdsim_EVA71_plot <- data.frame(time=date_decimal, obs=data$DAT_INFER$N.EVA71, pred=-out_det$Ival_EVA71_samp)

	## Prepare to plot forward sims for CVA16
	divisor <- tapply(data$DAT_INFER$WEEK, data$DAT_INFER$YEAR, max)
	divisor_rep <- as.numeric(rep(divisor, divisor))
	date_decimal <- data$DAT_INFER$YEAR + ((data$DAT_INFER$WEEK-1) / divisor_rep)
	fwdsim_CVA16 <- data.frame(time=date_decimal, obs=data$DAT_INFER$N.CVA16, pred=out_det$Ival_CVA16_samp)
	y_range_CVA16 <- max(fwdsim_CVA16$obs, abs(fwdsim_CVA16$pred))
	fwdsim_CVA16_plot <- data.frame(time=date_decimal, obs=data$DAT_INFER$N.CVA16, pred=-out_det$Ival_CVA16_samp)
	
	## Change to long format
	gather(fwdsim_EVA71_plot, type, value, -time) -> fwdsim_EVA71_plot
	gather(fwdsim_CVA16_plot, type, value, -time) -> fwdsim_CVA16_plot
	
	## Plot the forward sim for EVA71
	fwdsim_EVA71_plot %>%
		ggplot(aes(x=time, y=value, colour=type)) +
			theme_bw() +
			xlab("Time") +
			ylab("Cases") +
			geom_line(size=1) +
			scale_x_continuous(breaks=seq(start_year_infer,max_year+1,by=1), minor_breaks=seq(start_year_infer,max_year+1,by=1)) +
			scale_y_continuous(limits=c(-y_range_EVA71, y_range_EVA71), labels=abs) +
			scale_colour_manual(values=c("black", "darkseagreen"), labels=c("Observed", "Predicted")) +
			theme(
				plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
				legend.title=element_blank(),
				legend.position=c(0.925,0.115),
				legend.background=element_rect(linetype="solid", colour="black")) -> p2
	
	## Plot the forward sim for CVA16
	fwdsim_CVA16_plot %>%
		ggplot(aes(x=time, y=value, colour=type)) +
			theme_bw() +
			xlab("Time") +
			ylab("Cases") +
			geom_line(size=1) +
			scale_x_continuous(breaks=seq(start_year_infer,max_year+1,by=1), minor_breaks=seq(start_year_infer,max_year+1,by=1)) +
			scale_y_continuous(limits=c(-y_range_CVA16, y_range_CVA16), labels=abs) +
			scale_colour_manual(values=c("black", "firebrick"), labels=c("Observed", "Predicted")) +
			theme(
				plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
				legend.title=element_blank(),
				legend.position=c(0.925,0.115),
				legend.background=element_rect(linetype="solid", colour="black")) -> p5
	
	## Make xyplot for EVA71
	X <- unname(tapply(fwdsim_EVA71$obs, (seq_along(fwdsim_EVA71$obs)-1) %/% 4, sum))
	Y <- unname(tapply(fwdsim_EVA71$pred, (seq_along(fwdsim_EVA71$pred)-1) %/% 4, sum))
	Z <- 1:length(X)
	dat3 <- data.frame(x=X, y=Y, z=Z)
	
	ggplot(dat3, aes(x=x, y=y)) +
		geom_point(size=1.5) +
		theme_bw() +
		xlab("Observed") +
		ylab("Predicted") +
		xlim(c(0, max(dat3$x,dat3$y))) +
		ylim(c(0, max(dat3$x,dat3$y))) +
		geom_segment(x=0, xend=max(dat3$x,dat3$y), y=0, yend=max(dat3$x,dat3$y), colour="black", size=0.3, linetype="dotted") +
		geom_smooth(method="lm", formula="y~-1+x", se=T, colour="darkseagreen", fill="darkseagreen", alpha=0.05, size=1) +
		theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm")) -> p3
	
	## Make xyplot for CVA16
	X <- unname(tapply(fwdsim_CVA16$obs, (seq_along(fwdsim_CVA16$obs)-1) %/% 4, sum))
	Y <- unname(tapply(fwdsim_CVA16$pred, (seq_along(fwdsim_CVA16$pred)-1) %/% 4, sum))
	Z <- 1:length(X)
	dat6 <- data.frame(x=X, y=Y, z=Z)
	
	ggplot(dat6, aes(x=x, y=y)) +
		geom_point(size=1.5) +
		theme_bw() +
		xlab("Observed") +
		ylab("Predicted") +
		xlim(c(0, max(dat6$x,dat6$y))) +
		ylim(c(0, max(dat6$x,dat6$y))) +
		geom_segment(x=0, xend=max(dat6$x,dat6$y), y=0, yend=max(dat6$x,dat6$y), colour="black", size=0.3, linetype="dotted") +
		geom_smooth(method="lm", formula="y~-1+x", se=T, colour="firebrick", fill="firebrick", alpha=0.05, size=1) +
		theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm")) -> p6
	
	##########################
	## To fit out-of-sample ##
	##########################
	
	## Fit the joint regression for EVA71
	glmfit.EVA71 <- lm(data_train$log.I.new.EVA71 ~ -1 + as.factor(data_train$seas_old) + data_train$D.old.EVA71 + offset(alpha * data_train$log.I.old.EVA71))
	alpha_EVA71 <- alpha
	sbar_EVA71 <- (1/glmfit.EVA71$coef["data_train$D.old.EVA71"]) / data_train$Nbar
	S_EVA71 <- (1/glmfit.EVA71$coef["data_train$D.old.EVA71"]) + data_train$D.old.EVA71
	log_beta_EVA71 <- glmfit.EVA71$coef[1:53]
	beta_EVA71 <- exp(log_beta_EVA71)
	
	## Fit the joint regression for CVA16
	glmfit.CVA16 <- lm(data_train$log.I.new.CVA16 ~ -1 + as.factor(data_train$seas_old) + data_train$D.old.CVA16 + offset(alpha * data_train$log.I.old.CVA16))
	alpha_CVA16 <- alpha
	sbar_CVA16 <- (1/glmfit.CVA16$coef["data_train$D.old.CVA16"]) / data_train$Nbar
	S_CVA16 <- (1/glmfit.CVA16$coef["data_train$D.old.CVA16"]) + data_train$D.old.CVA16
	log_beta_CVA16 <- glmfit.CVA16$coef[1:53]
	beta_CVA16 <- exp(log_beta_CVA16)
	
	# Run 1 deterministic simulation
	reset_week = length(data_train$seas_old)
	out_det_oos <- Sim_TSIR_2sero_reset(startS_EVA71=S_EVA71[1], startS_CVA16=S_CVA16[1], startI_EVA71_with_offset=data_train$DAT_WHOLE$N.EVA71[1:data_train$start_week_infer] * data_train$rho1.start, startI_CVA16_with_offset=data_train$DAT_WHOLE$N.CVA16[1:data_train$start_week_infer] * data_train$rho2.start, alpha_EVA71=alpha_EVA71, alpha_CVA16=alpha_CVA16, beta_EVA71=beta_EVA71/sbar_EVA71, beta_CVA16=beta_CVA16/sbar_CVA16, delta_EVA71=1, delta_CVA16=1, k_EVA71=k_EVA71, k_CVA16=k_CVA16, Bt=data_test$DAT_INFER$BIRTHS, Nt=data_test$DAT_INFER$POP, seas=data_test$DAT_INFER$WEEK, start_week_infer=data_test$start_week_infer, start_week_infer_minus_1=data_test$start_week_infer_minus_1, length_infer=data_test$length_infer, deterministic=TRUE, reset_week=reset_week, reset_It_EVA71_with_lag=data_train$I_EVA71[(reset_week-k_CVA16):reset_week], reset_It_CVA16_with_lag=data_train$I_CVA16[(reset_week-k_EVA71):reset_week], reset_St_EVA71=S_EVA71[reset_week], reset_St_CVA16=S_CVA16[reset_week])
	
	## Down-sample
	out_det_oos$Ival_EVA71_samp <- round(out_det_oos$Ival_EVA71 * 1/data_train$rho1.start)
	out_det_oos$Ival_CVA16_samp <- round(out_det_oos$Ival_CVA16 * 1/data_train$rho2.start)
	
	## Which time-points to include?
	which_post <- (length(data_train$seas_old)+2):length(data_test$seas_old)
	
	## Get the out-of-sample model fits for EVA71
	extsim_EVA71 <- data.frame(time=date_decimal[which_post], obs=data_test$DAT_INFER$N.EVA71[which_post], pred=out_det_oos$Ival_EVA71_samp[which_post])
	y_range_extsim_EVA71 <- max(extsim_EVA71$obs, extsim_EVA71$pred)
	
	## Plot the forward sim for EVA71
	ggplot(extsim_EVA71, aes(x=time, y=obs)) +
		theme_bw() +
		geom_line(col="black", size=1) +
		xlab("Time") +
		ylab("Cases") +
		geom_line(aes(x=time, y=-pred), col="cornflowerblue", size=1) +
		theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm")) +
		scale_x_continuous(breaks=seq(last_year_included+1, max_year+1, by=1), minor_breaks=seq(last_year_included, max_year+1, by=1)) +
		scale_y_continuous(limits=c(-y_range_extsim_EVA71, y_range_extsim_EVA71), labels=abs) -> p7

	## Get the out-of-sample model fits for CVA16
	extsim_CVA16 <- data.frame(time=date_decimal[which_post], obs=data_test$DAT_INFER$N.CVA16[which_post], pred=out_det_oos$Ival_CVA16_samp[which_post])
	y_range_extsim_CVA16 <- max(extsim_CVA16$obs, extsim_CVA16$pred)
	
	## Plot the forward sim for CVA16
	ggplot(extsim_CVA16, aes(x=time, y=obs)) +
	theme_bw() +
	geom_line(col="black", size=1) +
	xlab("Time") +
	ylab("Cases") +
	geom_line(aes(x=time, y=-pred), col="purple", size=1) +
	theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm")) +
	scale_x_continuous(breaks=seq(last_year_included+1, max_year+1, by=1), minor_breaks=seq(last_year_included, max_year+1, by=1)) +
	scale_y_continuous(limits=c(-y_range_extsim_CVA16, y_range_extsim_CVA16), labels=abs) -> p8
	
	grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, layout_matrix = rbind(c(1,1,2,2,2,2,3,3), c(1,1,2,2,2,2,7,7), c(4,4,5,5,5,5,6,6), c(4,4,5,5,5,5,8,8)))
	
}
