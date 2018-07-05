library(Hmisc)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

source("Code/Functions_TSIR.R")

data_hfmd <- read.csv("Data/data_raw.csv")

###################
## Make Figure 6 ##
###################

Run_TSIR_2sero_log_periodogram <- function(data, which_proportion, k_EVA71, k_CVA16, alpha, min_year, start_year_infer, max_year, deterministic=TRUE) {
	
	## Estimate under-reporting
	out <- Est_UR_2sero_linear(data=data, which_proportion=which_proportion, min_year=min_year, start_year_infer=start_year_infer, max_year=max_year, k_EVA71=k_EVA71, k_CVA16=k_CVA16)
	
	# Postulate a reasonable range of candidate values
	S.mean.EVA71 <- seq(0.02, 0.5, by=0.001)
	S.mean.CVA16 <- seq(0.02, 0.5, by=0.001)
	
	# Log-likelihoods
	log.lik.EVA71 <- log.lik.CVA16 <- rep(NA, length(S.mean.EVA71))
	
	# Loop over all values of S.mean: EVA71
	for (m in 1:length(S.mean.EVA71)) {
		
		log.S.old <- log(S.mean.EVA71[m] * out$mean.N.old.EVA71 + out$D.old.EVA71)
		log.S.old[which(log.S.old==-Inf)] <- NA
		glmfit_tmp <- lm(out$log.I.new.EVA71 ~ -1 + as.factor(out$seas_old) + offset(alpha*out$log.I.old.EVA71 + log.S.old - out$log.N.old.EVA71))
		log.lik.EVA71[m] <- deviance(glmfit_tmp)
		rm(log.S.old)
		
	}
	
	# Loop over all values of S.mean: CVA16
	for (m in 1:length(S.mean.CVA16)) {
		
		log.S.old <- log(S.mean.CVA16[m] * out$mean.N.old.CVA16 + out$D.old.CVA16)
		log.S.old[which(log.S.old==-Inf)] <- NA
		glmfit_tmp <- lm(out$log.I.new.CVA16 ~ -1 + as.factor(out$seas_old) + offset(alpha*out$log.I.old.CVA16 + log.S.old - out$log.N.old.CVA16))
		log.lik.CVA16[m] <- deviance(glmfit_tmp)
		rm(log.S.old)
		
	}
	
	# Get MLE of sm
	sm.EVA71 <- S.mean.EVA71[which(log.lik.EVA71 == min(log.lik.EVA71))] * out$mean.N.old.EVA71
	log.S.old.EVA71 <- log(sm.EVA71 + out$D.old.EVA71)
	glmfit.EVA71 <- lm(out$log.I.new.EVA71 ~ -1 + as.factor(out$seas_old) + offset(alpha*out$log.I.old.EVA71 + log.S.old.EVA71 - out$log.N.old.EVA71))
	alpha_EVA71 <- alpha
	beta_EVA71 <- exp(glmfit.EVA71$coef[1:53])
	beta_EVA71_CI <- exp(confint(glmfit.EVA71))[1:53,]
	
	sm.CVA16 <- S.mean.CVA16[which(log.lik.CVA16 == min(log.lik.CVA16))] * out$mean.N.old.CVA16
	log.S.old.CVA16 <- log(sm.CVA16 + out$D.old.CVA16)
	glmfit.CVA16 <- lm(out$log.I.new.CVA16 ~ -1 + as.factor(out$seas_old) + offset(alpha*out$log.I.old.CVA16 + log.S.old.CVA16 - out$log.N.old.CVA16))
	alpha_CVA16 <- alpha
	beta_CVA16 <- exp(glmfit.CVA16$coef[1:53])
	beta_CVA16_CI <- exp(confint(glmfit.CVA16))[1:53,]
	
	# Reconstructed susceptibles
	S_EVA71 <- sm.EVA71 + out$DEV.EVA71
	S_CVA16 <- sm.CVA16 + out$DEV.CVA16
	
	# Run 1 deterministic simulation
	out_det <- Sim_TSIR_2sero(startS_EVA71=S_EVA71[1], startS_CVA16=S_CVA16[1], startI_EVA71_with_offset=out$DAT_WHOLE$N.EVA71[1:out$start_week_infer] * out$rho1.start, startI_CVA16_with_offset=out$DAT_WHOLE$N.CVA16[1:out$start_week_infer] * out$rho2.start, alpha_EVA71=alpha_EVA71, alpha_CVA16=alpha_CVA16, beta_EVA71=beta_EVA71, beta_CVA16=beta_CVA16, delta_EVA71=1, delta_CVA16=1, k_EVA71=k_EVA71, k_CVA16=k_CVA16, Bt=out$DAT_INFER$BIRTHS, Nt=out$DAT_INFER$POP, seas=out$DAT_INFER$WEEK, start_week_infer=out$start_week_infer, start_week_infer_minus_1=out$start_week_infer_minus_1, length_infer=out$length_infer, deterministic=TRUE)
	
	## Down-sample
	out_det$Ival_EVA71_samp <- round(out_det$Ival_EVA71 * 1/data$rho1.start)
	out_det$Ival_CVA16_samp <- round(out_det$Ival_CVA16 * 1/data$rho2.start)
	
	## Run sims for longer
	XX_years <- 50
	start_year <- 30
	XX_weeks <- XX_years * 52
	Bt_long <- c(out$DAT_INFER$BIRTHS, rep(tail(out$DAT_INFER$BIRTHS,1), times=XX_weeks))
	Nt_long <- c(out$DAT_INFER$POP, rep(tail(out$DAT_INFER$POP,1), times=XX_weeks))
	seas_long <- c(out$DAT_INFER$WEEK, rep(1:52, times=XX_years))
	length_infer_long <- out$length_infer + XX_weeks
	post_transients <- ((start_year*52)+1):XX_weeks
	# time_post_transients <- seq(start_year, XX_years, length=(XX_years-start_year)*52+1)[1:((XX_years-start_year)*52)]
	
	# Run 1 deterministic simulation
	out_det <- Sim_TSIR_2sero(startS_EVA71=S_EVA71[1], startS_CVA16=S_CVA16[1], startI_EVA71_with_offset=out$DAT_WHOLE$N.EVA71[1:out$start_week_infer] * out$rho1.start, startI_CVA16_with_offset=out$DAT_WHOLE$N.CVA16[1:out$start_week_infer] * out$rho2.start, alpha_EVA71=alpha_EVA71, alpha_CVA16=alpha_CVA16, beta_EVA71=beta_EVA71, beta_CVA16=beta_CVA16, delta_EVA71=1, delta_CVA16=1, k_EVA71=k_EVA71, k_CVA16=k_CVA16, Bt=Bt_long, Nt=Nt_long, seas=seas_long, start_week_infer=out$start_week_infer, start_week_infer_minus_1=out$start_week_infer_minus_1, length_infer=length_infer_long, deterministic=TRUE)
	
	## Down-sample
	out_det$Ival_EVA71_samp <- round(out_det$Ival_EVA71 * 1/data$rho1.start)
	out_det$Ival_CVA16_samp <- round(out_det$Ival_CVA16 * 1/data$rho2.start)
	
	# Classic periodogram plot for EV-A71
	spec_EVA71 <- spectrum(log(out_det$Ival_EVA71[post_transients]), plot=FALSE)
	
	# Classic periodogram plot for CV-A16
	spec_CVA16 <- spectrum(log(out_det$Ival_CVA16[post_transients]), plot=FALSE)
	
	## Object to return
	spec_obj <- data.frame(k_EVA71=k_EVA71, period_EVA71=rev(1/spec_EVA71$freq/52), spectralamp_EVA71=rev(spec_EVA71$spec), k_CVA16=k_CVA16, period_CVA16=rev(1/spec_CVA16$freq/52), spectralamp_CVA16=rev(spec_CVA16$spec))
	tmp_EVA71 <- approx(spec_obj$period_EVA71, spec_obj$spectralamp_EVA71, n=5000)
	tmp_CVA16 <- approx(spec_obj$period_CVA16, spec_obj$spectralamp_CVA16, n=5000)
	
	## Shuffle each of the 2 series many times to generate a 'white noise spectrum'
	n_shuffle <- 100000
	
	## Shuffle the EV-A71
	spec_shuffle_save_EVA71 <- matrix(0, nrow=n_shuffle, ncol=length(spec_EVA71$freq))
	for(k in 1:n_shuffle) {
		
		data_shuffle_tmp <- sample(out_det$Ival_EVA71[post_transients])
		spec_shuffle_tmp <- spectrum(log(data_shuffle_tmp), plot=FALSE)
		spec_shuffle_save_EVA71[k,] <- rev(spec_shuffle_tmp$spec)
		
	}
	
	EVA71_x_lb <- rev(1/spec_EVA71$freq/52)
	EVA71_y_lb_0.025 <- apply(spec_shuffle_save_EVA71, 2, quantile, probs=0.025)
	EVA71_y_lb_0.05 <- apply(spec_shuffle_save_EVA71, 2, quantile, probs=0.05)
	EVA71_xy_lb_0.025 <- approx(EVA71_x_lb, EVA71_y_lb_0.025, n=5000)
	EVA71_xy_lb_0.05 <- approx(EVA71_x_lb, EVA71_y_lb_0.05, n=5000)
	
	## Shuffle the CV-A16
	spec_shuffle_save_CVA16 <- matrix(0, nrow=n_shuffle, ncol=length(spec_CVA16$freq))
	for(k in 1:n_shuffle) {
		
		data_shuffle_tmp <- sample(out_det$Ival_CVA16[post_transients])
		spec_shuffle_tmp <- spectrum(log(data_shuffle_tmp), plot=FALSE)
		spec_shuffle_save_CVA16[k,] <- rev(spec_shuffle_tmp$spec)
		
	}
	
	CVA16_x_lb <- rev(1/spec_CVA16$freq/52)
	CVA16_y_lb_0.025 <- apply(spec_shuffle_save_CVA16, 2, quantile, probs=0.025)
	CVA16_y_lb_0.05 <- apply(spec_shuffle_save_CVA16, 2, quantile, probs=0.05)
	CVA16_xy_lb_0.025 <- approx(CVA16_x_lb, CVA16_y_lb_0.025, n=5000)
	CVA16_xy_lb_0.05 <- approx(CVA16_x_lb, CVA16_y_lb_0.05, n=5000)
	
	## Object to return
	spec_obj <- data.frame(k_EVA71=k_EVA71, period_EVA71=tmp_EVA71$x, spectralamp_EVA71=tmp_EVA71$y, k_CVA16=k_CVA16, period_CVA16=tmp_CVA16$x, spectralamp_CVA16=tmp_CVA16$y, shuffle_spectralamp_EVA71_0.025=EVA71_xy_lb_0.025$y, shuffle_spectralamp_EVA71_0.05=EVA71_xy_lb_0.05$y, shuffle_spectralamp_CVA16_0.025=CVA16_xy_lb_0.025$y, shuffle_spectralamp_CVA16_0.05=CVA16_xy_lb_0.05$y)
	
	return(spec_obj)
	
}

## EVA71's bifurcation diagram
k_EVA71_fix <- 8
alpha <- 0.975
start_year <- 1997

## Fit at k_CVA16=0
save_EVA71 <- NULL
out <- Run_TSIR_2sero_log_periodogram(data=data_hfmd, which_proportion="WMA", k_EVA71=k_EVA71_fix, k_CVA16=0, alpha=alpha, min_year=1982, start_year_infer=start_year, max_year=2015, deterministic=TRUE)
save_EVA71 <- out

## Sweep over values of k_CVA16
for(i in 1:52) {
	
	out <- Run_TSIR_2sero_log_periodogram(data=data_hfmd, which_proportion="WMA", k_EVA71=k_EVA71_fix, k_CVA16=i, alpha=alpha, min_year=1982, start_year_infer=start_year, max_year=2015, deterministic=TRUE)
	save_EVA71 <- rbind(save_EVA71, out)
	print(i)
	
}

## Plot
cols <- colorRampPalette(c("white", "cornflowerblue"))(20)
cols_max_value <- max(c(max(save_EVA71$spectralamp_EVA71), max(save_EVA71$spectralamp_CVA16)))
cols_breaks <- c(1, 1e1, 1e2, 1e3, 2e3)
max_per <- 6
point_y <- 0.11
point_size <- 3.5
point_shape <- 18

ggplot(subset(save_EVA71, period_EVA71<=max_per & spectralamp_EVA71 >= shuffle_spectralamp_EVA71_0.025), aes(x=k_CVA16, y=period_EVA71, fill=spectralamp_EVA71)) +
	geom_raster() +
	scale_fill_gradientn(colours=cols, trans="sqrt", limits=c(0, cols_max_value), breaks=cols_breaks) +
	labs(x="Cross-protection after CV-A16 infection (weeks)", y="Period of EV-A71 (years)", fill="Log EV-A71\nSpectral density") +
	theme_bw() +
	scale_x_continuous(breaks=seq(0,52,by=2), labels=seq(0,52,by=2), expand=c(0,0)) +
	scale_y_continuous(breaks=0:max_per, labels=0:max_per, expand=c(0,0)) +
	guides(fill=guide_colorbar(barwidth=1.5, barheight=9)) +
	geom_tile(data=subset(save_EVA71, period_EVA71<=max_per & spectralamp_EVA71 < shuffle_spectralamp_EVA71_0.025), aes(x=k_CVA16, y=period_EVA71), colour="mistyrose") +
	geom_point(aes(x=k_EVA71_fix, y=point_y), colour="darkseagreen", size=point_size, shape=point_shape) +
	geom_point(aes(x=39, y=point_y), colour="firebrick", size=point_size, shape=point_shape) -> p1

ggplot(subset(save_EVA71, period_CVA16<=max_per & spectralamp_CVA16 >= shuffle_spectralamp_CVA16_0.025), aes(x=k_CVA16, y=period_CVA16, fill=spectralamp_CVA16)) +
	geom_raster() +
	scale_fill_gradientn(colours=cols, trans="sqrt", limits=c(0, cols_max_value), breaks=cols_breaks) +
	labs(x="Cross-protection after CV-A16 infection (weeks)", y="Period of CV-A16 (years)", fill="Log CV-A16\nSpectral density") +
	theme_bw() +
	scale_x_continuous(breaks=seq(0,52,by=2), labels=seq(0,52,by=2), expand=c(0,0)) +
	scale_y_continuous(breaks=0:max_per, labels=0:max_per, expand=c(0,0)) +
	guides(fill=guide_colorbar(barwidth=1.5, barheight=9)) +
	geom_tile(data=subset(save_EVA71, period_CVA16<=max_per & spectralamp_CVA16 < shuffle_spectralamp_CVA16_0.025), aes(x=k_CVA16, y=period_CVA16), colour="mistyrose") +
	geom_point(aes(x=k_EVA71_fix, y=point_y), colour="darkseagreen", size=point_size, shape=point_shape) +
	geom_point(aes(x=39, y=point_y), colour="firebrick", size=point_size, shape=point_shape) -> p2

grid.arrange(p1, p2, ncol=1)
