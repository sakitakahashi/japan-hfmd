library(Rwave)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(boot)
library(moments)
library(mgcv)

data_hfmd <- read.csv("Data/data_raw.csv")

###################
## Make Figure 1 ##
###################

## For wavelet
no <- 8
nv <- 16
a <- 2^seq(1, no+1-1/nv, by=1/nv)

## For plotting
divisor <- tapply(data_hfmd$Week_of_year, data_hfmd$Year, max)
divisor_rep <- as.numeric(rep(divisor, divisor))
date_decimal <- data_hfmd$Year + ((data_hfmd$Week_of_year-1)/divisor_rep)

par(mfrow=c(6,1), mar=c(0.5,7,0.5,0.75), oma=c(1.5,0.5,1.5,0.5))

plot(x=date_decimal, data_hfmd$EVA71, type="b", xaxt='n', xlab="", ylab="Counts of EV-A71", col="darkseagreen", xaxs="i", xlim=c(1982,2016))
axis(side=1, at=1982:2016, labels=F)
axis(side=3, at=1982:2016, labels=1982:2016)
abline(v=seq(1984,2014,by=2), col="grey80", lty=2)

wfit <- cwt(sqrt(data_hfmd$EVA71), no, nv, plot=F)
wspec <- Mod(wfit)
image(x=date_decimal, wspec, y=a/52, ylim=c(0,5), xaxt="n", xlab="", ylab="Period of EV-A71 (years)", main="", xlim=c(1982,2016))
contour(x=date_decimal, wspec, y=a/52, ylim=c(0,5), add=T, drawlabels=F)
axis(side=1, at=1982:2016, labels=F)
abline(v=seq(1984,2014,by=2), col="grey80", lty=2)

plot(x=date_decimal, data_hfmd$CVA16, type="b", xaxt='n', xlab="", ylab="Counts of CV-A16", col="firebrick", xaxs="i", xlim=c(1982,2016))
axis(side=1, at=1982:2016, labels=F)
abline(v=seq(1984,2014,by=2), col="grey80", lty=2)

wfit <- cwt(sqrt(data_hfmd$CVA16), no, nv, plot=F)
wspec <- Mod(wfit)
image(x=date_decimal, wspec, y=a/52, ylim=c(0,5), xaxt="n", xlab="", ylab="Period of CV-A16 (years)", main="", xlim=c(1982,2016))
contour(x=date_decimal, wspec, y=a/52, ylim=c(0,5), add=T, drawlabels=F)
axis(side=1, at=1982:2016, labels=F)
abline(v=seq(1984,2014,by=2), col="grey80", lty=2)

plot(x=date_decimal, data_hfmd$HFMD_per_sentinel * data_hfmd$Sentinels, type="b", xaxt='n', xlab="", ylab="Counts of HFMD", xaxs="i", xlim=c(1982,2016))
axis(side=1, at=1982:2016, labels=F)
abline(v=seq(1984,2014,by=2), col="grey80", lty=2)

plot(x=date_decimal, data_hfmd$Other_CVA, type="b", xaxt='n', xlab="", ylab="Counts of CV-A6", col="cornflowerblue", xaxs="i", xlim=c(1982,2016))
axis(side=1, at=1982:2016, labels=1982:2016)
abline(v=seq(1984,2014,by=2), col="grey80", lty=2)

###################
## Make Figure 4 ##
###################

## Data set: 1982-2015
data_hfmd_1982_2015 <- data_hfmd[,c("Year", "Week_of_year", "EVA71", "CVA16")]
rownames(data_hfmd_1982_2015) <- 1:nrow(data_hfmd_1982_2015)
data_hfmd_1982_2015$Year <- as.numeric(as.character(data_hfmd_1982_2015$Year))
data_hfmd_1982_2015$Week_of_year <- as.numeric(as.character(data_hfmd_1982_2015$Week_of_year))
data_hfmd_1982_2015$EVA71 <- as.numeric(as.character(data_hfmd_1982_2015$EVA71))
data_hfmd_1982_2015$CVA16 <- as.numeric(as.character(data_hfmd_1982_2015$CVA16))

## Go year by year
years_1982_2015 <- unique(data_hfmd_1982_2015$Year)
n_years_1982_2015 <- length(years_1982_2015)

## CV-A16 size
reds <- brewer.pal(4, "Spectral")
CVA16_size_continuous_1982_2015 <- data.frame(year=years_1982_2015, size=NA, quartile=NA, col=NA)
for(i in 1:n_years_1982_2015) {
	
	CVA16_size_continuous_1982_2015$size[i] <- sum(data_hfmd_1982_2015$CVA16[which(data_hfmd_1982_2015$Year==years_1982_2015[i])])
	
}
CVA16_size_continuous_1982_2015$quartile <- cut(CVA16_size_continuous_1982_2015$size, breaks=quantile(CVA16_size_continuous_1982_2015$size, probs=seq(0,1,by=0.25)), labels=1:4, include.lowest=TRUE)
CVA16_size_continuous_1982_2015$col <- reds[CVA16_size_continuous_1982_2015$quartile]

## EV-A71 size
greens <- brewer.pal(4, "Spectral")
EVA71_size_continuous_1982_2015 <- data.frame(year=years_1982_2015, size=NA, quartile=NA, col=NA)
for(i in 1:n_years_1982_2015) {
	
	EVA71_size_continuous_1982_2015$size[i] <- sum(data_hfmd_1982_2015$EVA71[which(data_hfmd_1982_2015$Year==years_1982_2015[i])])
	
}
EVA71_size_continuous_1982_2015$quartile <- cut(EVA71_size_continuous_1982_2015$size, breaks=quantile(EVA71_size_continuous_1982_2015$size, probs=seq(0,1,by=0.25)), labels=1:4, include.lowest=TRUE)
EVA71_size_continuous_1982_2015$col <- greens[EVA71_size_continuous_1982_2015$quartile]

## Function to obtain COG from the data
COG_function <- function(data, indices) {
	
	data_incl <- data[indices]
	return(sum(data_incl)/length(data_incl))
	
}

## Function to obtain skewness from the data
skew_function <- function(data, indices) {
	
	data_incl <- data[indices]
	return(skewness(data_incl))
	
}

## Function to obtain all measures
stats_estimate <- function(strain, years_1982_2015, n_years_1982_2015) {

	names_tmp <- c("year", "COG_mean", "COG_lb", "COG_ub", "COG_se", "sum_thatstrain", "raw_mean", "raw_median", "raw_sd", "raw_skewness", "raw_kurtosis", "skewness", "skewness_lb", "skewness_ub", "skewness_se")
	COG <- as.data.frame(matrix(NA, ncol=length(names_tmp), nrow=n_years_1982_2015))
	names(COG) <- names_tmp
	
	for(i in 1:n_years_1982_2015) {
		
		print(i)
		
		## Which year?
		COG$year[i] <- years_1982_2015[i]
		data_year <- data_hfmd_1982_2015[which(data_hfmd_1982_2015$Year==years_1982_2015[i]),]
		
		## Re-arrange data to facilitate bootstrap interval
		data_boot <- rep(data_year[,"Week_of_year"], times=data_year[,strain])
		
		## Compute mean center of gravity in weeks
		COG$COG_mean[i] <- COG_function(data_boot)
		COG$raw_mean[i] <- mean(data_boot)
		COG$raw_median[i] <- median(data_boot)
		COG$raw_sd[i] <- sd(data_boot)
		COG$raw_skewness[i] <- moments::skewness(data_boot, na.rm=TRUE)
		COG$raw_kurtosis[i] <- moments::kurtosis(data_boot, na.rm=TRUE)
	
		## Bootstrapping with 10000 replications 
		results <- boot(data=data_boot, statistic=COG_function, R=10000)
		COG$COG_se[i] <- sd(results$t)
		output <- boot.ci(results, type="norm")
		
		## Get 95% CI
		COG$COG_lb[i] <- ifelse(!is.null(output[4]$normal[2]), output[4]$normal[2], NA)
		COG$COG_ub[i] <- ifelse(!is.null(output[4]$normal[3]), output[4]$normal[3], NA)
		
		## Get # of cases that year
		COG$sum_thatstrain[i] <- sum(data_year[,strain])

		## Compute mean skewness in weeks
		COG$skewness[i] <- skew_function(data_boot)
		
		## Bootstrapping with 10000 replications
		results2 <- tryCatch({ boot(data=data_boot, statistic=skew_function, R=10000) }, error=function(e) { return(NA) })
		COG$skewness_se[i] <- sd(results2$t)
		output2 <- tryCatch({ boot.ci(results2, type="norm") }, error=function(e) { return(NA) })
		
		## Get 95% CI
		COG$skewness_lb[i] <- tryCatch({ ifelse(!is.null(output2[4]$normal[2]), output2[4]$normal[2], NA) }, error=function(e) { return(NA) })
		COG$skewness_ub[i] <- tryCatch({ ifelse(!is.null(output2[4]$normal[3]), output2[4]$normal[3], NA) }, error=function(e) { return(NA) })
		
	}
	
	return(COG)
	
}

## Estimate everything
COG_EVA71 <- stats_estimate(strain="EVA71", years_1982_2015=years_1982_2015, n_years_1982_2015=n_years_1982_2015)
COG_CVA16 <- stats_estimate(strain="CVA16", years_1982_2015=years_1982_2015, n_years_1982_2015=n_years_1982_2015)

## Pull data together
data_plot_COG_EVA71 <- data.frame(COG_EVA71, quartile=CVA16_size_continuous_1982_2015$quartile)
data_plot_COG_CVA16 <- data.frame(COG_CVA16, quartile=EVA71_size_continuous_1982_2015$quartile)

## Which vars to plot?
x_var <- "factor(quartile)"
y1_var <- "raw_mean"
y2_var <- "raw_sd"
y3_var <- "raw_skewness"
y4_var <- "raw_kurtosis"
y5_var <- "raw_mean-raw_median"
fill_var <- "quartile"

data_plot_1982_2015 <- data.frame(CVA16=CVA16_size_continuous_1982_2015$size, EVA71=EVA71_size_continuous_1982_2015$size, logCVA16=log(CVA16_size_continuous_1982_2015$size), logEVA71=log(EVA71_size_continuous_1982_2015$size))

## Plot
p1 <- ggplot(data_plot_1982_2015, aes(CVA16, EVA71)) +
	theme_bw() +
	xlab("CV-A16") +
	ylab("EV-A71") +
	geom_point(size=2) +
	theme(plot.title = element_text(hjust=0.5)) +
	ggtitle("Annual detections")

p2 <- ggplot(data_plot_COG_EVA71, aes_string(x=x_var, y=y1_var, fill=fill_var)) +
	geom_hline(yintercept=30, size=2, colour="lightgray") +
	geom_boxplot() +
	scale_fill_manual(values=brewer.pal(4, "Spectral")) +
	scale_x_discrete(labels=c("Lowest 25%", "Q2", "Q3", "Highest 25%")) +
	theme_bw() +
	theme(legend.position="none", plot.title=element_text(hjust=0.5)) +
	xlab("Size of CV-A16") +
	ylab("Week") +
	ggtitle("COG of EV-A71")

p3 <- ggplot(data_plot_COG_CVA16, aes_string(x=x_var, y=y1_var, fill=fill_var)) +
	geom_hline(yintercept=30, size=2, colour="lightgray") +
	geom_boxplot() +
	scale_fill_manual(values=brewer.pal(4, "Spectral")) +
	scale_x_discrete(labels=c("Lowest 25%", "Q2", "Q3", "Highest 25%")) + 
	theme_bw() +
	theme(legend.position="none", plot.title=element_text(hjust=0.5)) +
	xlab("Size of EV-A71") +
	ylab("Week") +
	ggtitle("COG of CV-A16")

p4 <- ggplot(data_plot_COG_EVA71, aes_string(x=x_var, y=y3_var, fill=fill_var)) +
	geom_hline(yintercept=0, size=2, colour="lightgray") +
	geom_boxplot() +
	scale_fill_manual(values=brewer.pal(4, "Spectral")) +
	scale_x_discrete(labels=c("Lowest 25%", "Q2", "Q3", "Highest 25%")) +
	theme_bw() +
	theme(legend.position="none", plot.title=element_text(hjust=0.5)) +
	xlab("Size of CV-A16") +
	ylab("Skewness") +
	ggtitle("Skewness of EV-A71")

p5 <- ggplot(data_plot_COG_CVA16, aes_string(x=x_var, y=y3_var, fill=fill_var)) +
	geom_hline(yintercept=0, size=2, colour="lightgray") + 
	geom_boxplot() +
	scale_fill_manual(values=brewer.pal(4, "Spectral")) +
	scale_x_discrete(labels=c("Lowest 25%", "Q2", "Q3", "Highest 25%")) +
	theme_bw() +
	theme(legend.position="none", plot.title=element_text(hjust=0.5)) +
	xlab("Size of EV-A71") +
	ylab("Skewness") +
	ggtitle("Skewness of CV-A16")

grid.arrange(p1, grid.rect(gp=gpar(col="white")), p2, p3, p4, p5, ncol=2)
