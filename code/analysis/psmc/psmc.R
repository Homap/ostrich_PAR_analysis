#********************************************************************************
# Analysis of PSMC results of three ostrich subspecies
#********************************************************************************
# Written by Homa Papoli Yazdi - July 15 2020
#********************************************************************************

#********************************************************************************
# Function to read PSMC result and write out a table with years and Ne
#********************************************************************************
# Rescale the ith iteration result of PSMC, and make ready for plotting
# file: result file from PSMC
# i.iteration: the ith iteration
# mu: mutation rate
# s: bin size
# g: years per generation
# This function is written by Shenglin Liu, Apr 3, 2016
psmc.result<-function(file,i.iteration=25,mu=4.63e-09,s=100,g=4)
{
  X<-scan(file=black,what="",sep="\n",quiet=TRUE)
  
  START<-grep("^RD",X)
  END<-grep("^//",X)
  
  X<-X[START[i.iteration+1]:END[i.iteration+1]]
  
  TR<-grep("^TR",X,value=TRUE)
  RS<-grep("^RS",X,value=TRUE)
  
  write(TR,"temp.psmc.result")
  theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
  N0<-theta0/4/mu/s
  
  write(RS,"temp.psmc.result")
  a<-read.table("temp.psmc.result")
  Generation<-as.numeric(2*N0*a[,3])
  Ne<-as.numeric(N0*a[,4])
  
  file.remove("temp.psmc.result")
  
  n.points<-length(Ne)
  YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
              Generation[n.points])*g
  Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
        Ne[n.points])
  
  data.frame(YearsAgo,Ne)
}

#********************************************************************************
# Function to plot PSMC for each individual with 100 bootstraps
#********************************************************************************
# for the first run, run psmc.result() on the output and create the plot.
# for each file that starts with "round", we should run psmc.result() and then
# plot it with lines() for rounds 1 to 100. 

# Plot PSMC curve with bootstraps for P1878_107
black = '../data/psmc/P1878_110.autosome.psmc.out'
psmc.result.black <- psmc.result(black)

pdf("black.pdf")
plot(1,1,ylab = "Effective population size (x10^4)",
     xlab = "Years(g=4, mu=0.5e-8)", xlim = c(4, 7), ylim = c(0,90))
files<-grep("\\.psmc$",dir(),value=TRUE)
ne_data <- list()
year_data <- list()
for(i in files[1:99]) {
  print(i)
  d <- psmc.result(i,i.iteration=25,mu=4.63e-09,s=100,g=4)
  ne_data[[i]] <- d$Ne
  year_data[[i]] <- d$YearsAgo
  lines(log10(d$YearsAgo), d$Ne/10000,
        type="l",lwd=0.5, col="grey")  
}
lines(log10(psmc.result.black$YearsAgo), psmc.result.black$Ne/10000, type = "l", 
     col = "black", lwd=3)
dev.off()

# Plot PSMC curve with bootstraps for P1878_117
blue = '../data/psmc/P1878_121.autosome.psmc.out'
psmc.result.blue <- psmc.result(blue)
pdf("blue.pdf")
plot(1,1,ylab = "Effective population size (x10^4)",
     xlab = "Years(g=4, mu=0.5e-8)", xlim = c(4, 7), ylim = c(0,125))
files<-grep("\\.psmc$",dir("../../blue/P1878_118"),value=TRUE)
ne_data <- list()
year_data <- list()
for(i in files[1:98]) {
  print(i)
  d <- psmc.result(i,i.iteration=25,mu=4.63e-09,s=100,g=4)
  ne_data[[i]] <- d$Ne
  year_data[[i]] <- d$YearsAgo
  lines(log10(d$YearsAgo), d$Ne/10000,
        type="l",lwd=0.5, col="lightblue")
}
lines(log10(psmc.result.blue$YearsAgo), psmc.result.blue$Ne/10000, type = "l", 
      col = "blue", lwd=3)
dev.off()

# Plot PSMC curve with bootstraps for P1878_127
red = '../data/psmc/P1878_128.autosome.psmc.out'
psmc.result.red <- psmc.result(red)
pdf("red.pdf")
plot(1,1,ylab = "Effective population size (x10^4)",
     xlab = "Years(g=4, mu=0.5e-8)", xlim = c(4, 7), ylim = c(0,40))

files<-grep("\\.psmc$",dir(),value=TRUE)
ne_data <- list()
year_data <- list()
for(i in files[1:97]) {
  print(i)
  d <- psmc.result(i,i.iteration=25,mu=4.63e-09,s=100,g=4)
  ne_data[[i]] <- d$Ne
  year_data[[i]] <- d$YearsAgo
  lines(log10(d$YearsAgo), d$Ne/10000,
        type="l",lwd=0.5, col="pink")  
}
lines(log10(psmc.result.red$YearsAgo), psmc.result.red$Ne/10000, type = "l", 
      col = "red", lwd=3)
dev.off()

#********************************************************************************
# Function to plot PSMC for the three subspecies in one plot
#********************************************************************************
black = 'black/P1878_107.autosome.psmcfa.out'
blue = 'blue/P1878_119.autosome.psmcfa.out'
red = 'red/P1878_127.autosome.psmcfa.out'

psmc.result.black <- psmc.result(black)
psmc.result.blue <- psmc.result(blue)
psmc.result.red <- psmc.result(red)

pdf("three_species.pdf")
plot(log10(psmc.result.blue$YearsAgo), psmc.result.blue$Ne/10000, type = "l", 
     col = "blue", ylab = "Effective population size (x10^4)",
     xlab = "Years(g=4, mu=0.5e-8)", lwd=3)
lines(log10(psmc.result.black$YearsAgo), psmc.result.black$Ne/10000, col = "black", lwd=3)
lines(log10(psmc.result.red$YearsAgo), psmc.result.red$Ne/10000, col = "red", lwd=3)
dev.off()


#********************************************************************************
# Function to calculate mean and confidence interval (CI) of Ne 
#********************************************************************************
psmc_bootstrap <- function(ne_data, year_data) {
  files<-grep("\\.psmc$",dir(),value=TRUE)
  ne_data <- list()
  year_data <- list()
  for(i in files[1:1001]) {
    print(i)
    d <- psmc.result(i,i.iteration=25,mu=4.63e-09,s=100,g=4)
    ne_data[[i]] <- d$Ne
    year_data[[i]] <- d$YearsAgo
  }
  ne_year_list <- list(ne_data, year_data)
}

# turn your list into a dataframe
black_CI <- psmc_bootstrap(ne_black, year_black)
blue_CI <- psmc_bootstrap(ne_blue, year_blue)
red_CI <- psmc_bootstrap(ne_red, year_red)
# turn your list into a dataframe
neDF_black <- data.frame(black_CI[[1]])
neDF_blue <- data.frame(blue_CI[[1]])
neDF_red <- data.frame(red_CI[[1]])

# To get the 95% confidence interval
# Lower bound for 95% CI is the 2.5th quantile:
quantile(bstrap,.025) # For each time point
quantile(bstrap,.975) # For each time point

outputCI <- data.frame(resampling_mean_year, resampling_mean_v, CI_lower_v, CI_upper_v)
write.table(outputCI, "black_out.csv", col.names = colnames(outputCI), 
            row.names = FALSE, sep = "\t", quote = FALSE)

# Check here to get 95% confidence interval
# http://www.stat.ucla.edu/~rgould/110as02/bsci


#*************************************************************************************************
#*************************************************************************************************
# PSMC Ne
f <- read.table("../../../Figures_2_4_5/Figure_5/PSMC/Black_P1878_110_Ne_year.txt", header = T)
head(f)
library(psych)
harmonic.mean(f$Ne, na.rm = T)
f <- read.table("../../../Figures_2_4_5/Figure_5/PSMC/Blue_P1878_121_Ne_year.txt", header = T)
head(f)
library(psych)
harmonic.mean(f$Ne, na.rm = T)
f <- read.table("../../../Figures_2_4_5/Figure_5/PSMC/Red_P1878_128_Ne_year.txt", header = T)
head(f)
library(psych)
harmonic.mean(f$Ne, na.rm = T)
