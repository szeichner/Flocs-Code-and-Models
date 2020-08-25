#Script to generate log-normal distributed grain size distribution
#Written by: J. A. Nghiem
#Last edited: August 24, 2020

#Summary: This script generates a log-normal grain size distribution based on input distribution parameters.
#It creates grain size and concentration data in a user-specified number of grain size classes for input into sedimentation model.

#Load packages
library(dplyr)
library(data.table)
library(ggplot2)
library(here)

#Inputs
#Discretization parameters
nclass <- 28 #number of discrete grain size classes
grain_size_min <- 0.00002 #mm, lower grain size bound of smallest grain size class
grain_size_max <- 2 #mm, upper grain size bound of largest grain size class
#Set expectation and standard deviation of grain size distribution
ex <- 0.024 #mm, mean of distribution is 24 microns
sd <- 0.1 #mm, standard deviation of distribution is 100 microns
ssc <- 0.0002082547 #suspended sediment volumetric concentration from Mississippi River depth profiles from Jordan, 1965 
#ssc has units vol sediment/vol water; mass concentration is about 0.55 g/L (or kg/m^3) and divide by sediment density (assumed 2650 kg/m^3) to obtain ssc
#Path to csv file that will be created by the script
output_path <- here("Model", "parametric_gsd.csv")

#Calculations and processing below
var <- sd^2 #convert standard deviation to variance

#Calculate the log-normal distribution parameters (mu and sigma squared)
mu <- log(ex^2/sqrt(var+ex^2))
sgsq <- log((var/ex^2)+1)

#Generate a vector of grain sizes to evaluate the probability density
grain_sizes <- exp(seq(log(0.000001), log(2), length.out=1000)) #log-spaced from 10^-6 to 2 mm
#Create a data frame to calculate probability densities
gs_df <- data.frame(size=grain_sizes) %>%
  mutate(log_size=log(grain_sizes)) %>% #take the log of grain sizes
  mutate(density=dnorm(log_size, mean=mu, sd=sqrt(sgsq))) #log grain size is normally distributed

#Plot the grain size distribution in semilog space
ggplot(gs_df, aes(size, density))+
  geom_line()+
  #plot vertical lines to mark grain size transitions
  geom_vline(xintercept=c(0.002, 0.0625, 2), lty=2, col="red")+ #size classification of clay, silt, and sand
  geom_vline(xintercept=0.04, col="blue")+ #approximate floc threshold according to Lamb et al., 2020
  #make the horizontal axis log scale
  scale_x_log10(breaks=10^(-10:10), labels=scales::trans_format("log10", scales::math_format(10^.x)))+
  annotation_logticks(sides="b")+
  labs(x="grain size (mm)", y="probability density")+
  theme_bw()+
  theme(panel.grid=element_blank())

#Create data frame with discretized grain size classes and concentration partitioned into each class
#Generate log-spaced grain size bin boundaries
bbounds <- seq(log(grain_size_min), log(grain_size_max), length.out=nclass+1) #values in log space
lower <- bbounds[-length(bbounds)]; upper <- bbounds[-1]

discretized_df <- data.frame(lower, upper) %>%
  mutate(center=(lower+upper)/2) %>% #calculate bin centers as simple average of lower and upper bounds in log space
  mutate(frac=pnorm(upper, mean=mu, sd=sqrt(sgsq))-pnorm(lower, mean=mu, sd=sqrt(sgsq))) %>% #calculate the fraction in each bin
  mutate(frac=frac/sum(frac)) %>% #normalize fractions so that their sum is 1
  mutate(gsc=ssc*frac) %>% #partition total suspended sediment concentration into each grain size class
  mutate_at(vars("lower", "upper", "center"), exp) #take exponential of grain sizes to convert them back to mm

#Save grain size class data to a csv file
fwrite(discretized_df, output_path)
#Explanations of variables:
#lower: minimum grain size (mm) contained in a grain size class
#upper: maximum grain size (mm) contained in a grain size class
#center: representative grain size (mm) for a grain size class, calculated as the average of the lower and upper grain sizes in log space
#gsc: volumetric suspended sediment concentration in a grain size class
