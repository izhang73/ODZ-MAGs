library(ggplot2) ## package for data visualization 
library(readxl) ##package for reading in excel files
library(cowplot) ##package for saving multiple plots on one page

##set working directory
setwd("/Users/irene/Desktop/Bess_Xin_Amal_Krona") 

##set theme
mytheme <- theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=10), # setting for the x axis text (color, size, ...)
        axis.text.y = element_text(colour="black", size=10), # setting for the y axis text
        axis.line   = element_line(colour="black", size=0.1),# setting for the axis line 
        axis.title  = element_text(size=12),     # setting for the title 
        panel.border = element_rect(colour = "black", fill=NA, size=0.4), # setting for the panel border
        axis.ticks  = element_line(colour="black", size=0.4),# setting for the axis ticks 
        panel.grid.major = element_blank(),                  # setting for the major grid of the plot panel 
        panel.grid.minor = element_blank()                   # setting for the minor grid of the plot panel
  )

##read in data from datasheet as AMAL_profile_data
AMAL_profile_data <- read_excel("/Users/irene/Desktop/Bess_Xin_Amal_Krona/Data_for_R_profiles.xlsx", 
                                     sheet="AMAL")

##plot nitrate and nitrite data on one plot
AMAL_NOx <- ggplot()+
  geom_point(data=AMAL_profile_data, aes(x=Nitrate, y=Depth, colour=Location_Year), pch=17, size= 2.5) +
  geom_path(data=AMAL_profile_data, aes(x=Nitrate, y=Depth, colour=Location_Year, group = Location_Year), linetype = "dashed")+
  geom_point(data=AMAL_profile_data, aes(x=Nitrite, y=Depth, colour=Location_Year), pch=19, size= 2.5, show.legend = FALSE) +
  geom_path(data=AMAL_profile_data, aes(x=Nitrite, y=Depth, colour=Location_Year, group = Location_Year), linetype = "dashed")+
  scale_x_continuous(limits=c(0, 30), position = "bottom")+
  scale_y_reverse(limits=c(400,0), breaks=c(0, 50, 100, 150, 200, 250, 300, 350, 400))+
  labs(title="AMAL NOx")+ ylab("Depth (m)")+ xlab("NOx (uM)")+ mytheme


##plot nitrate data only
AMAL_NO3 <- ggplot()+
  geom_point(data=AMAL_profile_data, aes(x=Nitrate, y=Depth, colour=Location_Year), pch=17, size= 2.5, show.legend = FALSE) +
  scale_x_continuous(limits=c(0, 30), position = "bottom")+
  scale_y_reverse(limits=c(400,0), breaks=c(0, 50, 100, 150, 200, 250, 300, 350, 400))+
  labs(title="AMAL NO3")+ ylab("Depth (m)")+ xlab("NO3 (uM)")+ mytheme

##plot nitrite data only
AMAL_NO2 <- ggplot()+
  geom_point(data=AMAL_profile_data, aes(x=Nitrite, y=Depth, colour=Location_Year), pch=19, size= 2.5) +
  scale_x_continuous(limits=c(0, 10), position = "bottom")+
  scale_y_reverse(limits=c(400,0), breaks=c(0, 50, 100, 150, 200, 250, 300, 350, 400))+
  labs(title="AMAL NO2")+ ylab("Depth (m)")+ xlab("NO2 (uM)")+ mytheme

##plot oxygen data
AMAL_O2 <- ggplot()+
  geom_point(data=AMAL_profile_data, aes(x=Oxygen, y=Depth, colour=Location_Year), pch=19, size= 2.5, show.legend = FALSE) +
  scale_x_continuous(limits=c(0, 3), position = "bottom")+
  scale_y_reverse(limits=c(400,0), breaks=c(0, 50, 100, 150, 200, 250, 300, 350, 400))+
  labs(title="Oxygen")+ ylab("Depth (m)")+ xlab("Probe O2")+ mytheme


#plot salinity data
AMAL_Salinity <- ggplot()+
  geom_point(data=AMAL_profile_data, aes(x=Salinity, y=Depth, colour=Location_Year), pch=19, size= 2.5, show.legend = FALSE) +
  scale_x_continuous(limits=c(32, 37), position = "bottom")+
  scale_y_reverse(limits=c(400,0), breaks=c(0, 50, 100, 150, 200, 250, 300, 350, 400))+
  labs(title="Salinity")+ ylab("Depth (m)")+ xlab("Salinity")+ mytheme

##plot temperature data
AMAL_Temperature <- ggplot()+
  geom_point(data=AMAL_profile_data, aes(x=Temperature, y=Depth, colour=Location_Year), pch=19, size= 2.5) +
  scale_x_continuous(limits=c(10, 22), position = "bottom")+
  scale_y_reverse(limits=c(400,0), breaks=c(0, 50, 100, 150, 200, 250, 300, 350, 400))+
  labs(title="Temperature")+ ylab("Depth (m)")+ xlab("Temperature (C)")+ mytheme

##save various combinations of the plots as pdf files
ggsave(filename="AMAL Profiles.pdf", width=8, height=5, 
       useDingbats=FALSE,
       plot_grid(AMAL_NOx, 
                 #AMAL_NO3, 
                 #AMAL_NO2, 
                 AMAL_O2,
                 #AMAL_Salinity, 
                 #AMAL_Temperature,
                 ncol = 2))

##save various combinations of the plots as pdf files
ggsave(filename="AMAL TS Profiles.pdf", width=5, height=5, 
       useDingbats=FALSE,
       plot_grid(#AMAL_NOx, 
                 #AMAL_NO3, 
                 #AMAL_NO2, 
                 #AMAL_O2,
                 AMAL_Salinity, 
                 AMAL_Temperature,
                 ncol = 2))

##save various combinations of the plots as pdf files
ggsave(filename="AMAL Grouped Profiles.pdf", width=16, height=5, 
       useDingbats=FALSE,
       plot_grid(AMAL_NOx, 
         #AMAL_NO3, 
         #AMAL_NO2, 
         AMAL_O2,
         AMAL_Salinity, 
         AMAL_Temperature,
         ncol = 4))
