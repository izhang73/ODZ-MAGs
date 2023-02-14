library(tidyverse)
library(readxl)
library(ggplot2)
library(tidyr)
library(MetBrewer)
library(xlsx)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets") 

##set theme
mytheme <- theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=8), # setting for the x axis text (color, size, ...)
        axis.text.y = element_text(colour="black", size=8), # setting for the y axis text
        axis.line   = element_line(colour="black", size=0.1),# setting for the axis line 
        axis.title  = element_text(size=10),     # setting for the title 
        panel.border = element_rect(colour = "black", fill=NA, size=0.4), # setting for the panel border
        axis.ticks  = element_line(colour="black", size=0.4),# setting for the axis ticks 
        panel.grid.major = element_blank(),                  # setting for the major grid of the plot panel 
        panel.grid.minor = element_blank()                   # setting for the minor grid of the plot panel
  )

#read in Ganesh data and simplify column names
Ganesh_nuts <- read_excel("all_metadata.xlsx", 
                         sheet="Ganesh_nuts")
#add a column for NO3 and substract NO2 from NO3_NO2 to get NO3 only
#change BDL to 0 
Ganesh_nuts$NO2[Ganesh_nuts$NO2 == 'BDL'] <- 0
Ganesh_nuts$NO2 <- (as.numeric(Ganesh_nuts$NO2))
Ganesh_nuts$NO3_NO2[Ganesh_nuts$NO3_NO2 == 'BDL'] <- 0
Ganesh_nuts$NH4[Ganesh_nuts$NH4 == 'BDL'] <- 0
Ganesh_nuts$depth <- (as.numeric(Ganesh_nuts$depth))
Ganesh_nuts$NO3 <- (as.numeric(Ganesh_nuts$NO3_NO2) - as.numeric(Ganesh_nuts$NO2))

#calculate N* and change form of data from wide to long
Ganesh_nuts$Nstar <- (as.numeric(Ganesh_nuts$NO3_NO2)) - 16*(as.numeric(Ganesh_nuts$PO4)) + 2.9
Ganesh_df <- gather(Ganesh_nuts, key = nutrient, value = uM, 
                    c("NO3", "NO2", "Nstar"))
Ganesh_df$uM <- (as.numeric(Ganesh_df$uM))

#plot NO3, NO2, and N*
Ganesh_nuts_plot<-ggplot(Ganesh_df, aes(x=uM, y = depth, colour = nutrient))  + 
  geom_point(size=0.5) + geom_vline(xintercept=0,linetype=3)+
  scale_x_continuous(limits=c(-30, 50), position = "top")+
  scale_y_reverse(limits=c(1200,0))+
  labs(title="ETSP Nutrients with Depth")+ ylab("Depth (m)")+ xlab("Nutrients (µM)")+ mytheme+
  theme(legend.position = 'right')
Ganesh_nuts_plot
Ganesh_smooth<-Ganesh_nuts_plot + geom_smooth(aes(group=nutrient), orientation = "y", method = "loess", span=0.2, se = FALSE,size=0.75)+
  theme(legend.position = 'none')
Ganesh_smooth

#read in Stewart data and simplify column names
Stewart_nuts <- read_excel("all_metadata.xlsx", 
                           sheet="Stewart_nuts_O2")
colnames(Stewart_nuts)[3] ="depth"
colnames(Stewart_nuts)[5] ="NO3"
colnames(Stewart_nuts)[6] ="NO2"
colnames(Stewart_nuts)[9] ="PO4"

#calculate N* and change form of data from wide to long
Stewart_nuts$Nstar <- ((as.numeric(Stewart_nuts$NO3))+(as.numeric(Stewart_nuts$NO2))) - 16*(as.numeric(Stewart_nuts$PO4)) + 2.9
Stewart_df <- gather(Stewart_nuts, key = nutrient, value = uM, 
                    c("NO3", "NO2", "Nstar"))
Stewart_df$uM <- (as.numeric(Stewart_df$uM))

#plot NO3, NO2, and N*
Stewart_nuts_plot<-ggplot(Stewart_df, aes(x=uM, y = depth, colour = nutrient))  + geom_point(size=0.5) +
  scale_x_continuous(limits=c(-30, 50), position = "top")+
  geom_vline(xintercept=0,linetype=3)+
  scale_y_reverse(limits=c(1000,0))+
  labs(title="ETSP Nutrients with Depth")+ ylab("Depth (m)")+ xlab("Nutrients (µM)")+ mytheme+
  theme(legend.position = 'right')
Stewart_nuts_plot
Stewart_smooth<-Stewart_nuts_plot + geom_smooth(aes(group=nutrient), orientation = "y", se = FALSE, size=0.75)+theme(legend.position = 'none')
Stewart_smooth

#read in Fuchsman data and simplify column names
Fuchsman_nuts <- read_excel("all_metadata.xlsx", 
                          sheet="Fuchsman_nuts")
colnames(Fuchsman_nuts)[7] ="depth"
colnames(Fuchsman_nuts)[12] ="NO3"
colnames(Fuchsman_nuts)[13] ="NO2"
colnames(Fuchsman_nuts)[14] ="PO4"

#calculate N* and change form of data from wide to long
Fuchsman_nuts$Nstar <- ((as.numeric(Fuchsman_nuts$NO3))+(as.numeric(Fuchsman_nuts$NO2))) - 16*(as.numeric(Fuchsman_nuts$PO4)) + 2.9
Fuchsman_df <- gather(Fuchsman_nuts, key = nutrient, value = uM, 
                    c("NO3", "NO2", "Nstar"))
Fuchsman_df$uM <- (as.numeric(Fuchsman_df$uM))

#plot NO3, NO2, and N*
Fuchsman_nuts_plot<-ggplot(Fuchsman_df, aes(x=uM, y = depth, colour = nutrient))  + geom_point(size=0.5) +
  scale_x_continuous(limits=c(-30, 50), position = "top")+
  scale_y_reverse(limits=c(1200,0))+
  geom_vline(xintercept=0,linetype=3)+
  labs(title="ETNP Nutrients with Depth")+ ylab("Depth (m)")+ xlab("Nutrients (µM)")+ mytheme+
  theme(legend.position = 'right')
Fuchsman_nuts_plot
Fuchsman_smooth<-Fuchsman_nuts_plot + geom_smooth(aes(group=nutrient), orientation = "y", method = "loess", span=0.2, se = FALSE,size=0.75)+theme(legend.position = 'none')
Fuchsman_smooth

#read in Amal Arabian Sea data and simplify column names
Amal_AS_nuts <- read_excel("all_metadata.xlsx", 
                            sheet="Amal_Arabian_nuts_O2")
colnames(Amal_AS_nuts)[1] ="depth"
colnames(Amal_AS_nuts)[2] ="NO3"
colnames(Amal_AS_nuts)[3] ="NO2"
colnames(Amal_AS_nuts)[4] ="PO4"

#calculate N* and change form of data from wide to long
Amal_AS_nuts$Nstar <- ((as.numeric(Amal_AS_nuts$NO3))+(as.numeric(Amal_AS_nuts$NO2))) - 16*(as.numeric(Amal_AS_nuts$PO4)) + 2.9
Amal_AS_df <- gather(Amal_AS_nuts, key = nutrient, value = uM, 
                      c("NO3", "NO2", "Nstar"))
Amal_AS_df$uM <- (as.numeric(Amal_AS_df$uM))

#plot NO3, NO2, and N*
Amal_AS_nuts_plot<-ggplot(Amal_AS_df, aes(x=uM, y = depth, colour = nutrient))  + geom_point(size=0.5) +
  scale_x_continuous(limits=c(-30, 50), position = "top")+
  scale_y_reverse(limits=c(1200,0))+
  geom_vline(xintercept=0,linetype=3)+
  labs(title="Arabian Sea Nutrients with Depth")+ ylab("Depth (m)")+ xlab("Nutrients (µM)")+ mytheme+
  theme(legend.position = 'right')
Amal_AS_nuts_plot
Amal_AS_smooth<-Amal_AS_nuts_plot + geom_smooth(aes(group=nutrient), orientation = "y", method = "loess", span=0.2, se = FALSE,size=0.75)
Amal_AS_smooth

#read in Amal ETNP data and simplify column names
Amal_ETNP_nuts <- read_excel("all_metadata.xlsx", 
                           sheet="AMAL_ETNP_nuts")

#calculate N* and change form of data from wide to long
Amal_ETNP_nuts$Nstar <- ((as.numeric(Amal_ETNP_nuts$NO3))+(as.numeric(Amal_ETNP_nuts$NO2))) - 16*(as.numeric(Amal_ETNP_nuts$PO4)) + 2.9
Amal_ETNP_df <- gather(Amal_ETNP_nuts, key = nutrient, value = uM, 
                     c("NO3", "NO2", "Nstar"))
Amal_ETNP_df$uM <- (as.numeric(Amal_ETNP_df$uM))

#plot NO3, NO2, and N*
Amal_ETNP_nuts_plot<-ggplot(Amal_ETNP_df, aes(x=uM, y = depth, colour = nutrient))  + geom_point(size=0.5) +
  scale_x_continuous(limits=c(-20, 50), position = "top")+
  scale_y_reverse(limits=c(1000,0))+
  geom_vline(xintercept=0,linetype=3)+
  labs(title="ETNP Nutrients with Depth")+ ylab("Depth (m)")+ xlab("Nutrients (µM)")+ mytheme+
  theme(legend.position = 'right')
Amal_ETNP_nuts_plot

#read in Glass data and simplify names
#Glass data doesn't have NO2? so cannot calculate N*
Glass_nuts <- read_excel("all_metadata.xlsx", 
                             sheet="Glass_nuts")
colnames(Glass_nuts)[2] ="depth"
colnames(Glass_nuts)[4] ="NO3"
colnames(Glass_nuts)[5] ="PO4"
Glass_nuts$NO3 <- (as.numeric(Glass_nuts$NO3))

#plot NO3 only
Glass_nuts_plot<-ggplot(Glass_nuts, aes(x=NO3, y = depth))  + geom_point(size=0.25) +
  scale_x_continuous(limits=c(0, 50), position = "top")+
  scale_y_reverse(limits=c(1000,0))+
  labs(title="ETNP Nutrients with Depth")+ ylab("Depth (m)")+ xlab("Nutrients (µM)")+ mytheme+
  theme(legend.position = 'right')
Glass_nuts_plot

#read in Tsementzi data and simplify names
#skip plotting Tsementzi data as there are only 2 datapoints
Tsementzi_nuts <- read_excel("all_metadata.xlsx", 
                         sheet="Tsementzi_nuts")
colnames(Tsementzi_nuts)[2] ="depth"
colnames(Tsementzi_nuts)[4] ="NO3"
colnames(Tsementzi_nuts)[5] ="PO4"

nutrient_plots<-plot_grid(Ganesh_smooth,Stewart_nuts_plot,Fuchsman_smooth,
                          Amal_AS_smooth,Amal_ETNP_nuts_plot,Glass_nuts_plot,
                labels = "AUTO",
                ncol = 3,
                nrow = 2)
nutrient_plots

Ganesh_add_line<-Ganesh_smooth + geom_hline(yintercept=32, linetype="dashed", color = "dark grey")
Ganesh_add_line

Ganesh_annote <- Ganesh_add_line + annotate("rect", ymin = 64, ymax = 390, xmax = -30, xmin=50,
         alpha = .4,fill = "grey")
Ganesh_annote

Fuchsman_add_line<-Fuchsman_smooth + geom_hline(yintercept=52, linetype="dashed", color = "dark grey")
Fuchsman_add_line

Fuchsman_annote <- Fuchsman_add_line + annotate("rect", ymin = 93, ymax = 993, xmax = -30, xmin=50,
                                            alpha = .4,fill = "grey")
Fuchsman_annote

Amal_AS_add_line<-Amal_AS_smooth + geom_hline(yintercept=95, linetype="dashed", color = "dark grey")
Amal_AS_add_line

Amal_AS_annote <- Amal_AS_add_line + annotate("rect", ymin = 125, ymax = 900, xmax = -30, xmin=50,
                                                alpha = .4,fill = "grey")
Amal_AS_annote

representative_nutrient_plots <-plot_grid(Ganesh_annote,Fuchsman_annote, Amal_AS_annote,
                                                          rel_widths = c(1, 1,1.5),
                                                          ncol = 3)
representative_nutrient_plots

setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="representative_nutrient_plots_ODZ.svg", plot=representative_nutrient_plots)
ggsave(filename="representative_nutrient_plots_ODZ.png", plot=representative_nutrient_plots)
