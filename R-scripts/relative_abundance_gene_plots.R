library(tidyverse)
library(readxl)
library(ggplot2)
library(cowplot)
library(tidyr)
library(MetBrewer)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets") 

AMAL_profile_data <- read_excel("Data_for_R_profiles.xlsx", 
                                sheet="AMAL")
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

denite_df <- gather(AMAL_profile_data, key = gene, value = Percent, 
             c("narG_percent", "napA_percent", "nirK_percent","nirS_percent","nor_percent","nosZ_percent"))

denite_plot<-ggplot(denite_df, aes(x=Percent, y = Depth, group = gene, colour = gene))  + geom_point(size=2)+
  scale_color_manual(values=met.brewer("Klimt", 6))+
  scale_x_continuous(limits=c(0, 50), position = "bottom")+
  scale_y_reverse(limits=c(400,0), breaks=c(0, 50, 100, 150, 200, 250, 300, 350, 400))+
  labs(title="Denitrification Genes (%) with Depth")+ ylab("Depth (m)")+ xlab("Denitrification Genes (%)")+ mytheme+
  theme(legend.position = 'none')
print(denite_plot)
denite_plot_line<-denite_plot + geom_line(aes(group=gene),orientation = "y")
denite_plot_line
smooth_denite_plot<-denite_plot + geom_smooth(aes(group=gene), orientation = "y", method ="loess", span=0.75, se = FALSE,size=0.75)
smooth_denite_plot

other_df <- gather(AMAL_profile_data, key = gene, value = Percent, 
                    c("nxrA_percent", "amoA_percent", "nrfA_percent","nifH_percent","hzo_percent"))

other_plot<-ggplot(other_df, aes(x=Percent, y = Depth, group = gene, colour = gene))  + geom_point(size=2)+
  scale_color_manual(values=met.brewer("Juarez", 5))+
  scale_x_continuous(limits=c(0, 40), position = "bottom")+
  scale_y_reverse(limits=c(400,0), breaks=c(0, 50, 100, 150, 200, 250, 300, 350, 400))+
  labs(title="Other Nitrogen Genes (%) with Depth")+ ylab("Depth (m)")+ xlab("Other Nitrogen Genes (%)")+ mytheme +
  theme(legend.position = 'none')
print(other_plot)
other_plot_line<-other_plot + geom_line(aes(group=gene),orientation = "y")
other_plot_line
smooth_other_plot<-other_plot + geom_smooth(aes(group=gene), orientation = "y", method = "loess", span=0.75, se = FALSE,size=0.75)
smooth_other_plot

all_plots <- plot_grid(denite_plot,other_plot,
                              labels = "AUTO",
                              ncol = 2)
print(all_plots)

smooth_plots <- plot_grid(smooth_denite_plot,smooth_other_plot,
                       labels = "AUTO",
                       ncol = 2)
print(smooth_plots)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="Relative_abundance_gene_plots.svg", plot=all_plots)
ggsave(filename="Relative_abundance_gene_plots_smooth.svg", plot=smooth_plots)
