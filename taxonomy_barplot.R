library(cowplot)
library(phyloseq)
library(MetBrewer)
library(ggplot2)
library(vegan)
library(dplyr)
library(xlsx)
library(paletteer)
library(Polychrome)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets")

#read in all MAG tables from CoverM
MAG_abund_table_raw <- read.xlsx2("CoverM_GTDB_joined.xlsx", 
                                   sheetIndex=1,row.names = 1)

MAG_abund_table_raw <- mutate_all(MAG_abund_table_raw, function(x) as.numeric(as.character(x)))
MAG_matrix <- as.matrix(MAG_abund_table_raw)

tax_mat_val <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=5, row.names=1)
tax_mat_val <- as.matrix(tax_mat_val)
meta <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=3, row.names=1)

OTU = otu_table(MAG_matrix, taxa_are_rows = TRUE)
TAX2 = tax_table(tax_mat_val)
samples = sample_data(meta)

#make phyloseq object
MAG_phyloseq = phyloseq(OTU, TAX2, samples)

#created a color palette in taxonomy_barplot_denitrification_v3.R, need to read it in
palette_34 <- read.csv("P34.csv", header=TRUE)
P34 <- c(palette_34[,2])
PhylumList <- c(palette_34[,1])
names(P34) <- PhylumList

#subset phyloseq object by Location
Arabian <- subset_samples(MAG_phyloseq, Location=="Arabian Sea")
ETNP <- subset_samples(MAG_phyloseq, Location=="Eastern Tropical North Pacific Ocean")
ETSP <- subset_samples(MAG_phyloseq, Location=="Eastern Tropical South Pacific Ocean")

#for datasets with particle fractions, remove all particle and filtered fractions
ETNP_bulk <- subset_samples(ETNP, type=="untreated seawater")
ETSP_bulk <- subset_samples(ETSP, type=="untreated seawater")

#for datasets with particle fractions, subset all particle and filtered fractions
ETNP_particle <- subset_samples(ETNP, type=="particle")
ETSP_particle <- subset_samples(ETSP, type=="particle")

desired_AS_order <- list("AMAL9", "AMAL10", "AMAL11", "AMAL12")

desired_ETNP_order <- list("AMAL5","AMAL17","Glass_SRR1509790","Glass_SRR1509797","AMAL18","AMAL1","AMAL13",
                           "Fuchsman_SRR4465037","Tsementzi_SRR3718412","Fuchsman_SRR4465034","Glass_SRR1509798",
                           "Glass_SRR1509792","Fuchsman_SRR4465025","AMAL14","Fuchsman_SRR4465024","Fuchsman_SRR4465031",
                           "Glass_SRR1509793","Fuchsman_SRR4465027","AMAL2","Fuchsman_SRR4465026","Fuchsman_SRR4465030",
                           "Fuchsman_SRR4465032","Tsementzi_SRR3718413","Glass_SRR1509794","Glass_SRR1509799",
                           "Fuchsman_SRR4465029","Fuchsman_SRR4465033","Fuchsman_SRR4465028","Fuchsman_SRR4465036","AMAL7",
                           "AMAL3","AMAL15","AMAL8","AMAL20","Fuchsman_SRR4465035","Glass_SRR1509796","Glass_SRR1509800")

desired_ETNP_bulk_order <- list("AMAL5","AMAL17","Glass_SRR1509790","Glass_SRR1509797","AMAL18","AMAL1","AMAL13",
                           "Fuchsman_SRR4465037","Tsementzi_SRR3718412","Fuchsman_SRR4465034","Glass_SRR1509798",
                           "Glass_SRR1509792","Fuchsman_SRR4465025","AMAL14","Fuchsman_SRR4465024",
                           "Glass_SRR1509793","Fuchsman_SRR4465027","AMAL2","Fuchsman_SRR4465026",
                           "Tsementzi_SRR3718413","Glass_SRR1509794","Glass_SRR1509799",
                           "Fuchsman_SRR4465029","Fuchsman_SRR4465028","Fuchsman_SRR4465036","AMAL7",
                           "AMAL3","AMAL15","AMAL8","AMAL20","Fuchsman_SRR4465035","Glass_SRR1509796","Glass_SRR1509800")

desired_ETSP_order <- list("Stewart_SRR304684","Stewart_SRR304671","Stewart_SRR064444","Stewart_SRR304672","Stewart_SRR304674",
                            "Stewart_SRR304656","Stewart_SRR070081","Ganesh_SRR960580","Ganesh_SRR961671","Stewart_SRR064448",
                            "Stewart_SRR304673","Stewart_SRR304680","Ganesh_SRR961673","Ganesh_SRR961675","Stewart_SRR070084",
                            "Stewart_SRR064450","Stewart_SRR070082","Ganesh_SRR961676","Ganesh_SRR961677","Stewart_SRR304668",
                            "Stewart_SRR304683","Ganesh_SRR961679","Ganesh_SRR961680")

desired_ETSP_bulk_order <- list("Stewart_SRR304684","Stewart_SRR304671","Stewart_SRR064444","Stewart_SRR304672","Stewart_SRR304674",
                                "Stewart_SRR304656","Stewart_SRR070081","Ganesh_SRR960580","Stewart_SRR064448",
                                "Stewart_SRR304673","Stewart_SRR304680","Ganesh_SRR961673","Stewart_SRR070084",
                                "Stewart_SRR064450","Stewart_SRR070082","Ganesh_SRR961676","Stewart_SRR304668",
                                "Stewart_SRR304683","Ganesh_SRR961679")

relabel_ETSP_bulk <- c('St_2008_15m','St_2008_35m','St_2008_50m','St_2008_50m','St_2008_50m',
                       'St_2008_65m','St_2008_70m','Ga_2010_70m','St_2008_110m','St_2008_110m','St_2008_110m',
                       'Ga_2010_110m','St_2008_150m','St_2008_200m','St_2008_200m','Ga_2010_200m','St_2008_500m',
                       'St_2008_800m','Ga_2010_1000m','Ga_2010_1000m')
relabel_ETNP_bulk <- c('2016_10m','2018_16m','Gl_2013_30m','Gl_2013_30m','2018_45m','2016_53m','2018_60m','Fu_2012_60m','Ts_2013_68m','Fu_2012_70m',
                       'Gl_2013_80m','Gl_2013_85m','Fu_2012_90m','2018_95m','Fu_2012_100m','Gl_2013_100m','Fu_2012_110m','2016_120m','Fu_2012_120m',
                       'Ts_2013_120m','Gl_2013_125m','Gl_2013_125m','Fu_2012_140m','Fu_2012_160m','Fu_2012_180m','2016_185m',
                       '2016_200m','2018_200m','2016_215m','2018_250m','Fu_2012_300m','Gl_2013_300m','Gl_2013_300m')
relabel_AS <- c('2007_130m','2007_150m','2007_200m','2007_400m')

#plot for Arabian Sea
all_AS_plot <- plot_bar(Arabian, fill = "Phylum") +  scale_fill_manual(values=P34) + 
  xlab("Arabian Sea") + ylab("Relative Abundance (%)")+
  geom_bar( stat="identity", position="stack") +scale_x_discrete(labels = relabel_AS)
all_AS_plot$data$Sample <- factor(all_AS_plot$data$Sample, levels = desired_AS_order)
all_AS_plot <- (all_AS_plot) + theme(legend.position = "none")+
  theme(axis.text.x=element_text(angle=45,hjust=1)) + theme(text = element_text(size = 7))  
print(all_AS_plot)

all_ETNP_plot <- plot_bar(ETNP, fill = "Phylum") +  scale_fill_manual(values=P34) + 
  geom_bar( stat="identity", position="stack") +
  xlab("ETNP") + ylab("Relative Abundance (%)")
all_ETNP_plot$data$Sample <- factor(all_ETNP_plot$data$Sample, levels = desired_ETNP_order)
#print(all_ETNP_plot)

all_ETNP_bulk_plot <- plot_bar(ETNP_bulk, fill = "Phylum") +  
  scale_fill_manual(values=P34) + xlab("ETNP") + ylab("Relative Abundance (%)") +
  geom_bar( stat="identity", position="stack")+scale_x_discrete(labels = relabel_ETNP_bulk)
all_ETNP_bulk_plot$data$Sample <- factor(all_ETNP_bulk_plot$data$Sample, levels = desired_ETNP_bulk_order)
all_ETNP_bulk_plot <- (all_ETNP_bulk_plot) + theme(legend.position = "none")+
  theme(axis.text.x=element_text(angle=45,hjust=1)) + theme(text = element_text(size = 7))  
#print(all_ETNP_bulk_plot)

all_ETSP_plot <- plot_bar(ETSP, fill = "Phylum") +  scale_fill_manual(values=P34) + 
  geom_bar( stat="identity", position="stack")+ xlab("ETSP") + ylab("Relative Abundance (%)") +
  scale_x_discrete(labels = relabel_ETSP_bulk)
all_ETSP_plot$data$Sample <- factor(all_ETSP_plot$data$Sample, levels = desired_ETSP_order)
#print(all_ETSP_plot)

all_ETSP_bulk_plot <- plot_bar(ETSP_bulk, fill = "Phylum") + xlab("ETSP") +
  ylab("Relative Abundance (%)") + 
  scale_fill_manual(values=P34) + scale_x_discrete(labels = relabel_ETSP_bulk) +
  geom_bar( stat="identity", position="stack")
all_ETSP_bulk_plot$data$Sample <- factor(all_ETSP_bulk_plot$data$Sample, levels = desired_ETSP_bulk_order)
all_ETSP_bulk_plot <- (all_ETSP_bulk_plot) +theme(legend.key.size = unit(0.4, 'cm'),legend.position = "right")+ 
  theme(axis.text.x=element_text(angle=45,hjust=1)) + theme(text = element_text(size = 7))  
#print(all_ETSP_bulk_plot)

all_plots_phylum <- plot_grid(all_AS_plot,all_ETNP_bulk_plot, all_ETSP_bulk_plot, 
          align = 'h',
          ncol = 3,
          rel_widths = c(0.6,2.2,2.9))
all_plots_phylum

setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="Overall_CoverM_bulk_phylum.svg", plot=all_plots_phylum)
ggsave(filename="Overall_CoverM_bulk_phylum.png", plot=all_plots_phylum)

#test out some co-occurrence plots
FR_MAG_plot <- plot_net(MAG_phyloseq_HMMrun2, "bray", type="taxa", color = "Phylum",maxdist = 0.3, laymeth="fruchterman.reingold", point_label = NULL, point_size = 3, point_alpha = 1,
)+ scale_color_paletteer_d("pals::glasbey")
FR_MAG_plot

FR_ETNP_plot <- plot_net(ETNP_bulk1, "bray", type="taxa", color = "Phylum",maxdist = 0.3, laymeth="fruchterman.reingold", point_label = NULL, point_size = 3, point_alpha = 1,
)+ scale_color_paletteer_d("pals::glasbey")
FR_ETNP_plot

FR_ETSP_plot <- plot_net(ETSP_bulk1, "bray", type="taxa", color = "Phylum",maxdist = 0.3, laymeth="fruchterman.reingold", point_label = NULL, point_size = 3, point_alpha = 1,
)+ scale_color_paletteer_d("pals::glasbey")
FR_ETSP_plot

FR_Arabian_plot <- plot_net(Arabian1, "bray", type="taxa", color = "Phylum",maxdist = 0.3, laymeth="fruchterman.reingold", point_label = NULL, point_size = 3, point_alpha = 1,
)+ scale_color_paletteer_d("pals::glasbey")
FR_Arabian_plot

#try co-occurrence plots with nitrogen data
FR_MAG_plot2 <- plot_net(MAG_phyloseq_HMMrun2, "bray", type="taxa", color = "Phylum",maxdist = 0.3, laymeth="fruchterman.reingold", point_label = NULL, point_size = 3, point_alpha = 1,
)+ scale_color_paletteer_d("pals::glasbey")
FR_MAG_plot2


plot_net(ETNP_bulk1, "bray", type="taxa", color = "nirS_HMMrun2",shape = "specialization_HMMrun2", rescale=TRUE,
         maxdist = 0.3, laymeth="fruchterman.reingold", point_label = NULL, point_size = 3, point_alpha = 1,
)+ scale_color_paletteer_d("pals::glasbey")

plot_net(ETSP_bulk1, "bray", type="taxa", color = "nirS_HMMrun2",shape = "number_denitrification_genes_HMMrun2", rescale=TRUE,
         maxdist = 0.3, laymeth="fruchterman.reingold", point_label = NULL, point_size = 3, point_alpha = 1,
)+ scale_color_paletteer_d("pals::glasbey")

plot_net(Arabian1, "bray", type="taxa", color = "nifH",shape = "number_denitrification_genes_HMMrun2", rescale=TRUE,
         maxdist = 0.2, laymeth="fruchterman.reingold", point_label = NULL, point_size = 3, point_alpha = 1,
)+ scale_color_paletteer_d("pals::glasbey")


plot_net(ETNP_bulk2, "bray", type="taxa", color = "nosZ_val",shape = "specialization_validated", rescale=TRUE,
         maxdist = 0.3, laymeth="fruchterman.reingold", point_label = NULL, point_size = 3, point_alpha = 1,
)+ scale_color_paletteer_d("pals::glasbey")

plot_net(ETSP_bulk2, "bray", type="taxa", color = "nirS_val",shape = "number_denitrification_genes_validated", rescale=TRUE,
         maxdist = 0.3, laymeth="fruchterman.reingold", point_label = NULL, point_size = 3, point_alpha = 1,
)+ scale_color_paletteer_d("pals::glasbey")

plot_net(Arabian2, "bray", type="taxa", color = "specialization_validated",shape = "number_denitrification_genes_validated", rescale=TRUE,
         maxdist = 0.2, laymeth="fruchterman.reingold", point_label = NULL, point_size = 3, point_alpha = 1,
)+ scale_color_paletteer_d("pals::glasbey")
