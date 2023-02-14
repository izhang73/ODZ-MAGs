library(cowplot)
library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(xlsx)
library(tidyverse)
library(RColorBrewer)
library(Polychrome)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets")

#read in all MAG tables from CoverM
MAG_abund_table_raw <- read.xlsx2("CoverM_GTDB_joined.xlsx", 
                                  sheetIndex=1,row.names = 1)

MAG_abund_table_raw <- mutate_all(MAG_abund_table_raw, function(x) as.numeric(as.character(x)))
MAG_matrix <- as.matrix(MAG_abund_table_raw)

tax_mat_val <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=8, row.names=1)
tax_mat_val <- as.matrix(tax_mat_val)

meta <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=3, row.names=1)

OTU = otu_table(MAG_matrix, taxa_are_rows = TRUE)
TAX2 = tax_table(tax_mat_val)
samples = sample_data(meta)

#make phyloseq object
MAG_phyloseq_val = phyloseq(OTU, TAX2, samples)

#read in the hard-coded Phylum palette 
palette_34 <- read.csv("P34.csv", header=TRUE)
P34 <- c(palette_34[,2])
PhylumList <- c(palette_34[,1])
names(P34) <- PhylumList

#Gene plots showing motility categories, run only once
#P10 = createPalette(9,  c("#048A81", "#E0ACD5", "#E4CC37"))
#names(P10) <- NULL
#MotList = unique(tax_table(MAG_phyloseq_val)[,"chem_mot_aer"])
#MotList <- c(MotList[,1])
#MotList[9]<-"NA"
#names(P10) <- MotList
#write.csv(P10, file ="P10_chemot.csv")

#read in hard-coded Chemot palette
Palette_10 <- read.csv("P10_chemot.csv", header=TRUE)
P10 <- c(Palette_10[,2])
MotList <- c(Palette_10[,1])
names(P10) <- MotList

#subset phyloseq object by Location
Arabian1 <- subset_samples(MAG_phyloseq_val, Location=="Arabian Sea")
ETNP1 <- subset_samples(MAG_phyloseq_val, Location=="Eastern Tropical North Pacific Ocean")
ETSP1 <- subset_samples(MAG_phyloseq_val, Location=="Eastern Tropical South Pacific Ocean")

#for datasets with particle fractions, remove all particle and filtered fractions
ETNP_bulk1 <- subset_samples(ETNP1, type=="untreated seawater")
ETSP_bulk1 <- subset_samples(ETSP1, type=="untreated seawater")

#for datasets with particle fractions, subset all particle and filtered fractions
ETNP_particle1 <- subset_samples(ETNP1, type=="particle")
ETSP_particle1 <- subset_samples(ETSP1, type=="particle")

#separate out denitrifiers from non-denitrifiers
denitrifier1 <- subset_taxa(MAG_phyloseq_val, denitrifier_val=="yes")
nondenitrifier1 <- subset_taxa(MAG_phyloseq_val, denitrifier_val=="no")

#subset only denitrifiers from different regions
denite_AS1 <- subset_taxa(Arabian1, denitrifier_val=="yes")
denite_ETNP1 <- subset_taxa(ETNP1, denitrifier_val=="yes")
denite_ETNP_bulk1 <- subset_taxa(ETNP_bulk1, denitrifier_val=="yes")
denite_ETSP1 <- subset_taxa(ETSP1, denitrifier_val=="yes")
denite_ETSP_bulk1 <- subset_taxa(ETSP_bulk1, denitrifier_val=="yes")

#subset non-denitrifiers
nondenite_AS1 <- subset_taxa(Arabian1, denitrifier_val=="no")
nondenite_ETNP1 <- subset_taxa(ETNP1, denitrifier_val=="no")
nondenite_ETNP_bulk1 <- subset_taxa(ETNP_bulk1, denitrifier_val=="no")
nondenite_ETSP1 <- subset_taxa(ETSP1, denitrifier_val=="no")
nondenite_ETSP_bulk1 <- subset_taxa(ETSP_bulk1, denitrifier_val=="no")

#subset only denitrifiers from particles
denite_ETNP_particle1 <- subset_taxa(ETNP_particle1, denitrifier_val=="yes")
denite_ETSP_particle1 <- subset_taxa(ETSP_particle1, denitrifier_val=="yes")

#set desired sample orders for each subset
desired_ETSP_particle_order <- list("Ganesh_SRR961671","Ganesh_SRR961675","Ganesh_SRR961677","Ganesh_SRR961680")
desired_ETNP_particle_order <- list("Fuchsman_SRR4465031","Fuchsman_SRR4465030","Fuchsman_SRR4465033")
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

#particle plots
denite_ETNP_particle_plot <- plot_bar(denite_ETNP_particle1, fill = "chem_mot_aer") +  
  scale_fill_manual(values=P10)  + 
  geom_bar( stat="identity", position="stack")+ xlab("ETNP particles")
denite_ETNP_particle_plot$data$Sample <- factor(denite_ETNP_particle_plot$data$Sample, levels = desired_ETNP_particle_order)
#denite_ETNP_particle_plot <- (denite_ETNP_particle_plot) + theme(legend.position = "right")
denite_ETNP_particle_plot <- (denite_ETNP_particle_plot) + theme(legend.position = "none")
print(denite_ETNP_particle_plot)

denite_ETSP_particle_plot <- plot_bar(denite_ETSP_particle1, fill = "chem_mot_aer") + 
  geom_bar( stat="identity", position="stack") + xlab("ETSP particles")+scale_fill_manual(values=P10)
denite_ETSP_particle_plot$data$Sample <- factor(denite_ETSP_particle_plot$data$Sample, levels = desired_ETSP_particle_order)
#denite_ETSP_particle_plot <- (denite_ETSP_particle_plot) + theme(legend.position = "right")
print(denite_ETSP_particle_plot) 

denite_particle_taxonomy <- plot_grid(denite_ETNP_particle_plot,denite_ETSP_particle_plot, 
                                   labels = "AUTO",
                                   align = 'h',
                                   ncol = 2,
                                   rel_widths = c(1,3.5))
print(denite_particle_taxonomy) 

ETNP_particle_plot <- plot_bar(ETNP_particle1, fill = "chem_mot_aer") +  
  scale_fill_manual(values=P10)  + 
  geom_bar( stat="identity", position="stack")+ xlab("ETNP particles")
ETNP_particle_plot$data$Sample <- factor(ETNP_particle_plot$data$Sample, levels = desired_ETNP_particle_order)
#ETNP_particle_plot <- (ETNP_particle_plot) + theme(legend.position = "right")
ETNP_particle_plot <- (ETNP_particle_plot) + theme(legend.position = "none")
print(ETNP_particle_plot)

ETSP_particle_plot <- plot_bar(ETSP_particle1, fill = "chem_mot_aer") + 
  geom_bar( stat="identity", position="stack") + xlab("ETSP particles")+scale_fill_manual(values=P10)
ETSP_particle_plot$data$Sample <- factor(ETSP_particle_plot$data$Sample, levels = desired_ETSP_particle_order)
#denite_ETSP_particle_plot <- (ETSP_particle_plot) + theme(legend.position = "right")
print(ETSP_particle_plot) 

all_particle_taxonomy <- plot_grid(ETNP_particle_plot,ETSP_particle_plot, 
                                   labels = "AUTO",
                                   align = 'h',
                                   ncol = 2,
                                   rel_widths = c(1,3.5))
print(all_particle_taxonomy) 

#if needed to save particle plots - probably don't need this
ggsave(filename="Denitrifier_CoverM_particle_all.svg", plot=all_particle_taxonomy)

#overall plots and bulk plots
denite_AS_plot <- plot_bar(denite_AS1, fill = "chem_mot_aer")  +
  geom_bar( stat="identity", position="stack")+ xlab("Arabian Sea") +
  ylab("Relative Abundance (%)")+ ylim(0, 25)+
  scale_x_discrete(labels = relabel_AS)+scale_fill_manual(values=P10)
desired_AS_order <- list("AMAL9", "AMAL10", "AMAL11", "AMAL12")
denite_AS_plot$data$Sample <- factor(denite_AS_plot$data$Sample, levels = desired_AS_order)
denite_AS_plot <- (denite_AS_plot) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  theme(text = element_text(size = 7))  
print(denite_AS_plot)

nondenite_AS_plot <- plot_bar(nondenite_AS1, fill = "chem_mot_aer")  +
  geom_bar( stat="identity", position="stack")+ xlab("Arabian Sea") +
  ylab("Relative Abundance (%)")+ ylim(0, 40)+
  scale_x_discrete(labels = relabel_AS)+scale_fill_manual(values=P10)
desired_AS_order <- list("AMAL9", "AMAL10", "AMAL11", "AMAL12")
nondenite_AS_plot$data$Sample <- factor(nondenite_AS_plot$data$Sample, levels = desired_AS_order)
nondenite_AS_plot <- (nondenite_AS_plot) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  theme(text = element_text(size = 7))  
print(nondenite_AS_plot)

denite_ETNP_bulk_plot <- plot_bar(denite_ETNP_bulk1, fill = "chem_mot_aer") + 
  geom_bar( stat="identity", position="stack")+ xlab("ETNP ODZ") + scale_fill_manual(values=P10)+
  ylab("Relative Abundance (%)")+ ylim(0, 25) +scale_x_discrete(labels = relabel_ETNP_bulk)
denite_ETNP_bulk_plot$data$Sample <- factor(denite_ETNP_bulk_plot$data$Sample, levels = desired_ETNP_bulk_order)
denite_ETNP_bulk_plot <- (denite_ETNP_bulk_plot) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  theme(text = element_text(size = 7))  
print(denite_ETNP_bulk_plot)

nondenite_ETNP_bulk_plot <- plot_bar(nondenite_ETNP_bulk1, fill = "chem_mot_aer") + 
  geom_bar( stat="identity", position="stack")+ xlab("ETNP ODZ") + scale_fill_manual(values=P10)+
  ylab("Relative Abundance (%)")+ ylim(0, 40) +scale_x_discrete(labels = relabel_ETNP_bulk)
nondenite_ETNP_bulk_plot$data$Sample <- factor(nondenite_ETNP_bulk_plot$data$Sample, levels = desired_ETNP_bulk_order)
nondenite_ETNP_bulk_plot <- (nondenite_ETNP_bulk_plot) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  theme(text = element_text(size = 7))  
print(nondenite_ETNP_bulk_plot)

denite_ETSP_bulk_plot <- plot_bar(denite_ETSP_bulk1, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ xlab("ETSP ODZ") + ylab("Relative Abundance (%)")+ ylim(0, 25)+scale_x_discrete(labels = relabel_ETSP_bulk)
denite_ETSP_bulk_plot$data$Sample <- factor(denite_ETSP_bulk_plot$data$Sample, levels = desired_ETSP_bulk_order) 
denite_ETSP_bulk_plot <- (denite_ETSP_bulk_plot) + theme(legend.key.size = unit(0.4, 'cm'),legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  theme(text = element_text(size = 7))  
print(denite_ETSP_bulk_plot)

nondenite_ETSP_bulk_plot <- plot_bar(nondenite_ETSP_bulk1, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ xlab("ETSP ODZ") + ylab("Relative Abundance (%)")+ ylim(0, 40)+scale_x_discrete(labels = relabel_ETSP_bulk)
nondenite_ETSP_bulk_plot$data$Sample <- factor(nondenite_ETSP_bulk_plot$data$Sample, levels = desired_ETSP_bulk_order) 
nondenite_ETSP_bulk_plot <- (nondenite_ETSP_bulk_plot) + theme(legend.key.size = unit(0.4, 'cm'),legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  theme(text = element_text(size = 7))  
print(nondenite_ETSP_bulk_plot)

all_bulk_plots <- plot_grid(denite_AS_plot,denite_ETNP_bulk_plot, denite_ETSP_bulk_plot, 
                            align = 'h',
                            ncol = 3,
                            rel_widths = c(0.55,2.2,2.9))
print(all_bulk_plots)
setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="Denitrifier_CoverM_bulk_chemot.svg", plot=all_bulk_plots)

all_nondenite_bulk_plots <- plot_grid(nondenite_AS_plot,nondenite_ETNP_bulk_plot, nondenite_ETSP_bulk_plot, 
                            align = 'h',
                            ncol = 3,
                            rel_widths = c(0.55,2.2,2.9))
print(all_nondenite_bulk_plots)
ggsave(filename="Nondenitrifier_CoverM_bulk_chemot.svg", plot=all_nondenite_bulk_plots)


#Make motility and chemotaxis barplots for each gene individually for the ETSP
napA_ETSP_bulk <- subset_taxa(ETSP_bulk1, napA_val=="1")
denite_ETSP_bulk_napA <- plot_bar(napA_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP napA")
denite_ETSP_bulk_napA$data$Sample <- factor(denite_ETSP_bulk_napA$data$Sample, levels = desired_ETSP_bulk_order)
denite_ETSP_bulk_napA <- (denite_ETSP_bulk_napA) + coord_cartesian(ylim = c(0, 5))+ theme(text = element_text(size = 7))
denite_ETSP_bulk_napA <- (denite_ETSP_bulk_napA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETSP_bulk_napA <- (denite_ETSP_bulk_napA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETSP_bulk_napA)

narG_ETSP_bulk <- subset_taxa(ETSP_bulk1, narG_val=="1")
denite_ETSP_bulk_narG <- plot_bar(narG_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP narG")
denite_ETSP_bulk_narG$data$Sample <- factor(denite_ETSP_bulk_narG$data$Sample, levels = desired_ETSP_bulk_order)
denite_ETSP_bulk_narG <- (denite_ETSP_bulk_narG) + coord_cartesian(ylim = c(0, 10))+ theme(text = element_text(size = 7))
#denite_ETSP_bulk_narG <- (denite_ETSP_bulk_narG) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETSP_bulk_narG <- (denite_ETSP_bulk_narG) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETSP_bulk_narG)

nirK_ETSP_bulk <- subset_taxa(ETSP_bulk1, nirK_val=="1")
denite_ETSP_bulk_nirK <- plot_bar(nirK_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP nirK")
denite_ETSP_bulk_nirK$data$Sample <- factor(denite_ETSP_bulk_nirK$data$Sample, levels = desired_ETSP_bulk_order)
denite_ETSP_bulk_nirK <- (denite_ETSP_bulk_nirK) + coord_cartesian(ylim = c(0, 7))+ theme(text = element_text(size = 7))
#denite_ETSP_bulk_nirK <- (denite_ETSP_bulk_nirK) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETSP_bulk_nirK <- (denite_ETSP_bulk_nirK) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETSP_bulk_nirK)

nirS_ETSP_bulk <- subset_taxa(ETSP_bulk1, nirS_val=="1")
denite_ETSP_bulk_nirS <- plot_bar(nirS_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP nirS")
denite_ETSP_bulk_nirS$data$Sample <- factor(denite_ETSP_bulk_nirS$data$Sample, levels = desired_ETSP_bulk_order)
denite_ETSP_bulk_nirS <- (denite_ETSP_bulk_nirS) + coord_cartesian(ylim = c(0, 0.5))+ theme(text = element_text(size = 7))
#denite_ETSP_bulk_nirS <- (denite_ETSP_bulk_nirS) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETSP_bulk_nirS <- (denite_ETSP_bulk_nirS) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETSP_bulk_nirS)

norB_ETSP_bulk <- subset_taxa(ETSP_bulk1, nor_new=="1")
denite_ETSP_bulk_norB <- plot_bar(norB_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP norB")
denite_ETSP_bulk_norB$data$Sample <- factor(denite_ETSP_bulk_norB$data$Sample, levels = desired_ETSP_bulk_order)
denite_ETSP_bulk_norB <- (denite_ETSP_bulk_norB) + coord_cartesian(ylim = c(0, 13))+ theme(text = element_text(size = 7))
#denite_ETSP_bulk_norB <- (denite_ETSP_bulk_norB) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETSP_bulk_norB <- (denite_ETSP_bulk_norB) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETSP_bulk_norB)

nosZ_ETSP_bulk <- subset_taxa(ETSP_bulk1, nosZ_val=="1")
denite_ETSP_bulk_nosZ <- plot_bar(nosZ_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP nosZ")
denite_ETSP_bulk_nosZ$data$Sample <- factor(denite_ETSP_bulk_nosZ$data$Sample, levels = desired_ETSP_bulk_order)
denite_ETSP_bulk_nosZ <- (denite_ETSP_bulk_nosZ) + coord_cartesian(ylim = c(0, 2))+ theme(text = element_text(size = 7))
#denite_ETSP_bulk_nosZ <- (denite_ETSP_bulk_nosZ) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETSP_bulk_nosZ <- (denite_ETSP_bulk_nosZ) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETSP_bulk_nosZ)

denite_ETSP_all <- plot_grid(denite_ETSP_bulk_napA, denite_ETSP_bulk_narG, denite_ETSP_bulk_nirK, denite_ETSP_bulk_nirS,
                             denite_ETSP_bulk_norB, denite_ETSP_bulk_nosZ,
                                   labels = "AUTO",
                                   align = 'h',
                                   nrow = 6)

print(denite_ETSP_all)
#setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures")
#ggsave(filename="Denitrifier_CoverM_genes_chemot_all.svg", plot=denite_ETSP_all)

nxrA_ETSP_bulk <- subset_taxa(ETSP_bulk1, nxrA=="1")
denite_ETSP_bulk_nxrA <- plot_bar(nxrA_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP nxrA")
denite_ETSP_bulk_nxrA$data$Sample <- factor(denite_ETSP_bulk_nxrA$data$Sample, levels = desired_ETSP_bulk_order)
denite_ETSP_bulk_nxrA <- (denite_ETSP_bulk_nxrA) + coord_cartesian(ylim = c(0, 5))+ theme(text = element_text(size = 7))
#denite_ETSP_bulk_nxrA <- (denite_ETSP_bulk_nxrA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETSP_bulk_nxrA <- (denite_ETSP_bulk_nxrA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETSP_bulk_nxrA)

nrfA_ETSP_bulk <- subset_taxa(ETSP_bulk1, nrfA=="1")
denite_ETSP_bulk_nrfA <- plot_bar(nrfA_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP nrfA")
denite_ETSP_bulk_nrfA$data$Sample <- factor(denite_ETSP_bulk_nrfA$data$Sample, levels = desired_ETSP_bulk_order)
denite_ETSP_bulk_nrfA <- (denite_ETSP_bulk_nrfA) + coord_cartesian(ylim = c(0, 5))+ theme(text = element_text(size = 7))
#denite_ETSP_bulk_nrfA <- (denite_ETSP_bulk_nrfA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETSP_bulk_nrfA <- (denite_ETSP_bulk_nrfA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETSP_bulk_nrfA)

nifH_ETSP_bulk <- subset_taxa(ETSP_bulk1, nifH=="1")
denite_ETSP_bulk_nifH <- plot_bar(nifH_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP nifH")
denite_ETSP_bulk_nifH$data$Sample <- factor(denite_ETSP_bulk_nifH$data$Sample, levels = desired_ETSP_bulk_order)
denite_ETSP_bulk_nifH <- (denite_ETSP_bulk_nifH) + coord_cartesian(ylim = c(0, 2))+ theme(text = element_text(size = 7))
#denite_ETSP_bulk_nifH <- (denite_ETSP_bulk_nifH) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETSP_bulk_nifH <- (denite_ETSP_bulk_nifH) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETSP_bulk_nifH)

amoA_ETSP_bulk <- subset_taxa(ETSP_bulk1, amoA=="1")
denite_ETSP_bulk_amoA <- plot_bar(amoA_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP amoA")
denite_ETSP_bulk_amoA$data$Sample <- factor(denite_ETSP_bulk_amoA$data$Sample, levels = desired_ETSP_bulk_order)
denite_ETSP_bulk_amoA <- (denite_ETSP_bulk_amoA) + coord_cartesian(ylim = c(0, 2.5))+ theme(text = element_text(size = 7))
#denite_ETSP_bulk_amoA <- (denite_ETSP_bulk_amoA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETSP_bulk_amoA <- (denite_ETSP_bulk_amoA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETSP_bulk_amoA)

hzo_ETSP_bulk <- subset_taxa(ETSP_bulk1, hzo=="1")
denite_ETSP_bulk_hzo <- plot_bar(hzo_ETSP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP hzo")
denite_ETSP_bulk_hzo$data$Sample <- factor(denite_ETSP_bulk_hzo$data$Sample, levels = desired_ETSP_bulk_order)
denite_ETSP_bulk_hzo <- (denite_ETSP_bulk_hzo) + coord_cartesian(ylim = c(0, 5))+ theme(text = element_text(size = 7))
#denite_ETSP_bulk_hzo <- (denite_ETSP_bulk_hzo) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETSP_bulk_hzo <- (denite_ETSP_bulk_hzo) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETSP_bulk_hzo)

other_genes_ETSP <- plot_grid(denite_ETSP_bulk_nxrA, denite_ETSP_bulk_nrfA, denite_ETSP_bulk_nifH,
                              denite_ETSP_bulk_amoA, denite_ETSP_bulk_hzo,
                              labels = "AUTO",
                              align = 'h',
                              nrow = 5)
print(other_genes_ETSP)
ggsave(filename="CoverM_othergenes_chemot_ETSP.svg", plot=other_genes_ETSP)


#Make taxa barplots for each gene individually for the ETNP
napA_ETNP_bulk <- subset_taxa(ETNP_bulk1, napA_val=="1")
denite_ETNP_bulk_napA <- plot_bar(napA_ETNP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10) +
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP napA")
denite_ETNP_bulk_napA$data$Sample <- factor(denite_ETNP_bulk_napA$data$Sample, levels = desired_ETNP_bulk_order)
denite_ETNP_bulk_napA <- (denite_ETNP_bulk_napA) + coord_cartesian(ylim = c(0, 7))+ theme(text = element_text(size = 7))
#denite_ETNP_bulk_napA <- (denite_ETNP_bulk_napA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETNP_bulk_napA <- (denite_ETNP_bulk_napA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETNP_bulk_napA)

narG_ETNP_bulk <- subset_taxa(ETNP_bulk1, narG_val=="1")
denite_ETNP_bulk_narG <- plot_bar(narG_ETNP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10)  +
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP narG")
denite_ETNP_bulk_narG$data$Sample <- factor(denite_ETNP_bulk_narG$data$Sample, levels = desired_ETNP_bulk_order)
denite_ETNP_bulk_narG <- (denite_ETNP_bulk_narG) + coord_cartesian(ylim = c(0, 12))+ theme(text = element_text(size = 7))
#denite_ETNP_bulk_narG <- (denite_ETNP_bulk_narG) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETNP_bulk_narG <- (denite_ETNP_bulk_narG) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETNP_bulk_narG)

nirK_ETNP_bulk <- subset_taxa(ETNP_bulk1, nirK_val=="1")
denite_ETNP_bulk_nirK <- plot_bar(nirK_ETNP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+ 
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP nirK")
denite_ETNP_bulk_nirK$data$Sample <- factor(denite_ETNP_bulk_nirK$data$Sample, levels = desired_ETNP_bulk_order)
denite_ETNP_bulk_nirK <- (denite_ETNP_bulk_nirK) + coord_cartesian(ylim = c(0, 7))+ theme(text = element_text(size = 7))
#denite_ETNP_bulk_nirK <- (denite_ETNP_bulk_nirK) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETNP_bulk_nirK <- (denite_ETNP_bulk_nirK) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETNP_bulk_nirK)

nirS_ETNP_bulk <- subset_taxa(ETNP_bulk1, nirS_val=="1")
denite_ETNP_bulk_nirS <- plot_bar(nirS_ETNP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10) +  
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP nirS")
denite_ETNP_bulk_nirS$data$Sample <- factor(denite_ETNP_bulk_nirS$data$Sample, levels = desired_ETNP_bulk_order)
denite_ETNP_bulk_nirS <- (denite_ETNP_bulk_nirS) + coord_cartesian(ylim = c(0, 0.6))+ theme(text = element_text(size = 7))
#denite_ETNP_bulk_nirS <- (denite_ETNP_bulk_nirS) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETNP_bulk_nirS <- (denite_ETNP_bulk_nirS) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETNP_bulk_nirS)

norB_ETNP_bulk <- subset_taxa(ETNP_bulk1, norB_val=="1")
denite_ETNP_bulk_norB <- plot_bar(norB_ETNP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP norB")
denite_ETNP_bulk_norB$data$Sample <- factor(denite_ETNP_bulk_norB$data$Sample, levels = desired_ETNP_bulk_order)
denite_ETNP_bulk_norB <- (denite_ETNP_bulk_norB) + coord_cartesian(ylim = c(0, 4))+ theme(text = element_text(size = 7))
#denite_ETNP_bulk_norB <- (denite_ETNP_bulk_norB) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETNP_bulk_norB <- (denite_ETNP_bulk_norB) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETNP_bulk_norB)

nosZ_ETNP_bulk <- subset_taxa(ETNP_bulk1, nosZ_val=="1")
denite_ETNP_bulk_nosZ <- plot_bar(nosZ_ETNP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP nosZ")
denite_ETNP_bulk_nosZ$data$Sample <- factor(denite_ETNP_bulk_nosZ$data$Sample, levels = desired_ETNP_bulk_order)
denite_ETNP_bulk_nosZ <- (denite_ETNP_bulk_nosZ) + coord_cartesian(ylim = c(0, 2.5))+ theme(text = element_text(size = 7))
#denite_ETNP_bulk_nosZ <- (denite_ETNP_bulk_nosZ) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETNP_bulk_nosZ <- (denite_ETNP_bulk_nosZ) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETNP_bulk_nosZ)

denite_ETNP_all <- plot_grid(denite_ETNP_bulk_napA, denite_ETNP_bulk_narG, denite_ETNP_bulk_nirK, denite_ETNP_bulk_nirS,
                             denite_ETNP_bulk_norB, denite_ETNP_bulk_nosZ,
                                   labels = "AUTO",
                                   align = 'h',
                                   nrow = 6)
print(denite_ETNP_all)
#ggsave(filename="CoverM_denite_chemot_ETNP_all.svg", plot=denite_ETNP_all)

nxrA_ETNP_bulk <- subset_taxa(ETNP_bulk1, nxrA=="1")
denite_ETNP_bulk_nxrA <- plot_bar(nxrA_ETNP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+  
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP nxrA")
denite_ETNP_bulk_nxrA$data$Sample <- factor(denite_ETNP_bulk_nxrA$data$Sample, levels = desired_ETNP_bulk_order)
denite_ETNP_bulk_nxrA <- (denite_ETNP_bulk_nxrA) + coord_cartesian(ylim = c(0, 10))+ theme(text = element_text(size = 7))
#denite_ETNP_bulk_nxrA <- (denite_ETNP_bulk_nxrA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETNP_bulk_nxrA <- (denite_ETNP_bulk_nxrA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETNP_bulk_nxrA)

nrfA_ETNP_bulk <- subset_taxa(ETNP_bulk1, nrfA=="1")
denite_ETNP_bulk_nrfA <- plot_bar(nrfA_ETNP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10) +
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP nrfA")
denite_ETNP_bulk_nrfA$data$Sample <- factor(denite_ETNP_bulk_nrfA$data$Sample, levels = desired_ETNP_bulk_order)
denite_ETNP_bulk_nrfA <- (denite_ETNP_bulk_nrfA) + coord_cartesian(ylim = c(0, 7))+ theme(text = element_text(size = 7))
#denite_ETNP_bulk_nrfA <- (denite_ETNP_bulk_nrfA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETNP_bulk_nrfA <- (denite_ETNP_bulk_nrfA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETNP_bulk_nrfA)

nifH_ETNP_bulk <- subset_taxa(ETNP_bulk1, nifH=="1")
denite_ETNP_bulk_nifH <- plot_bar(nifH_ETNP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+  
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP nifH")
denite_ETNP_bulk_nifH$data$Sample <- factor(denite_ETNP_bulk_nifH$data$Sample, levels = desired_ETNP_bulk_order)
denite_ETNP_bulk_nifH <- (denite_ETNP_bulk_nifH) + coord_cartesian(ylim = c(0, 4))+ theme(text = element_text(size = 7))
#denite_ETNP_bulk_nifH <- (denite_ETNP_bulk_nifH) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETNP_bulk_nifH <- (denite_ETNP_bulk_nifH) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETNP_bulk_nifH)

amoA_ETNP_bulk <- subset_taxa(ETNP_bulk1, amoA=="1")
denite_ETNP_bulk_amoA <- plot_bar(amoA_ETNP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP amoA")
denite_ETNP_bulk_amoA$data$Sample <- factor(denite_ETNP_bulk_amoA$data$Sample, levels = desired_ETNP_bulk_order)
denite_ETNP_bulk_amoA <- (denite_ETNP_bulk_amoA) + coord_cartesian(ylim = c(0, 4))+ theme(text = element_text(size = 7))
#denite_ETNP_bulk_amoA <- (denite_ETNP_bulk_amoA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETNP_bulk_amoA <- (denite_ETNP_bulk_amoA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETNP_bulk_amoA)

hzo_ETNP_bulk <- subset_taxa(ETNP_bulk1, hzo=="1")
denite_ETNP_bulk_hzo <- plot_bar(hzo_ETNP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP hzo")
denite_ETNP_bulk_hzo$data$Sample <- factor(denite_ETNP_bulk_hzo$data$Sample, levels = desired_ETNP_bulk_order)
denite_ETNP_bulk_hzo <- (denite_ETNP_bulk_hzo) + coord_cartesian(ylim = c(0, 10))+ theme(text = element_text(size = 7))
#denite_ETNP_bulk_hzo <- (denite_ETNP_bulk_hzo) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_ETNP_bulk_hzo <- (denite_ETNP_bulk_hzo) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_ETNP_bulk_hzo)

other_genes_ETNP <- plot_grid(denite_ETNP_bulk_nxrA, denite_ETNP_bulk_nrfA, denite_ETNP_bulk_nifH,
                              denite_ETNP_bulk_amoA, denite_ETNP_bulk_hzo,
                              labels = "AUTO",
                              align = 'h',
                              nrow = 5)
print(other_genes_ETNP)
ggsave(filename="CoverM_othergenes_chemot_ETNP_all.svg", plot=other_genes_ETNP)


#Make motility and chemotaxis barplots for each gene individually for the Arabian Sea
napA_Arabian <- subset_taxa(Arabian1, napA_val=="1")
denite_Arabian_napA <- plot_bar(napA_Arabian, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+  
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian napA")
denite_Arabian_napA$data$Sample <- factor(denite_Arabian_napA$data$Sample, levels = desired_AS_order)
denite_Arabian_napA <- (denite_Arabian_napA) + coord_cartesian(ylim = c(0, 5))+ theme(text = element_text(size = 7))
#denite_Arabian_napA <- (denite_Arabian_napA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_Arabian_napA <- (denite_Arabian_napA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_Arabian_napA)

narG_Arabian <- subset_taxa(Arabian1, narG_val=="1")
denite_Arabian_narG <- plot_bar(narG_Arabian, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+ 
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian narG")
denite_Arabian_narG$data$Sample <- factor(denite_Arabian_narG$data$Sample, levels = desired_AS_order)
denite_Arabian_narG <- (denite_Arabian_narG) + coord_cartesian(ylim = c(0, 15))+ theme(text = element_text(size = 7))
#denite_Arabian_narG <- (denite_Arabian_narG) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_Arabian_narG <- (denite_Arabian_narG) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_Arabian_narG)

nirK_Arabian <- subset_taxa(Arabian1, nirK_val=="1")
denite_Arabian_nirK <- plot_bar(nirK_Arabian, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+ 
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian nirK")
denite_Arabian_nirK$data$Sample <- factor(denite_Arabian_nirK$data$Sample, levels = desired_AS_order)
denite_Arabian_nirK <- (denite_Arabian_nirK) + coord_cartesian(ylim = c(0, 5))+ theme(text = element_text(size = 7))
#denite_Arabian_nirK <- (denite_Arabian_nirK) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_Arabian_nirK <- (denite_Arabian_nirK) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_Arabian_nirK)

nirS_Arabian <- subset_taxa(Arabian1, nirS_val=="1")
denite_Arabian_nirS <- plot_bar(nirS_Arabian, fill = "chem_mot_aer") + scale_fill_manual(values=P10) +
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian nirS")
denite_Arabian_nirS$data$Sample <- factor(denite_Arabian_nirS$data$Sample, levels = desired_AS_order)
denite_Arabian_nirS <- (denite_Arabian_nirS) + coord_cartesian(ylim = c(0, 1.5))+ theme(text = element_text(size = 7))
#denite_Arabian_nirS <- (denite_Arabian_nirS) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_Arabian_nirS <- (denite_Arabian_nirS) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_Arabian_nirS)

norB_Arabian <- subset_taxa(Arabian1, nor_new=="1")
denite_Arabian_norB <- plot_bar(norB_Arabian, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+ 
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian norB")
denite_Arabian_norB$data$Sample <- factor(denite_Arabian_norB$data$Sample, levels = desired_AS_order)
denite_Arabian_norB <- (denite_Arabian_norB) + coord_cartesian(ylim = c(0, 1.5))+ theme(text = element_text(size = 7))
#denite_Arabian_norB <- (denite_Arabian_norB) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_Arabian_norB <- (denite_Arabian_norB) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_Arabian_norB)

nosZ_Arabian <- subset_taxa(Arabian1, nosZ_val=="1")
denite_Arabian_nosZ <- plot_bar(nosZ_Arabian, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian nosZ")
denite_Arabian_nosZ$data$Sample <- factor(denite_Arabian_nosZ$data$Sample, levels = desired_AS_order)
denite_Arabian_nosZ <- (denite_Arabian_nosZ) + coord_cartesian(ylim = c(0, 3.5))+ theme(text = element_text(size = 7))
#denite_Arabian_nosZ <- (denite_Arabian_nosZ) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_Arabian_nosZ <- (denite_Arabian_nosZ) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_Arabian_nosZ)

all_AS <- plot_grid(denite_Arabian_napA, denite_Arabian_narG, denite_Arabian_nirK,
                    denite_Arabian_nirS, denite_Arabian_norB,denite_Arabian_nosZ,
                              labels = "AUTO",
                              align = 'h',
                              nrow = 6)
print(all_AS)
ggsave(filename="CoverM_othergenes_chemot_Arabian_all.svg", plot=all_AS)


nxrA_Arabian <- subset_taxa(Arabian1, nxrA=="1")
denite_Arabian_nxrA <- plot_bar(nxrA_Arabian, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+ 
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian nxrA")
denite_Arabian_nxrA$data$Sample <- factor(denite_Arabian_nxrA$data$Sample, levels = desired_AS_order)
denite_Arabian_nxrA <- (denite_Arabian_nxrA) + coord_cartesian(ylim = c(0, 10))+ theme(text = element_text(size = 7))
#denite_Arabian_nxrA <- (denite_Arabian_nxrA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_Arabian_nxrA <- (denite_Arabian_nxrA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_Arabian_nxrA)

nrfA_Arabian <- subset_taxa(Arabian1, nrfA=="1")
denite_Arabian_nrfA <- plot_bar(nrfA_Arabian, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+  
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian nrfA")
denite_Arabian_nrfA$data$Sample <- factor(denite_Arabian_nrfA$data$Sample, levels = desired_AS_order)
denite_Arabian_nrfA <- (denite_Arabian_nrfA) + coord_cartesian(ylim = c(0, 7))+ theme(text = element_text(size = 7))
#denite_Arabian_nrfA <- (denite_Arabian_nrfA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_Arabian_nrfA <- (denite_Arabian_nrfA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_Arabian_nrfA)

nifH_Arabian <- subset_taxa(Arabian1, nifH=="1")
denite_Arabian_nifH <- plot_bar(nifH_Arabian, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+  
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian nifH")
denite_Arabian_nifH$data$Sample <- factor(denite_Arabian_nifH$data$Sample, levels = desired_AS_order)
denite_Arabian_nifH <- (denite_Arabian_nifH) + coord_cartesian(ylim = c(0, 1))+ theme(text = element_text(size = 7))
#denite_Arabian_nifH <- (denite_Arabian_nifH) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_Arabian_nifH <- (denite_Arabian_nifH) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_Arabian_nifH)

amoA_Arabian <- subset_taxa(Arabian1, amoA=="1")
denite_Arabian_amoA <- plot_bar(amoA_Arabian, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian amoA")
denite_Arabian_amoA$data$Sample <- factor(denite_Arabian_amoA$data$Sample, levels = desired_AS_order)
denite_Arabian_amoA <- (denite_Arabian_amoA) + coord_cartesian(ylim = c(0, 0.25))+ theme(text = element_text(size = 7))
#denite_Arabian_amoA <- (denite_Arabian_amoA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_Arabian_amoA <- (denite_Arabian_amoA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_Arabian_amoA)

hzo_Arabian <- subset_taxa(Arabian1, hzo=="1")
denite_Arabian_hzo <- plot_bar(hzo_Arabian, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian hzo")
denite_Arabian_hzo$data$Sample <- factor(denite_Arabian_hzo$data$Sample, levels = desired_AS_order)
denite_Arabian_hzo <- (denite_Arabian_hzo) + coord_cartesian(ylim = c(0, 5))+ theme(text = element_text(size = 7))
#denite_Arabian_hzo <- (denite_Arabian_hzo) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
denite_Arabian_hzo <- (denite_Arabian_hzo) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(denite_Arabian_hzo)

other_genes_AS <- plot_grid(denite_Arabian_nxrA, denite_Arabian_nrfA, denite_Arabian_nifH,denite_Arabian_amoA,denite_Arabian_hzo,
                            labels = "AUTO",
                            align = 'h',
                            nrow = 5)
print(other_genes_AS)
ggsave(filename="CoverM_othergenes_chemot_Arabian_all.svg", plot=other_genes_AS)


#transform into proportion diagrams - basically scale it all up to 100%
scaled_denite_ETSPbulk1 = transform_sample_counts(denite_ETSP_bulk1, function(x) (x / sum(x))*100 )
scaled_nondenite_ETSPbulk1 = transform_sample_counts(nondenite_ETSP_bulk1, function(x) (x / sum(x))*100 )
scaled_denite_ETNPbulk1 = transform_sample_counts(denite_ETNP_bulk1, function(x) (x / sum(x))*100 )
scaled_nondenite_ETNPbulk1 = transform_sample_counts(nondenite_ETNP_bulk1, function(x) (x / sum(x))*100 )
scaled_denite_AS = transform_sample_counts(denite_AS1, function(x) (x / sum(x))*100 )
scaled_nondenite_AS = transform_sample_counts(nondenite_AS1, function(x) (x / sum(x))*100 )

scaled_denite_AS_plot <- plot_bar(scaled_denite_AS, fill = "chem_mot_aer")  +
  geom_bar( stat="identity", position="stack")+ xlab("Arabian Sea") +
  ylab("Proportion of Denitrifier Community (%)")+ ylim(0, 102)+
  scale_x_discrete(labels = relabel_AS)+scale_fill_manual(values=P10)
desired_AS_order <- list("AMAL9", "AMAL10", "AMAL11", "AMAL12")
scaled_denite_AS_plot$data$Sample <- factor(scaled_denite_AS_plot$data$Sample, levels = desired_AS_order)
scaled_denite_AS_plot <- (scaled_denite_AS_plot) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  theme(text = element_text(size = 7))  
print(scaled_denite_AS_plot)

scaled_nondenite_AS_plot <- plot_bar(scaled_nondenite_AS, fill = "chem_mot_aer")  +
  geom_bar( stat="identity", position="stack")+ xlab("Arabian Sea") +
  ylab("Proportion of Non-denitrifier Community (%)")+ ylim(0, 102)+
  scale_x_discrete(labels = relabel_AS)+scale_fill_manual(values=P10)
desired_AS_order <- list("AMAL9", "AMAL10", "AMAL11", "AMAL12")
scaled_nondenite_AS_plot$data$Sample <- factor(scaled_nondenite_AS_plot$data$Sample, levels = desired_AS_order)
scaled_nondenite_AS_plot <- (scaled_nondenite_AS_plot) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  theme(text = element_text(size = 7))  
print(scaled_nondenite_AS_plot)

scaled_denite_ETNP_bulk_plot <- plot_bar(scaled_denite_ETNPbulk1, fill = "chem_mot_aer") + 
  geom_bar( stat="identity", position="stack")+ xlab("ETNP ODZ") + scale_fill_manual(values=P10)+
  ylab("Proportion of Denitrifier Community (%)")+ ylim(0, 102) +scale_x_discrete(labels = relabel_ETNP_bulk)
scaled_denite_ETNP_bulk_plot$data$Sample <- factor(scaled_denite_ETNP_bulk_plot$data$Sample, levels = desired_ETNP_bulk_order)
scaled_denite_ETNP_bulk_plot<- (scaled_denite_ETNP_bulk_plot) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  theme(text = element_text(size = 7))  
print(scaled_denite_ETNP_bulk_plot)

scaled_nondenite_ETNPbulk_plot <- plot_bar(scaled_nondenite_ETNPbulk1, fill = "chem_mot_aer") + 
  geom_bar( stat="identity", position="stack")+ xlab("ETNP ODZ") + scale_fill_manual(values=P10)+
  ylab("Proportion of Non-denitrifier Community (%)")+ ylim(0, 102) +scale_x_discrete(labels = relabel_ETNP_bulk)
scaled_nondenite_ETNPbulk_plot$data$Sample <- factor(scaled_nondenite_ETNPbulk_plot$data$Sample, levels = desired_ETNP_bulk_order)
scaled_nondenite_ETNPbulk_plot <- (scaled_nondenite_ETNPbulk_plot) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  theme(text = element_text(size = 7))  
print(scaled_nondenite_ETNPbulk_plot)

scaled_denite_ETSPbulk_plot <- plot_bar(scaled_denite_ETSPbulk1 , fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ xlab("ETSP ODZ") + ylab("Proportion of Denitrifier Community (%)")+ ylim(0, 102)+scale_x_discrete(labels = relabel_ETSP_bulk)
scaled_denite_ETSPbulk_plot$data$Sample <- factor(scaled_denite_ETSPbulk_plot$data$Sample, levels = desired_ETSP_bulk_order) 
scaled_denite_ETSPbulk_plot <- (scaled_denite_ETSPbulk_plot) + theme(legend.key.size = unit(0.4, 'cm'),legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  theme(text = element_text(size = 7))  
print(scaled_denite_ETSPbulk_plot)

scaled_nondenite_ETSPbulk_plot <- plot_bar(scaled_nondenite_ETSPbulk1, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ xlab("ETSP ODZ") + ylab("Proportion of Non-denitrifier Community (%)")+ ylim(0, 102)+scale_x_discrete(labels = relabel_ETSP_bulk)
scaled_nondenite_ETSPbulk_plot$data$Sample <- factor(scaled_nondenite_ETSPbulk_plot$data$Sample, levels = desired_ETSP_bulk_order) 
scaled_nondenite_ETSPbulk_plot <- (scaled_nondenite_ETSPbulk_plot) + theme(legend.key.size = unit(0.4, 'cm'),legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  theme(text = element_text(size = 7))  
print(scaled_nondenite_ETSPbulk_plot)

all_scaled_denite_plots <- plot_grid(scaled_denite_AS_plot,scaled_denite_ETNP_bulk_plot, scaled_denite_ETSPbulk_plot, 
                            align = 'h',
                            ncol = 3,
                            labels = c('A'),
                            rel_widths = c(0.55,2.4,2.7))
print(all_scaled_denite_plots)
setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="Denitrifier_CoverM_bulk_chemot_scaled.svg", plot=all_scaled_denite_plots)

all_scaled_nondenite_plots <- plot_grid(scaled_nondenite_AS_plot,scaled_nondenite_ETNPbulk_plot, scaled_nondenite_ETSPbulk_plot, 
                                      align = 'h',
                                      ncol = 3,
                                      labels = c('A'),
                                      rel_widths = c(0.55,2.4,2.7))
print(all_scaled_nondenite_plots)
ggsave(filename="Nondenitrifier_CoverM_bulk_chemot_scaled.svg", plot=all_scaled_nondenite_plots)


#Make scaled motility and chemotaxis barplots for each gene individually for the ETSP

napA_ETSP_bulk <- subset_taxa(ETSP_bulk1, napA_val=="1")
scaled_napA_ETSP_bulk = transform_sample_counts(napA_ETSP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETSP_bulk_napA <- plot_bar(scaled_napA_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP napA")
scaled_denite_ETSP_bulk_napA$data$Sample <- factor(scaled_denite_ETSP_bulk_napA$data$Sample, levels = desired_ETSP_bulk_order)
scaled_denite_ETSP_bulk_napA <- (scaled_denite_ETSP_bulk_napA) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETSP_bulk_napA <- (scaled_denite_ETSP_bulk_napA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETSP_bulk_napA <- (scaled_denite_ETSP_bulk_napA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETSP_bulk_napA)

narG_ETSP_bulk <- subset_taxa(ETSP_bulk1, narG_val=="1")
scaled_narG_ETSP_bulk = transform_sample_counts(narG_ETSP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETSP_bulk_narG <- plot_bar(scaled_narG_ETSP_bulk , fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP narG")
scaled_denite_ETSP_bulk_narG$data$Sample <- factor(scaled_denite_ETSP_bulk_narG$data$Sample, levels = desired_ETSP_bulk_order)
scaled_denite_ETSP_bulk_narG <- (scaled_denite_ETSP_bulk_narG) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETSP_bulk_narG <- (scaled_denite_ETSP_bulk_narG) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETSP_bulk_narG <- (scaled_denite_ETSP_bulk_narG) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETSP_bulk_narG)

nirK_ETSP_bulk <- subset_taxa(ETSP_bulk1, nirK_val=="1")
scaled_nirK_ETSP_bulk = transform_sample_counts(nirK_ETSP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETSP_bulk_nirK <- plot_bar(scaled_nirK_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP nirK")
scaled_denite_ETSP_bulk_nirK$data$Sample <- factor(scaled_denite_ETSP_bulk_nirK$data$Sample, levels = desired_ETSP_bulk_order)
scaled_denite_ETSP_bulk_nirK <- (scaled_denite_ETSP_bulk_nirK) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETSP_bulk_nirK <- (scaled_denite_ETSP_bulk_nirK) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETSP_bulk_nirK <- (scaled_denite_ETSP_bulk_nirK) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETSP_bulk_nirK)

nirS_ETSP_bulk <- subset_taxa(ETSP_bulk1, nirS_val=="1")
scaled_nirS_ETSP_bulk = transform_sample_counts(nirS_ETSP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETSP_bulk_nirS <- plot_bar(scaled_nirS_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP nirS")
scaled_denite_ETSP_bulk_nirS$data$Sample <- factor(scaled_denite_ETSP_bulk_nirS$data$Sample, levels = desired_ETSP_bulk_order)
scaled_denite_ETSP_bulk_nirS <- (scaled_denite_ETSP_bulk_nirS) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETSP_bulk_nirS <- (scaled_denite_ETSP_bulk_nirS) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETSP_bulk_nirS <- (scaled_denite_ETSP_bulk_nirS) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETSP_bulk_nirS)

norB_ETSP_bulk <- subset_taxa(ETSP_bulk1, nor_new=="1")
scaled_norB_ETSP_bulk = transform_sample_counts(norB_ETSP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETSP_bulk_norB <- plot_bar(scaled_norB_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP norB")
scaled_denite_ETSP_bulk_norB$data$Sample <- factor(scaled_denite_ETSP_bulk_norB$data$Sample, levels = desired_ETSP_bulk_order)
scaled_denite_ETSP_bulk_norB <- (scaled_denite_ETSP_bulk_norB) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
scaled_denite_ETSP_bulk_norB <- (scaled_denite_ETSP_bulk_norB) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),axis.title.y=element_blank())
#scaled_denite_ETSP_bulk_norB <- (scaled_denite_ETSP_bulk_norB) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETSP_bulk_norB)

nosZ_ETSP_bulk <- subset_taxa(ETSP_bulk1, nosZ_val=="1")
scaled_nosZ_ETSP_bulk = transform_sample_counts(nosZ_ETSP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETSP_bulk_nosZ <- plot_bar(scaled_nosZ_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP nosZ")
scaled_denite_ETSP_bulk_nosZ$data$Sample <- factor(scaled_denite_ETSP_bulk_nosZ$data$Sample, levels = desired_ETSP_bulk_order)
scaled_denite_ETSP_bulk_nosZ <- (scaled_denite_ETSP_bulk_nosZ) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETSP_bulk_nosZ <- (scaled_denite_ETSP_bulk_nosZ) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETSP_bulk_nosZ <- (scaled_denite_ETSP_bulk_nosZ) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETSP_bulk_nosZ)

plist = list(scaled_denite_ETSP_bulk_napA, scaled_denite_ETSP_bulk_narG, scaled_denite_ETSP_bulk_nirK, scaled_denite_ETSP_bulk_nirS,scaled_denite_ETSP_bulk_norB, scaled_denite_ETSP_bulk_nosZ)

plist[1:length(plist)] = plist[1:length(plist)] %>% 
  map(~ .x + theme(axis.text.y=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   plot.margin=margin(l=0, unit="cm"))) 

scaled_denite_ETSP_all <- plot_grid(plotlist=plist, nrow=3,ncol=2, align="hv")
print(scaled_denite_ETSP_all)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="Denitrifier_ETSP_CoverM_genes_chemot_scaled.svg", plot=scaled_denite_ETSP_all)

nxrA_ETSP_bulk <- subset_taxa(ETSP_bulk1, nxrA=="1")
scaled_nxrA_ETSP_bulk = transform_sample_counts(nxrA_ETSP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETSP_bulk_nxrA <- plot_bar(scaled_nxrA_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP nxrA")
scaled_denite_ETSP_bulk_nxrA$data$Sample <- factor(scaled_denite_ETSP_bulk_nxrA$data$Sample, levels = desired_ETSP_bulk_order)
scaled_denite_ETSP_bulk_nxrA <- (scaled_denite_ETSP_bulk_nxrA) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETSP_bulk_nxrA <- (scaled_denite_ETSP_bulk_nxrA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETSP_bulk_nxrA <- (scaled_denite_ETSP_bulk_nxrA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETSP_bulk_nxrA)

nrfA_ETSP_bulk <- subset_taxa(ETSP_bulk1, nrfA=="1")
scaled_nrfA_ETSP_bulk = transform_sample_counts(nrfA_ETSP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETSP_bulk_nrfA <- plot_bar(scaled_nrfA_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP nrfA")
scaled_denite_ETSP_bulk_nrfA$data$Sample <- factor(scaled_denite_ETSP_bulk_nrfA$data$Sample, levels = desired_ETSP_bulk_order)
scaled_denite_ETSP_bulk_nrfA <- (scaled_denite_ETSP_bulk_nrfA) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETSP_bulk_nrfA <- (scaled_denite_ETSP_bulk_nrfA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETSP_bulk_nrfA <- (scaled_denite_ETSP_bulk_nrfA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETSP_bulk_nrfA)

nifH_ETSP_bulk <- subset_taxa(ETSP_bulk1, nifH=="1")
scaled_nifH_ETSP_bulk = transform_sample_counts(nifH_ETSP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETSP_bulk_nifH <- plot_bar(scaled_nifH_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP nifH")
scaled_denite_ETSP_bulk_nifH$data$Sample <- factor(scaled_denite_ETSP_bulk_nifH$data$Sample, levels = desired_ETSP_bulk_order)
scaled_denite_ETSP_bulk_nifH <- (scaled_denite_ETSP_bulk_nifH) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETSP_bulk_nifH <- (scaled_denite_ETSP_bulk_nifH) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETSP_bulk_nifH <- (scaled_denite_ETSP_bulk_nifH) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETSP_bulk_nifH)

amoA_ETSP_bulk <- subset_taxa(ETSP_bulk1, amoA=="1")
scaled_amoA_ETSP_bulk = transform_sample_counts(amoA_ETSP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETSP_bulk_amoA <- plot_bar(scaled_amoA_ETSP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP amoA")
scaled_denite_ETSP_bulk_amoA$data$Sample <- factor(scaled_denite_ETSP_bulk_amoA$data$Sample, levels = desired_ETSP_bulk_order)
scaled_denite_ETSP_bulk_amoA <- (scaled_denite_ETSP_bulk_amoA) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETSP_bulk_amoA <- (scaled_denite_ETSP_bulk_amoA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETSP_bulk_amoA <- (scaled_denite_ETSP_bulk_amoA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETSP_bulk_amoA)

hzo_ETSP_bulk <- subset_taxa(ETSP_bulk1, hzo=="1")
scaled_hzo_ETSP_bulk = transform_sample_counts(hzo_ETSP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETSP_bulk_hzo <- plot_bar(scaled_hzo_ETSP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP hzo")
scaled_denite_ETSP_bulk_hzo$data$Sample <- factor(scaled_denite_ETSP_bulk_hzo$data$Sample, levels = desired_ETSP_bulk_order)
scaled_denite_ETSP_bulk_hzo <- (scaled_denite_ETSP_bulk_hzo) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETSP_bulk_hzo <- (scaled_denite_ETSP_bulk_hzo) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETSP_bulk_hzo <- (scaled_denite_ETSP_bulk_hzo) +theme(legend.position = "none")+theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETSP_bulk_hzo)

nod_ETSP_bulk <- subset_taxa(ETSP_bulk1, nod=="1")
scaled_nod_ETSP_bulk = transform_sample_counts(nod_ETSP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETSP_bulk_nod <- plot_bar(scaled_nod_ETSP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETSP_bulk) + xlab("ETSP nod")
scaled_denite_ETSP_bulk_nod$data$Sample <- factor(scaled_denite_ETSP_bulk_nod$data$Sample, levels = desired_ETSP_bulk_order)
scaled_denite_ETSP_bulk_nod <- (scaled_denite_ETSP_bulk_nod) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETSP_bulk_nod <- (scaled_denite_ETSP_bulk_nod) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETSP_bulk_nod <- (scaled_denite_ETSP_bulk_nod) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETSP_bulk_nod)

otherlist = list(scaled_denite_ETSP_bulk_nxrA, scaled_denite_ETSP_bulk_nrfA, scaled_denite_ETSP_bulk_nifH,scaled_denite_ETSP_bulk_amoA, scaled_denite_ETSP_bulk_hzo, scaled_denite_ETSP_bulk_nod)

otherlist[1:length(otherlist)] = otherlist[1:length(otherlist)] %>% 
  map(~ .x + theme(axis.text.y=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   plot.margin=margin(l=0, unit="cm"))) 

scaled_othergenes_ETSP_all <- plot_grid(plotlist=otherlist, nrow=3,ncol=2, align="hv")
print(scaled_othergenes_ETSP_all)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="Denitrifier_CoverM_othergenes_chemot_scaled.svg", plot=scaled_othergenes_ETSP_all)

#Make taxa barplots for each gene individually for the ETNP
napA_ETNP_bulk <- subset_taxa(ETNP_bulk1, napA_val=="1")
scaled_napA_ETNP_bulk = transform_sample_counts(napA_ETNP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETNP_bulk_napA <- plot_bar(scaled_napA_ETNP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP napA")
scaled_denite_ETNP_bulk_napA$data$Sample <- factor(scaled_denite_ETNP_bulk_napA$data$Sample, levels = desired_ETNP_bulk_order)
scaled_denite_ETNP_bulk_napA <- (scaled_denite_ETNP_bulk_napA) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETNP_bulk_napA <- (scaled_denite_ETNP_bulk_napA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETNP_bulk_napA <- (scaled_denite_ETNP_bulk_napA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETNP_bulk_napA)

narG_ETNP_bulk <- subset_taxa(ETNP_bulk1, narG_val=="1")
scaled_narG_ETNP_bulk = transform_sample_counts(narG_ETNP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETNP_bulk_narG <- plot_bar(scaled_narG_ETNP_bulk , fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP narG")
scaled_denite_ETNP_bulk_narG$data$Sample <- factor(scaled_denite_ETNP_bulk_narG$data$Sample, levels = desired_ETNP_bulk_order)
scaled_denite_ETNP_bulk_narG <- (scaled_denite_ETNP_bulk_narG) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETNP_bulk_narG <- (scaled_denite_ETNP_bulk_narG) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETNP_bulk_narG <- (scaled_denite_ETNP_bulk_narG) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETNP_bulk_narG)

nirK_ETNP_bulk <- subset_taxa(ETNP_bulk1, nirK_val=="1")
scaled_nirK_ETNP_bulk = transform_sample_counts(nirK_ETNP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETNP_bulk_nirK <- plot_bar(scaled_nirK_ETNP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP nirK")
scaled_denite_ETNP_bulk_nirK$data$Sample <- factor(scaled_denite_ETNP_bulk_nirK$data$Sample, levels = desired_ETNP_bulk_order)
scaled_denite_ETNP_bulk_nirK <- (scaled_denite_ETNP_bulk_nirK) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETNP_bulk_nirK <- (scaled_denite_ETNP_bulk_nirK) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETNP_bulk_nirK <- (scaled_denite_ETNP_bulk_nirK) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETNP_bulk_nirK)

nirS_ETNP_bulk <- subset_taxa(ETNP_bulk1, nirS_val=="1")
scaled_nirS_ETNP_bulk = transform_sample_counts(nirS_ETNP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETNP_bulk_nirS <- plot_bar(scaled_nirS_ETNP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP nirS")
scaled_denite_ETNP_bulk_nirS$data$Sample <- factor(scaled_denite_ETNP_bulk_nirS$data$Sample, levels = desired_ETNP_bulk_order)
scaled_denite_ETNP_bulk_nirS <- (scaled_denite_ETNP_bulk_nirS) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETNP_bulk_nirS <- (scaled_denite_ETNP_bulk_nirS) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETNP_bulk_nirS <- (scaled_denite_ETNP_bulk_nirS) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETNP_bulk_nirS)

norB_ETNP_bulk <- subset_taxa(ETNP_bulk1, nor_new=="1")
scaled_norB_ETNP_bulk = transform_sample_counts(norB_ETNP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETNP_bulk_norB <- plot_bar(scaled_norB_ETNP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP norB")
scaled_denite_ETNP_bulk_norB$data$Sample <- factor(scaled_denite_ETNP_bulk_norB$data$Sample, levels = desired_ETNP_bulk_order)
scaled_denite_ETNP_bulk_norB <- (scaled_denite_ETNP_bulk_norB) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
scaled_denite_ETNP_bulk_norB <- (scaled_denite_ETNP_bulk_norB) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),axis.title.y=element_blank())
#scaled_denite_ETNP_bulk_norB <- (scaled_denite_ETNP_bulk_norB) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETNP_bulk_norB)

nosZ_ETNP_bulk <- subset_taxa(ETNP_bulk1, nosZ_val=="1")
scaled_nosZ_ETNP_bulk = transform_sample_counts(nosZ_ETNP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETNP_bulk_nosZ <- plot_bar(scaled_nosZ_ETNP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP nosZ")
scaled_denite_ETNP_bulk_nosZ$data$Sample <- factor(scaled_denite_ETNP_bulk_nosZ$data$Sample, levels = desired_ETNP_bulk_order)
scaled_denite_ETNP_bulk_nosZ <- (scaled_denite_ETNP_bulk_nosZ) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETNP_bulk_nosZ <- (scaled_denite_ETNP_bulk_nosZ) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETNP_bulk_nosZ <- (scaled_denite_ETNP_bulk_nosZ) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETNP_bulk_nosZ)

qlist = list(scaled_denite_ETNP_bulk_napA, scaled_denite_ETNP_bulk_narG, scaled_denite_ETNP_bulk_nirK, scaled_denite_ETNP_bulk_nirS,scaled_denite_ETNP_bulk_norB, scaled_denite_ETNP_bulk_nosZ)

qlist[1:length(qlist)] = qlist[1:length(qlist)] %>% 
  map(~ .x + theme(axis.text.y=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   plot.margin=margin(l=0, unit="cm"))) 

scaled_denite_ETNP_all <- plot_grid(plotlist=qlist, nrow=3,ncol=2, align="hv")
print(scaled_denite_ETNP_all)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="Denitrifier_ETNP_CoverM_genes_chemot_scaled.svg", plot=scaled_denite_ETNP_all)

nxrA_ETNP_bulk <- subset_taxa(ETNP_bulk1, nxrA=="1")
scaled_nxrA_ETNP_bulk = transform_sample_counts(nxrA_ETNP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETNP_bulk_nxrA <- plot_bar(scaled_nxrA_ETNP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP nxrA")
scaled_denite_ETNP_bulk_nxrA$data$Sample <- factor(scaled_denite_ETNP_bulk_nxrA$data$Sample, levels = desired_ETNP_bulk_order)
scaled_denite_ETNP_bulk_nxrA <- (scaled_denite_ETNP_bulk_nxrA) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETNP_bulk_nxrA <- (scaled_denite_ETNP_bulk_nxrA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETNP_bulk_nxrA <- (scaled_denite_ETNP_bulk_nxrA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETNP_bulk_nxrA)

nrfA_ETNP_bulk <- subset_taxa(ETNP_bulk1, nrfA=="1")
scaled_nrfA_ETNP_bulk = transform_sample_counts(nrfA_ETNP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETNP_bulk_nrfA <- plot_bar(scaled_nrfA_ETNP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP nrfA")
scaled_denite_ETNP_bulk_nrfA$data$Sample <- factor(scaled_denite_ETNP_bulk_nrfA$data$Sample, levels = desired_ETNP_bulk_order)
scaled_denite_ETNP_bulk_nrfA <- (scaled_denite_ETNP_bulk_nrfA) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETNP_bulk_nrfA <- (scaled_denite_ETNP_bulk_nrfA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETNP_bulk_nrfA <- (scaled_denite_ETNP_bulk_nrfA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETNP_bulk_nrfA)

nifH_ETNP_bulk <- subset_taxa(ETNP_bulk1, nifH=="1")
scaled_nifH_ETNP_bulk = transform_sample_counts(nifH_ETNP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETNP_bulk_nifH <- plot_bar(scaled_nifH_ETNP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP nifH")
scaled_denite_ETNP_bulk_nifH$data$Sample <- factor(scaled_denite_ETNP_bulk_nifH$data$Sample, levels = desired_ETNP_bulk_order)
scaled_denite_ETNP_bulk_nifH <- (scaled_denite_ETNP_bulk_nifH) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETNP_bulk_nifH <- (scaled_denite_ETNP_bulk_nifH) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETNP_bulk_nifH <- (scaled_denite_ETNP_bulk_nifH) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETNP_bulk_nifH)

amoA_ETNP_bulk <- subset_taxa(ETNP_bulk1, amoA=="1")
scaled_amoA_ETNP_bulk = transform_sample_counts(amoA_ETNP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETNP_bulk_amoA <- plot_bar(scaled_amoA_ETNP_bulk, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP amoA")
scaled_denite_ETNP_bulk_amoA$data$Sample <- factor(scaled_denite_ETNP_bulk_amoA$data$Sample, levels = desired_ETNP_bulk_order)
scaled_denite_ETNP_bulk_amoA <- (scaled_denite_ETNP_bulk_amoA) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETNP_bulk_amoA <- (scaled_denite_ETNP_bulk_amoA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETNP_bulk_amoA <- (scaled_denite_ETNP_bulk_amoA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETNP_bulk_amoA)

hzo_ETNP_bulk <- subset_taxa(ETNP_bulk1, hzo=="1")
scaled_hzo_ETNP_bulk = transform_sample_counts(hzo_ETNP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETNP_bulk_hzo <- plot_bar(scaled_hzo_ETNP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP hzo")
scaled_denite_ETNP_bulk_hzo$data$Sample <- factor(scaled_denite_ETNP_bulk_hzo$data$Sample, levels = desired_ETNP_bulk_order)
scaled_denite_ETNP_bulk_hzo <- (scaled_denite_ETNP_bulk_hzo) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETNP_bulk_hzo <- (scaled_denite_ETNP_bulk_hzo) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETNP_bulk_hzo <- (scaled_denite_ETNP_bulk_hzo) +theme(legend.position = "none")+theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETNP_bulk_hzo)

nod_ETNP_bulk <- subset_taxa(ETNP_bulk1, nod=="1")
scaled_nod_ETNP_bulk = transform_sample_counts(nod_ETNP_bulk, function(x) (x / sum(x))*100 )
scaled_denite_ETNP_bulk_nod <- plot_bar(scaled_nod_ETNP_bulk, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_ETNP_bulk) + xlab("ETNP nod")
scaled_denite_ETNP_bulk_nod$data$Sample <- factor(scaled_denite_ETNP_bulk_nod$data$Sample, levels = desired_ETNP_bulk_order)
scaled_denite_ETNP_bulk_nod <- (scaled_denite_ETNP_bulk_nod) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_ETNP_bulk_nod <- (scaled_denite_ETNP_bulk_nod) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_ETNP_bulk_nod <- (scaled_denite_ETNP_bulk_nod) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_ETNP_bulk_nod)

rlist = list(scaled_denite_ETNP_bulk_nxrA, scaled_denite_ETNP_bulk_nrfA, scaled_denite_ETNP_bulk_nifH,scaled_denite_ETNP_bulk_amoA, scaled_denite_ETNP_bulk_hzo, scaled_denite_ETNP_bulk_nod)

rlist[1:length(rlist)] = rlist[1:length(rlist)] %>% 
  map(~ .x + theme(axis.text.y=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   plot.margin=margin(l=0, unit="cm"))) 

scaled_othergenes_ETNP_all <- plot_grid(plotlist=rlist, nrow=3,ncol=2, align="hv")
print(scaled_othergenes_ETNP_all)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="Denitrifier_ETNP_CoverM_othergenes_chemot_scaled.svg", plot=scaled_othergenes_ETNP_all)


#Make motility and chemotaxis barplots for each gene individually for the Arabian Sea
napA_Arabian <- subset_taxa(Arabian1, napA_val=="1")
scaled_napA_Arabian = transform_sample_counts(napA_Arabian, function(x) (x / sum(x))*100 )
scaled_denite_Arabian_napA <- plot_bar(scaled_napA_Arabian, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian napA")
scaled_denite_Arabian_napA$data$Sample <- factor(scaled_denite_Arabian_napA$data$Sample, levels = desired_AS_order)
scaled_denite_Arabian_napA <- (scaled_denite_Arabian_napA) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_Arabian_napA <- (scaled_denite_Arabian_napA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_Arabian_napA <- (scaled_denite_Arabian_napA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_Arabian_napA)

narG_Arabian <- subset_taxa(Arabian1, narG_val=="1")
scaled_narG_Arabian = transform_sample_counts(narG_Arabian, function(x) (x / sum(x))*100 )
scaled_denite_Arabian_narG <- plot_bar(scaled_narG_Arabian , fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian narG")
scaled_denite_Arabian_narG$data$Sample <- factor(scaled_denite_Arabian_narG$data$Sample, levels = desired_AS_order)
scaled_denite_Arabian_narG <- (scaled_denite_Arabian_narG) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_Arabian_narG <- (scaled_denite_Arabian_narG) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_Arabian_narG <- (scaled_denite_Arabian_narG) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_Arabian_narG)

nirK_Arabian <- subset_taxa(Arabian1, nirK_val=="1")
scaled_nirK_Arabian = transform_sample_counts(nirK_Arabian, function(x) (x / sum(x))*100 )
scaled_denite_Arabian_nirK <- plot_bar(scaled_nirK_Arabian, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian nirK")
scaled_denite_Arabian_nirK$data$Sample <- factor(scaled_denite_Arabian_nirK$data$Sample, levels = desired_AS_order)
scaled_denite_Arabian_nirK <- (scaled_denite_Arabian_nirK) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_Arabian_nirK <- (scaled_denite_Arabian_nirK) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_Arabian_nirK <- (scaled_denite_Arabian_nirK) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_Arabian_nirK)

nirS_Arabian <- subset_taxa(Arabian1, nirS_val=="1")
scaled_nirS_Arabian = transform_sample_counts(nirS_Arabian, function(x) (x / sum(x))*100 )
scaled_denite_Arabian_nirS <- plot_bar(scaled_nirS_Arabian, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian nirS")
scaled_denite_Arabian_nirS$data$Sample <- factor(scaled_denite_Arabian_nirS$data$Sample, levels = desired_AS_order)
scaled_denite_Arabian_nirS <- (scaled_denite_Arabian_nirS) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_Arabian_nirS <- (scaled_denite_Arabian_nirS) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_Arabian_nirS <- (scaled_denite_Arabian_nirS) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_Arabian_nirS)

norB_Arabian <- subset_taxa(Arabian1, nor_new=="1")
scaled_norB_Arabian = transform_sample_counts(norB_Arabian, function(x) (x / sum(x))*100 )
scaled_denite_Arabian_norB <- plot_bar(scaled_norB_Arabian, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian norB")
scaled_denite_Arabian_norB$data$Sample <- factor(scaled_denite_Arabian_norB$data$Sample, levels = desired_AS_order)
scaled_denite_Arabian_norB <- (scaled_denite_Arabian_norB) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
scaled_denite_Arabian_norB <- (scaled_denite_Arabian_norB) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),axis.title.y=element_blank())
#scaled_denite_Arabian_norB <- (scaled_denite_Arabian_norB) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_Arabian_norB)

nosZ_Arabian <- subset_taxa(Arabian1, nosZ_val=="1")
scaled_nosZ_Arabian = transform_sample_counts(nosZ_Arabian, function(x) (x / sum(x))*100 )
scaled_denite_Arabian_nosZ <- plot_bar(scaled_nosZ_Arabian, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian nosZ")
scaled_denite_Arabian_nosZ$data$Sample <- factor(scaled_denite_Arabian_nosZ$data$Sample, levels = desired_AS_order)
scaled_denite_Arabian_nosZ <- (scaled_denite_Arabian_nosZ) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_Arabian_nosZ <- (scaled_denite_Arabian_nosZ) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_Arabian_nosZ <- (scaled_denite_Arabian_nosZ) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_Arabian_nosZ)

slist = list(scaled_denite_Arabian_napA, scaled_denite_Arabian_narG, scaled_denite_Arabian_nirK, scaled_denite_Arabian_nirS,scaled_denite_Arabian_norB, scaled_denite_Arabian_nosZ)

slist[1:length(slist)] = slist[1:length(slist)] %>% 
  map(~ .x + theme(axis.text.y=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   plot.margin=margin(l=0, unit="cm"))) 

scaled_denite_Arabian_all <- plot_grid(plotlist=slist, nrow=3,ncol=2, align="hv")
print(scaled_denite_Arabian_all)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="Denitrifier_Arabian_CoverM_genes_chemot_scaled.svg", plot=scaled_denite_Arabian_all)

nxrA_Arabian <- subset_taxa(Arabian1, nxrA=="1")
scaled_nxrA_Arabian = transform_sample_counts(nxrA_Arabian, function(x) (x / sum(x))*100 )
scaled_denite_Arabian_nxrA <- plot_bar(scaled_nxrA_Arabian, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian nxrA")
scaled_denite_Arabian_nxrA$data$Sample <- factor(scaled_denite_Arabian_nxrA$data$Sample, levels = desired_AS_order)
scaled_denite_Arabian_nxrA <- (scaled_denite_Arabian_nxrA) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_Arabian_nxrA <- (scaled_denite_Arabian_nxrA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_Arabian_nxrA <- (scaled_denite_Arabian_nxrA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_Arabian_nxrA)

nrfA_Arabian <- subset_taxa(Arabian1, nrfA=="1")
scaled_nrfA_Arabian = transform_sample_counts(nrfA_Arabian, function(x) (x / sum(x))*100 )
scaled_denite_Arabian_nrfA <- plot_bar(scaled_nrfA_Arabian, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian nrfA")
scaled_denite_Arabian_nrfA$data$Sample <- factor(scaled_denite_Arabian_nrfA$data$Sample, levels = desired_AS_order)
scaled_denite_Arabian_nrfA <- (scaled_denite_Arabian_nrfA) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_Arabian_nrfA <- (scaled_denite_Arabian_nrfA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_Arabian_nrfA <- (scaled_denite_Arabian_nrfA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_Arabian_nrfA)

nifH_Arabian <- subset_taxa(Arabian1, nifH=="1")
scaled_nifH_Arabian = transform_sample_counts(nifH_Arabian, function(x) (x / sum(x))*100 )
scaled_denite_Arabian_nifH <- plot_bar(scaled_nifH_Arabian, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian nifH")
scaled_denite_Arabian_nifH$data$Sample <- factor(scaled_denite_Arabian_nifH$data$Sample, levels = desired_AS_order)
scaled_denite_Arabian_nifH <- (scaled_denite_Arabian_nifH) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_Arabian_nifH <- (scaled_denite_Arabian_nifH) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_Arabian_nifH <- (scaled_denite_Arabian_nifH) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_Arabian_nifH)

amoA_Arabian <- subset_taxa(Arabian1, amoA=="1")
scaled_amoA_Arabian = transform_sample_counts(amoA_Arabian, function(x) (x / sum(x))*100 )
scaled_denite_Arabian_amoA <- plot_bar(scaled_amoA_Arabian, fill = "chem_mot_aer") +  scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian amoA")
scaled_denite_Arabian_amoA$data$Sample <- factor(scaled_denite_Arabian_amoA$data$Sample, levels = desired_AS_order)
scaled_denite_Arabian_amoA <- (scaled_denite_Arabian_amoA) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_Arabian_amoA <- (scaled_denite_Arabian_amoA) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_Arabian_amoA <- (scaled_denite_Arabian_amoA) + theme(legend.position = "none") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_Arabian_amoA)

hzo_Arabian <- subset_taxa(Arabian1, hzo=="1")
scaled_hzo_Arabian = transform_sample_counts(hzo_Arabian, function(x) (x / sum(x))*100 )
scaled_denite_Arabian_hzo <- plot_bar(scaled_hzo_Arabian, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian hzo")
scaled_denite_Arabian_hzo$data$Sample <- factor(scaled_denite_Arabian_hzo$data$Sample, levels = desired_AS_order)
scaled_denite_Arabian_hzo <- (scaled_denite_Arabian_hzo) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_Arabian_hzo <- (scaled_denite_Arabian_hzo) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_Arabian_hzo <- (scaled_denite_Arabian_hzo) +theme(legend.position = "none")+theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_Arabian_hzo)

nod_Arabian <- subset_taxa(Arabian1, nod=="1")
scaled_nod_Arabian = transform_sample_counts(nod_Arabian, function(x) (x / sum(x))*100 )
scaled_denite_Arabian_nod <- plot_bar(scaled_nod_Arabian, fill = "chem_mot_aer") + scale_fill_manual(values=P10)+
  geom_bar( stat="identity", position="stack")+ 
  scale_x_discrete(labels = relabel_AS) + xlab("Arabian nod")
scaled_denite_Arabian_nod$data$Sample <- factor(scaled_denite_Arabian_nod$data$Sample, levels = desired_AS_order)
scaled_denite_Arabian_nod <- (scaled_denite_Arabian_nod) + coord_cartesian(ylim = c(0, 102))+ theme(text = element_text(size = 7))
#scaled_denite_Arabian_nod <- (scaled_denite_Arabian_nod) + theme(legend.position = "right") + theme(axis.text.x=element_text(angle=45,hjust=1)) 
scaled_denite_Arabian_nod <- (scaled_denite_Arabian_nod) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),axis.title.y=element_blank())
print(scaled_denite_Arabian_nod)

tlist = list(scaled_denite_Arabian_nxrA, scaled_denite_Arabian_nrfA, scaled_denite_Arabian_nifH,scaled_denite_Arabian_amoA, scaled_denite_Arabian_hzo, scaled_denite_Arabian_nod)

tlist[1:length(tlist)] = tlist[1:length(tlist)] %>% 
  map(~ .x + theme(axis.text.y=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   plot.margin=margin(l=0, unit="cm"))) 

scaled_othergenes_Arabian_all <- plot_grid(plotlist=tlist, nrow=3,ncol=2, align="hv")
print(scaled_othergenes_Arabian_all)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="Denitrifier_Arabian_CoverM_othergenes_chemot_scaled.svg", plot=scaled_othergenes_Arabian_all)


