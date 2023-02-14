library(eulerr)
library(dplyr)
library(xlsx)
library(phyloseq)
library(cowplot)
library(ggplot2)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets")

#need to read in abundance tables
#added columns for denitrification steps into abundance table taxonomy information 
#make sure it's Sheet 7 with combined napnar and nir information

#read in all MAG tables from CoverM
MAG_abund_table_raw <- read.xlsx2("CoverM_GTDB_joined.xlsx", 
                                  sheetIndex=1,row.names = 1)
MAG_abund_table_raw <- mutate_all(MAG_abund_table_raw, function(x) as.numeric(as.character(x)))
MAG_matrix <- as.matrix(MAG_abund_table_raw)

#read in tax table
tax_mat_val <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=7, row.names=1)
tax_mat_val <- as.matrix(tax_mat_val)
tax_mat_val = subset(tax_mat_val, select = -c(1:9,11:12,14,16:18,20:27) )

meta <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=3, row.names=1)

OTU = otu_table(MAG_matrix, taxa_are_rows = TRUE)
TAX2 = tax_table(tax_mat_val)
samples = sample_data(meta)

#make phyloseq object
MAG_phyloseq_val = phyloseq(OTU, TAX2, samples)
rank_names(MAG_phyloseq_val)
denitrifier1 <- subset_taxa(MAG_phyloseq_val, denitrifier_val=="yes")

#add up relative abundances for all MAGs in each of these sets for each sample
glomtest <- tax_glom(denitrifier1 , taxrank="specialization_val", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
melttest<-psmelt(glomtest)

#view and export as a text file
#edited outside of R to check accuracy and remove unneeded columns and info
View(melttest)
str(melttest)
write.table(as.data.frame(melttest), file = "melttest.txt", sep = "\t")

#plot euler diagrams with this information (saved as Excel)
euler_info <- read.xlsx2("melttest.xlsx", sheetIndex=1)

#subset out for each sample
euler_Amal1 <- subset(euler_info, Sample=="AMAL1")
euler_Amal1 = subset(euler_Amal1, select = -c(1,7) )
euler_Amal1[,1]<-as.numeric(euler_Amal1[,1])

samples <- unique(euler_info[,1])

subset_samples <- split(euler_info, as.factor(euler_info$Sample))
subset_genes <- split(euler_info, as.factor(euler_info$specialization_val))
sample_order <- subset_genes$`123`[,1]
genes_order <- names(subset_genes)

euler_Amal1[,1]<-as.numeric(euler_Amal1[,1])
df1 <-lapply(subset_samples[[1]][,2],as.numeric)


#look at how to index a list of dataframes: df[[listindex]][x,y]
#change types to numeric
for (sample in 1:length(sample_order)) {
  subset_samples[[sample]][,2]<-as.numeric(subset_samples[[sample]][,2])
}

#initialize empty list
euler_list <- list()
euler_list

#for loop to apply euler function across all samples
#this is a good for loop but unfortunately it doesn't output the right format for euler...
for (sample in 1:length(sample_order)) {
  tmp <- euler(c("nap_nar" = subset_samples[[sample]][1,2], 
        "nir" = subset_samples[[sample]][2,2],
        "nor" = subset_samples[[sample]][3,2],
        "nos" = subset_samples[[sample]][4,2], 
        "nap_nar&nir" = subset_samples[[sample]][5,2], 
        "nap_nar&nor" = subset_samples[[sample]][6,2],
        "nap_nar&nos" = subset_samples[[sample]][7,2],
        "nir&nor" = subset_samples[[sample]][8,2],
        "nir&nos" = subset_samples[[sample]][9,2],
        "nap_nar&nir&nor" = subset_samples[[sample]][10,2],
        "nap_nar&nir&nos" = subset_samples[[sample]][11,2],
        "nap_nar&nor&nos" = subset_samples[[sample]][12,2],
        "nir&nor&nos" = subset_samples[[sample]][13,2]),
        shape="ellipse",factor_names=TRUE
        
  )
  euler_list <- rbind(euler_list, tmp)
}

#changed name to make it easier to find and replace
#find and replace to run this 56 times ugh
splits <- subset_samples
Stewart_SRR070082 <-euler(c("nap_nar" = splits[[44]][1,2], 
        "nir" = splits[[44]][2,2],
        "nor" = splits[[44]][3,2],
        "nos" = splits[[44]][4,2], 
        "nap_nar&nir" = splits[[44]][5,2], 
        "nap_nar&nor" = splits[[44]][6,2],
        "nap_nar&nos" = splits[[44]][7,2],
        "nir&nor" = splits[[44]][8,2],
        "nir&nos" = splits[[44]][9,2],
        "nap_nar&nir&nor" = splits[[44]][10,2],
        "nap_nar&nir&nos" = splits[[44]][11,2],
        "nap_nar&nor&nos" = splits[[44]][12,2],
        "nir&nor&nos" = splits[[44]][13,2]),
      shape="ellipse",factor_names=TRUE)

plot_Glass_SRR1509799 <- plot(
  Glass_SRR1509799,
  fills = TRUE,
  edges = TRUE,
  legend = TRUE,
  labels = identical(legend, TRUE),
  quantities = FALSE,
  strips = NULL,
  main = "ETNP_Gl_2013_125m",
  n = 200L
)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")

all_surface <- plot_grid(plot_Amal18, plot_Amal1, 
                         plot_Amal13,plot_Glass_SRR1509797,plot_Stewart_SRR304671,
                         align = 'hv',
                         nrow =2,
                         ncol=3)
all_surface
ggsave(filename="surface_eulers.svg", plot=all_surface)
ggsave(filename="surface_eulers.png", plot=all_surface)

ETSP_surface <- plot_grid(plot_Stewart_SRR304684, plot_Stewart_SRR304671, 
                          align = 'hv',
                          nrow =1,
                          ncol=2)
ETSP_surface

ETNP_surface <- plot_grid(plot_Amal5, plot_Amal17, 
                          plot_Glass_SRR1509798,plot_Glass_SRR1509797,
                           align = 'hv',
                           nrow =2,
                           ncol=2)
ETNP_surface

ETNP_oxycline <- plot_grid(plot_Fuchsman_SRR4465037,
                          plot_Tsementzi_SRR3718412,plot_Fuchsman_SRR4465034,
                          plot_Glass_SRR1509798,plot_Glass_SRR1509792,
                          plot_Fuchsman_SRR4465025,plot_Amal14,
                          align = 'hv',
                          nrow =3,
                          ncol=3)
ETNP_oxycline
ggsave(filename="ETNP_oxycline_eulers.svg", plot=ETNP_oxycline)
ggsave(filename="ETNP_oxycline_eulers.png", plot=ETNP_oxycline)

ETNP_ODZ1 <- plot_grid(plot_Fuchsman_SRR4465024,
                      plot_Glass_SRR1509793,
                      plot_Fuchsman_SRR4465027,
                      plot_Amal2,
                      plot_Fuchsman_SRR4465026,
                      plot_Tsementzi_SRR3718413,
                      plot_Glass_SRR1509794,
                      plot_Glass_SRR1509799,
                      plot_Fuchsman_SRR4465029,
                      plot_Fuchsman_SRR4465028,
                      align = 'hv',
                      nrow =3,
                      ncol=4)
ETNP_ODZ1
ggsave(filename="ETNP_ODZ_shallower_eulers.svg", plot=ETNP_ODZ1)
ggsave(filename="ETNP_ODZ_shallower_eulers.png", plot=ETNP_ODZ1)

ETNP_ODZ2 <- plot_grid(plot_Fuchsman_SRR4465036,plot_Amal7,
                       plot_Amal3,plot_Amal15,plot_Amal8,
                       plot_Amal20,plot_Fuchsman_SRR4465035,
                       plot_Glass_SRR1509796,plot_Glass_SRR1509800,
                       align='hv',
                       nrow=3,
                       ncol=3)
ETNP_ODZ2
ggsave(filename="ETNP_ODZ_deeper_eulers.svg", plot=ETNP_ODZ2)
ggsave(filename="ETNP_ODZ_deeper_eulers.png", plot=ETNP_ODZ2)

ETSP_oxycline <- plot_grid(plot_Stewart_SRR064444,
                           plot_Stewart_SRR304672,
                           plot_Stewart_SRR304674,
                           plot_Stewart_SRR304656,
                           plot_Stewart_SRR070081,
                           plot_Ganesh_SRR960580,
                           plot_Stewart_SRR064448,
                           plot_Stewart_SRR304673,
                           plot_Stewart_SRR304680,
                           plot_Ganesh_SRR961673,
                           align = 'hv',
                           nrow =3,
                           ncol=4)
ETSP_oxycline
ggsave(filename="ETSP_oxycline_eulers.svg", plot=ETSP_oxycline)
ggsave(filename="ETSP_oxycline_eulers.png", plot=ETSP_oxycline)

ETSP_ODZ <- plot_grid(plot_Stewart_SRR070084,
                           plot_Stewart_SRR064450,
                           plot_Stewart_SRR070082,
                           plot_Ganesh_SRR961676,
                           plot_Stewart_SRR304668,
                           plot_Stewart_SRR304683,
                           plot_Ganesh_SRR961679,
                           align = 'hv',
                           nrow =3,
                           ncol=2)
ETSP_ODZ
ggsave(filename="ETSP_ODZ_eulers.svg", plot=ETSP_ODZ)
ggsave(filename="ETSP_ODZ_eulers.png", plot=ETSP_ODZ)

Arabian_plots <- plot_grid(plot_Amal9, plot_Amal10, plot_Amal11,
          plot_Amal12,
          align = 'hv',
          nrow =2,
          ncol=2)
Arabian_plots

ggsave(filename="Arabian_eulers.svg", plot=Arabian_plots)
ggsave(filename="Arabian_eulers.png", plot=Arabian_plots)

desired_AS_order <- list("AMAL9", "AMAL10", "AMAL11", "AMAL12")
relabel_AS <- list('2007_130m','2007_150m','2007_200m','2007_400m')
Arabian <- do.call(rbind, Map(data.frame, A=desired_AS_order, B=relabel_AS))
colnames(Arabian)[1] = "Sample"

desired_ETNP_bulk_order <- list("AMAL5","AMAL17","Glass_SRR1509790","Glass_SRR1509797","AMAL18","AMAL1","AMAL13",
                                "Fuchsman_SRR4465037","Tsementzi_SRR3718412","Fuchsman_SRR4465034","Glass_SRR1509798",
                                "Glass_SRR1509792","Fuchsman_SRR4465025","AMAL14","Fuchsman_SRR4465024",
                                "Glass_SRR1509793","Fuchsman_SRR4465027","AMAL2","Fuchsman_SRR4465026",
                                "Tsementzi_SRR3718413","Glass_SRR1509794","Glass_SRR1509799",
                                "Fuchsman_SRR4465029","Fuchsman_SRR4465028","Fuchsman_SRR4465036","AMAL7",
                                "AMAL3","AMAL15","AMAL8","AMAL20","Fuchsman_SRR4465035","Glass_SRR1509796","Glass_SRR1509800")
relabel_ETNP_bulk <- list('2016_10m','2018_16m','Gl_2013_30m','Gl_2013_30m','2018_45m','2016_53m','2018_60m','Fu_2012_60m','Ts_2013_68m','Fu_2012_70m',
                       'Gl_2013_80m','Gl_2013_85m','Fu_2012_90m','2018_95m','Fu_2012_100m','Gl_2013_100m','Fu_2012_110m','2016_120m','Fu_2012_120m',
                       'Ts_2013_120m','Gl_2013_125m','Gl_2013_125m','Fu_2012_140m','Fu_2012_160m','Fu_2012_180m','2016_185m',
                       '2016_200m','2018_200m','2016_215m','2018_250m','Fu_2012_300m','Gl_2013_300m','Gl_2013_300m')
ETNP <- do.call(rbind, Map(data.frame, A=desired_ETNP_bulk_order, B=relabel_ETNP_bulk))
colnames(ETNP)[1] = "Sample"

desired_ETSP_bulk_order <- list("Stewart_SRR304684","Stewart_SRR304671","Stewart_SRR064444","Stewart_SRR304672","Stewart_SRR304674",
                                "Stewart_SRR304656","Stewart_SRR070081","Ganesh_SRR960580","Stewart_SRR064448",
                                "Stewart_SRR304673","Stewart_SRR304680","Ganesh_SRR961673","Stewart_SRR070084",
                                "Stewart_SRR064450","Stewart_SRR070082","Ganesh_SRR961676","Stewart_SRR304668",
                                "Stewart_SRR304683","Ganesh_SRR961679")
relabel_ETSP_bulk <- list('St_2008_15m','St_2008_35m','St_2008_50m','St_2008_50m','St_2008_50m',
                       'St_2008_65m','St_2008_70m','Ga_2010_70m','St_2008_110m','St_2008_110m','St_2008_110m',
                       'Ga_2010_110m','St_2008_150m','St_2008_200m','St_2008_200m','Ga_2010_200m','St_2008_500m',
                       'St_2008_800m','Ga_2010_1000m')
ETSP <- do.call(rbind, Map(data.frame, A=desired_ETSP_bulk_order, B=relabel_ETSP_bulk))
colnames(ETSP)[1] = "Sample"

all_samples <- rbind(ETNP, ETSP, Arabian)

others <- plot_grid(plot_Amal18, plot_Amal1, plot_Amal11,
                    plot_Amal12,
                    align = 'hv',
                    nrow =1,
                    ncol=2)
others
setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="other_eulers.svg", plot=others)

euler_first<- subset(euler_info, specialization_val==1 | specialization_val==12| specialization_val==13| specialization_val==14| specialization_val==123| specialization_val==124| specialization_val==134)

euler_first$Abundance <- as.numeric(euler_first$Abundance)
euler_first$total <- euler_first %>%
  group_by(Sample) %>% 
  transmute(Total=sum(Abundance))

euler_first <- merge(
  euler_first, 
  all_samples, by="Sample", all = T) 

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets")
write.csv(euler_first, file = "euler_first.csv")

