library(readxl)
library(UpSetR)
library(xlsx)
library(ComplexUpset)
library(ggplot2)
library(pals)
library(Polychrome)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets") 

#read in hard-coded palette
palette_34 <- read.csv("P34.csv", header=TRUE)
P34 <- c(palette_34[,2])
PhylumList <- c(palette_34[,1])
names(P34) <- PhylumList

#read in gene information
all_MAGs <- read.xlsx("ODZ_MAG_all_info.xlsx", sheetName="all_ODZ_validated_newnors")%>%data.frame


#create a new column combining nirS and nirK and napA and narG, then change to binary
all_MAGs$nir_combined <- rowSums(all_MAGs[,c("nirS_val","nirK_val")], na.rm=F)
all_MAGs["nir_combined"][all_MAGs["nir_combined"]=="2"]<-"1"
all_MAGs["nir_combined"]<-as.numeric(all_MAGs$nir_combined)

all_MAGs$napnar_combined <- rowSums(all_MAGs[,c("napA_val","narG_val")], na.rm=F)
all_MAGs["napnar_combined"][all_MAGs["napnar_combined"]=="2"]<-"1"
all_MAGs["napnar_combined"]<-as.numeric(all_MAGs$napnar_combined)

#subset based on different things to plot
all_ETNP_MAGs <- subset(all_MAGs, Ocean == 'ETNP')%>%data.frame
all_Arabian_MAGs <- subset(all_MAGs, Ocean == 'Arabian Sea')%>%data.frame
all_dereplicated_MAGs <- subset(all_MAGs, dereplicated == 'dereplicated')%>%data.frame

all_MAGs_over_90 <- subset(all_MAGs,  completeness >90)%>%data.frame
all_MAGs_over_90_ETNP <- subset(all_MAGs, Ocean == 'ETNP' & completeness >90)
all_MAGs_over_90_Arabian <- subset(all_MAGs, Ocean == 'Arabian Sea' & completeness >90)

all_MAGs_surface <- subset(all_MAGs,Depth_bin == 'surface')%>%data.frame
all_MAGs_oxycline <- subset(all_MAGs,Depth_bin == 'oxycline')%>%data.frame
all_MAGs_ODZ <- subset(all_MAGs,Depth_bin == 'ODZ core')%>%data.frame

all_MAGs_over_90_ETNP_surface <- subset(all_MAGs, Ocean == 'ETNP' & completeness >90 & Depth_bin == 'surface')
all_MAGs_over_90_ETNP_oxycline <- subset(all_MAGs, Ocean == 'ETNP' & completeness >90 & Depth_bin == 'oxycline')
all_MAGs_over_90_ETNP_ODZ <- subset(all_MAGs, Ocean == 'ETNP' & completeness >90 & Depth_bin == 'ODZ core')

all_MAGs_over_90_derep <- subset(all_MAGs, completeness >90 & dereplicated == 'dereplicated' )
all_MAGs_over_90_ETNP_derep <- subset(all_MAGs, Ocean == 'ETNP' & completeness >90 & dereplicated == 'dereplicated' )
all_MAGs_over_90_Arabian_derep <- subset(all_MAGs, Ocean == 'Arabian Sea' & completeness >90 & dereplicated == 'dereplicated')

all_MAGs_except_surface <- subset(all_MAGs, Depth_bin == 'oxycline' | Depth_bin == 'ODZ core')
all_MAGs_except_surface_derep <- subset(all_MAGs, Depth_bin == 'oxycline' | Depth_bin == 'ODZ core' & dereplicated == 'dereplicated') 

#list genes to plot intersections of
genes = list("napA_val", "narG_val", "nirK_val", "nirS_val", "nor_new", "nosZ_val", "amoA","nifH","nrfA","nxrA")
as.character(genes)

denite_genes = list("napA_val", "narG_val", "nirK_val", "nirS_val", "nor_new", "nosZ_val")
as.character(denite_genes)

denite_steps = list("napnar_combined","nir_combined","nor_new","nosZ_val")
as.character(denite_steps)


MAGs_scaled_all <- upset(
  all_MAGs,
  genes,
  sort_intersections_by=c('degree'),
  sort_sets=FALSE,
  base_annotations = list(
    'Intersection size'=intersection_size(
      mapping=aes(fill=Phylum))
    + scale_fill_manual(values=P34)
    + ylab('Number of MAGs in Gene Set')
  ),
  annotations = list(
    'Completeness'=list(
      aes=aes(x=intersection, y=completeness),
      geom=list(
        geom_jitter(na.rm=TRUE, size = 0.2),
        geom_violin(alpha=0.5, na.rm=TRUE)
      )
  )
  ),
  min_size=5,
  width_ratio=0.2,
  min_degree=1,
  encode_sets=FALSE,
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=..count..), hjust=0.8, stat='count')
    + annotate(geom='text', label = '@', x='Number of MAGs', y=120, color='white', size=2)
    + theme(axis.text.x=element_text(angle=90))
  )
)+ patchwork::plot_layout(heights=c(0.5, 1.5, 1))


MAGs_scaled_denite_all <- upset(
  all_MAGs,
  denite_genes,
  sort_intersections_by=c('degree'),
  sort_sets=FALSE,
  base_annotations = list(
    'Intersection size'=intersection_size(
      mapping=aes(fill=Phylum))
    + scale_fill_manual(values=P34)
    + ylab('Number of MAGs in Gene Set')
  ),
  annotations = list(
    'Completeness'=list(
      aes=aes(x=intersection, y=completeness),
      geom=list(
        geom_jitter(na.rm=TRUE, size = 0.2),
        geom_violin(alpha=0.5, na.rm=TRUE)
      )
    )
  ),
  min_size=1,
  width_ratio=0.2,
  min_degree=1,
  encode_sets=FALSE,
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=..count..), hjust=0.8, stat='count')
    + annotate(geom='text', label = '@', x='Number of MAGs', y=120, color='white', size=2)
    + theme(axis.text.x=element_text(angle=90))
  )
)+ patchwork::plot_layout(heights=c(0.55, 1.5, 1))

MAGs_scaled_denite_steps <- upset(
  all_MAGs,
  denite_steps,
  sort_intersections_by=c('degree'),
  sort_sets=FALSE,
  base_annotations = list(
    'Intersection size'=intersection_size(
      mapping=aes(fill=Phylum))
    + scale_fill_manual(values=P34)
    + ylab('Number of MAGs in Gene Set')
  ),
  annotations = list(
    'Completeness'=list(
      aes=aes(x=intersection, y=completeness),
      geom=list(
        geom_jitter(na.rm=TRUE, size = 0.2),
        geom_violin(alpha=0.5, na.rm=TRUE)
      )
    )
  ),
  min_size=1,
  width_ratio=0.2,
  min_degree=1,
  encode_sets=FALSE,
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=..count..), hjust=0.8, stat='count')
    + annotate(geom='text', label = '@', x='Number of MAGs', y=120, color='white', size=2)
    + theme(axis.text.x=element_text(angle=90))
  )
)+ patchwork::plot_layout(heights=c(0.5, 1.5, 1))


MAGs_scaled_over_90 <- upset(
  all_MAGs_over_90,
  genes,
  sort_intersections_by=c('degree'),
  sort_sets=FALSE,
  base_annotations = list(
    'Intersection size'=intersection_size(
      mapping=aes(fill=Phylum))
    + scale_fill_manual(values=P34)
    + ylab('Number of MAGs in Gene Set')
  ),
  annotations = list(
    'Completeness'=list(
      aes=aes(x=intersection, y=completeness),
      geom=list(
        geom_jitter(na.rm=TRUE),
        geom_violin(alpha=0.5, na.rm=TRUE)
      )
    )
  ),
  min_size=1,
  width_ratio=0.1,
  min_degree=1,
  encode_sets=FALSE,
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=..count..), hjust=0.8, stat='count')
    + annotate(geom='text', label = '@', x='Number of MAGs', y=120, color='white', size=2)
    + theme(axis.text.x=element_text(angle=90))
  )+ patchwork::plot_layout(heights=c(0.5, 1.5, 1))
)

MAGs_scaled_over_90_denite <- upset(
  all_MAGs_over_90,
  denite_genes,
  sort_intersections_by=c('degree'),
  sort_sets=FALSE,
  base_annotations = list(
    'Intersection size'=intersection_size(
      mapping=aes(fill=Phylum))
    + scale_fill_manual(values=P34)
    + ylab('Number of MAGs in Gene Set')
  ),
  min_size=1,
  width_ratio=0.1,
  min_degree=1,
  encode_sets=FALSE,
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=..count..), hjust=0.8, stat='count')
    + annotate(geom='text', label = '@', x='Number of MAGs', y=120, color='white', size=2)
    + theme(axis.text.x=element_text(angle=90))
  )
)



MAGs_scaled_except_surface <- upset(
  all_MAGs_except_surface,
  genes,
  sort_intersections_by=c('degree'),
  sort_sets=FALSE,
  base_annotations = list(
    'Intersection size'=intersection_size(
      mapping=aes(fill=Phylum))
    + scale_fill_manual(values=P34)
    + ylab('Number of MAGs in Gene Set')
  ),
  annotations = list(
    'Completeness'=list(
      aes=aes(x=intersection, y=completeness),
      geom=list(
        geom_jitter(na.rm=TRUE),
        geom_violin(alpha=0.5, na.rm=TRUE)
      )
    )
  ),
  min_size=5,
  width_ratio=0.1,
  min_degree=1,
  encode_sets=FALSE,
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=..count..), hjust=0.8, stat='count')
    + annotate(geom='text', label = '@', x='Number of MAGs', y=120, color='white', size=2)
    + theme(axis.text.x=element_text(angle=90))
  )
)+ patchwork::plot_layout(heights=c(0.5, 1.5, 1))

MAGs_scaled_dereplicated <- upset(
  all_dereplicated_MAGs,
  genes,
  sort_intersections_by=c('degree'),
  sort_sets=FALSE,
  base_annotations = list(
    'Intersection size'=intersection_size(
      mapping=aes(fill=Phylum))
    + scale_fill_manual(values=P34)
    + ylab('Number of MAGs in Gene Set')
  ),
  annotations = list(
    'Completeness'=list(
      aes=aes(x=intersection, y=completeness),
      geom=list(
        geom_jitter(na.rm=TRUE),
        geom_violin(alpha=0.5, na.rm=TRUE, size=0.2)
      )
    )
  ),
  min_size=3,
  width_ratio=0.1,
  min_degree=1,
  encode_sets=FALSE,
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=..count..), hjust=0.8, stat='count')
    + annotate(geom='text', label = '@', x='Number of MAGs', y=120, color='white', size=2)
    + theme(axis.text.x=element_text(angle=90))
  )
)+ patchwork::plot_layout(heights=c(0.5, 1.5, 1))

MAGs_scaled_except_surface_dereplicated <- upset(
  all_MAGs_except_surface_derep,
  genes,
  sort_intersections_by=c('degree'),
  sort_sets=FALSE,
  base_annotations = list(
    'Intersection size'=intersection_size(
      mapping=aes(fill=Phylum))
    + scale_fill_manual(values=P34)
    + ylab('Number of MAGs in Gene Set')
  ),
  annotations = list(
    'Completeness'=list(
      aes=aes(x=intersection, y=completeness),
      geom=list(
        geom_jitter(na.rm=TRUE, size = 0.2),
        geom_violin(alpha=0.5, na.rm=TRUE)
      )
    )
  ),
  min_size=4,
  width_ratio=0.2,
  min_degree=1,
  encode_sets=FALSE,
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=..count..), hjust=0.8, stat='count')
    + annotate(geom='text', label = '@', x='Number of MAGs', y=120, color='white', size=2)
    + theme(axis.text.x=element_text(angle=90))
  )
)


print(MAGs_scaled_all)
print(MAGs_scaled_denite_all)
print(MAGs_scaled_except_surface)
print(MAGs_scaled_dereplicated)
print(MAGs_scaled_except_surface_dereplicated)
print(dereplicated_MAGs_over_90_plot_scaled)
print(MAGs_scaled_over_90)
print(MAGs_scaled_over_90_denite)


setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="Sets_with_taxonomy_all_newnors.svg", 
       MAGs_scaled_all)

#ggsave(filename="Sets_with_taxonomy_no_surface_new.svg", 
       #MAGs_scaled_except_surface)

ggsave(filename="Sets_with_taxonomy_derep_newnors.svg", 
       MAGs_scaled_dereplicated)
                 

#ggsave(filename="Sets_taxonomy_no_surface_derep.svg", 
       #MAGs_scaled_except_surface_dereplicated)

#ggsave(filename="Sets_taxonomy_over_90.svg", 
       #MAGs_scaled_over_90)

ggsave(filename="Sets_taxonomy_denite_only_newnors.svg", 
       MAGs_scaled_denite_all)

ggsave(filename="Sets_taxonomy_denite_over90_newnors.svg", 
       MAGs_scaled_over_90_denite)


#check out nor intersections
nors = list("nor_qc", "nor_weird")
as.character(nors)

MAGs_nors <- upset(
  all_MAGs,
  nors,
  sort_intersections_by=c('degree'),
  sort_sets=FALSE,
  base_annotations = list(
    'Intersection size'=intersection_size(
      mapping=aes(fill=Phylum))
    + scale_fill_manual(values=P34)
    + ylab('Number of MAGs in Nor Set')
  ),
  annotations = list(
    'Completeness'=list(
      aes=aes(x=intersection, y=completeness),
      geom=list(
        geom_jitter(na.rm=TRUE, size = 0.2),
        geom_violin(alpha=0.5, na.rm=TRUE)
      )
    )
  ),
  min_size=4,
  width_ratio=0.2,
  min_degree=1,
  encode_sets=FALSE,
  set_sizes=(
    upset_set_size()
    + geom_text(aes(label=..count..), hjust=0.8, stat='count')
    + annotate(geom='text', label = '@', x='Number of MAGs', y=120, color='white', size=2)
    + theme(axis.text.x=element_text(angle=90))
  )
)

MAGs_nors
setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")
ggsave(filename="Sets_with_taxonomy_nors.svg", 
       MAGs_nors)
