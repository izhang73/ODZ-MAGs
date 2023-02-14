library(readxl)
library(dplyr)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/all_ODZ_MAGs/all_MAGs_derep_99_all")

#read in all MAG tables from CoverM
AMAL_taxa <- read_excel("all_MAGs_derep_CoverM.xlsx",
                          sheet="all-AMAL")

#subset MAG tables by sample ID
AMAL1_CoverM <- subset(AMAL_taxa, Sample=="AMAL1")
AMAL2_CoverM <- subset(AMAL_taxa, Sample=="AMAL2")
AMAL3_CoverM <- subset(AMAL_taxa, Sample=="AMAL3")
AMAL5_CoverM <- subset(AMAL_taxa, Sample=="AMAL5")
AMAL7_CoverM <- subset(AMAL_taxa, Sample=="AMAL7")
AMAL8_CoverM <- subset(AMAL_taxa, Sample=="AMAL8")
AMAL9_CoverM <- subset(AMAL_taxa, Sample=="AMAL9")
AMAL10_CoverM <- subset(AMAL_taxa, Sample=="AMAL10")
AMAL11_CoverM <- subset(AMAL_taxa, Sample=="AMAL11")
AMAL12_CoverM <- subset(AMAL_taxa, Sample=="AMAL12")
AMAL13_CoverM <- subset(AMAL_taxa, Sample=="AMAL13")
AMAL14_CoverM <- subset(AMAL_taxa, Sample=="AMAL14")
AMAL15_CoverM <- subset(AMAL_taxa, Sample=="AMAL15")
AMAL17_CoverM <- subset(AMAL_taxa, Sample=="AMAL17")
AMAL18_CoverM <- subset(AMAL_taxa, Sample=="AMAL18")
AMAL20_CoverM <- subset(AMAL_taxa, Sample=="AMAL20")

Fuchsman_taxa <- read_excel("all_MAGs_derep_CoverM.xlsx",
                        sheet="all-Fuchsman")
Fuchsman_samples<-Fuchsman_taxa %>% distinct(Sample)
Fuchsman_SRR4465024_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465024")
Fuchsman_SRR4465025_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465025")
Fuchsman_SRR4465026_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465026")
Fuchsman_SRR4465027_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465027")
Fuchsman_SRR4465028_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465028")
Fuchsman_SRR4465029_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465029")
Fuchsman_SRR4465030_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465030")
Fuchsman_SRR4465031_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465031")
Fuchsman_SRR4465032_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465032")
Fuchsman_SRR4465033_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465033")
Fuchsman_SRR4465034_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465034")
Fuchsman_SRR4465035_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465035")
Fuchsman_SRR4465036_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465036")
Fuchsman_SRR4465037_CoverM <- subset(Fuchsman_taxa, Sample=="Fuchsman_SRR4465037")

Glass_taxa <- read_excel("all_MAGs_derep_CoverM.xlsx",
                        sheet="all-Glass")
Glass_samples<-Glass_taxa %>% distinct(Sample)
Glass_SRR1509790_CoverM <- subset(Glass_taxa, Sample=="Glass_SRR1509790")
Glass_SRR1509792_CoverM <- subset(Glass_taxa, Sample=="Glass_SRR1509792")
Glass_SRR1509793_CoverM <- subset(Glass_taxa, Sample=="Glass_SRR1509793")
Glass_SRR1509794_CoverM <- subset(Glass_taxa, Sample=="Glass_SRR1509794")
Glass_SRR1509796_CoverM <- subset(Glass_taxa, Sample=="Glass_SRR1509796")
Glass_SRR1509797_CoverM <- subset(Glass_taxa, Sample=="Glass_SRR1509797")
Glass_SRR1509798_CoverM <- subset(Glass_taxa, Sample=="Glass_SRR1509798")
Glass_SRR1509799_CoverM <- subset(Glass_taxa, Sample=="Glass_SRR1509799")
Glass_SRR1509800_CoverM <- subset(Glass_taxa, Sample=="Glass_SRR1509800")

Tsementzi_taxa <- read_excel("all_MAGs_derep_CoverM.xlsx",
                        sheet="all-Tsementzi")
Tsementzi_samples<-Tsementzi_taxa %>% distinct(Sample)
Tsementzi_SRR3718412_CoverM <- subset(Tsementzi_taxa, Sample=="Tsementzi_SRR3718412")
Tsementzi_SRR3718413_CoverM <- subset(Tsementzi_taxa, Sample=="Tsementzi_SRR3718413")

Stewart_taxa <- read_excel("all_MAGs_derep_CoverM.xlsx",
                        sheet="all-Stewart")
Stewart_samples<-Stewart_taxa %>% distinct(Sample)
Stewart_SRR064444_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR064444")
Stewart_SRR064448_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR064448")
Stewart_SRR064450_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR064450")
Stewart_SRR304656_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR304656")
Stewart_SRR304668_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR304668")
Stewart_SRR304671_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR304671")
Stewart_SRR304672_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR304672")
Stewart_SRR304673_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR304673")
Stewart_SRR304674_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR304674")
Stewart_SRR304680_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR304680")
Stewart_SRR304684_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR304684")
Stewart_SRR070081_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR070081")
Stewart_SRR070082_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR070082")
Stewart_SRR304683_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR304683")
Stewart_SRR070084_CoverM <- subset(Stewart_taxa, Sample=="Stewart_SRR070084")

Ganesh_taxa <- read_excel("all_MAGs_derep_CoverM.xlsx",
                        sheet="all-Ganesh")
Ganesh_samples<-Ganesh_taxa %>% distinct(Sample)
Ganesh_SRR960580_CoverM <- subset(Ganesh_taxa, Sample=="Ganesh_SRR960580")
Ganesh_SRR961673_CoverM <- subset(Ganesh_taxa, Sample=="Ganesh_SRR961673")
Ganesh_SRR961676_CoverM <- subset(Ganesh_taxa, Sample=="Ganesh_SRR961676")
Ganesh_SRR961679_CoverM <- subset(Ganesh_taxa, Sample=="Ganesh_SRR961679")
Ganesh_SRR961671_CoverM <- subset(Ganesh_taxa, Sample=="Ganesh_SRR961671")
Ganesh_SRR961675_CoverM <- subset(Ganesh_taxa, Sample=="Ganesh_SRR961675")
Ganesh_SRR961677_CoverM <- subset(Ganesh_taxa, Sample=="Ganesh_SRR961677")
Ganesh_SRR961680_CoverM <- subset(Ganesh_taxa, Sample=="Ganesh_SRR961680")

#make MAG table of relative abundances by doing a bunch of full joins
df1 <- AMAL1_CoverM %>% full_join(AMAL5_CoverM,by="Genome") %>% full_join(AMAL13_CoverM,by="Genome")%>% full_join(AMAL17_CoverM,by="Genome")
df2 <- AMAL9_CoverM %>% full_join(AMAL10_CoverM,by="Genome") %>% full_join(AMAL11_CoverM,by="Genome")%>% full_join(AMAL12_CoverM,by="Genome")
df3 <- AMAL18_CoverM %>% full_join(AMAL14_CoverM,by="Genome") %>% full_join(AMAL2_CoverM,by="Genome")%>% full_join(AMAL7_CoverM,by="Genome")%>% full_join(AMAL8_CoverM,by="Genome")%>% full_join(AMAL3_CoverM,by="Genome")%>% full_join(AMAL15_CoverM,by="Genome")
AMAL_df <- df1 %>% full_join(df2,by="Genome") %>% full_join(df3,by="Genome")

df4 <- Fuchsman_SRR4465024_CoverM %>% full_join(Fuchsman_SRR4465025_CoverM,by="Genome") %>% full_join(Fuchsman_SRR4465026_CoverM,by="Genome")%>% full_join(Fuchsman_SRR4465027_CoverM,by="Genome")
df5 <- Fuchsman_SRR4465028_CoverM %>% full_join(Fuchsman_SRR4465029_CoverM,by="Genome") %>% full_join(Fuchsman_SRR4465030_CoverM,by="Genome")%>% full_join(Fuchsman_SRR4465031_CoverM,by="Genome") %>% full_join(Fuchsman_SRR4465037_CoverM,by="Genome")
df6 <- Fuchsman_SRR4465032_CoverM %>% full_join(Fuchsman_SRR4465033_CoverM,by="Genome") %>% full_join(Fuchsman_SRR4465034_CoverM,by="Genome")%>% full_join(Fuchsman_SRR4465035_CoverM,by="Genome") %>% full_join(Fuchsman_SRR4465036_CoverM,by="Genome")
Fuchsman_df <- df4 %>% full_join(df5,by="Genome") %>% full_join(df6,by="Genome")

df7 <- Glass_SRR1509790_CoverM %>% full_join(Glass_SRR1509792_CoverM,by="Genome") %>% full_join(Glass_SRR1509793_CoverM,by="Genome")%>% full_join(Glass_SRR1509794_CoverM,by="Genome")
df8 <- Glass_SRR1509796_CoverM %>% full_join(Glass_SRR1509797_CoverM,by="Genome") %>% full_join(Glass_SRR1509798_CoverM,by="Genome")%>% full_join(Glass_SRR1509799_CoverM,by="Genome") %>% full_join(Glass_SRR1509800_CoverM,by="Genome")
Glass_df <- df7 %>% full_join(df8,by="Genome")

Tsementzi_df <- Tsementzi_SRR3718412_CoverM %>% full_join(Tsementzi_SRR3718413_CoverM,by="Genome")

df9 <- Stewart_SRR064444_CoverM %>% full_join(Stewart_SRR064448_CoverM,by="Genome") %>% full_join(Stewart_SRR064450_CoverM,by="Genome")%>% full_join(Stewart_SRR304656_CoverM,by="Genome")
df10 <- Stewart_SRR304668_CoverM %>% full_join(Stewart_SRR304671_CoverM,by="Genome") %>% full_join(Stewart_SRR304672_CoverM,by="Genome")%>% full_join(Stewart_SRR304673_CoverM,by="Genome") %>% full_join(Stewart_SRR304674_CoverM,by="Genome")
df11 <- Stewart_SRR304680_CoverM %>% full_join(Stewart_SRR304684_CoverM,by="Genome") %>% full_join(Stewart_SRR070081_CoverM,by="Genome")%>% full_join(Stewart_SRR070082_CoverM,by="Genome") %>% full_join(Stewart_SRR304683_CoverM,by="Genome")
Stewart_df <- df9 %>% full_join(df10,by="Genome") %>% full_join(df11,by="Genome") %>% full_join( Stewart_SRR070084_CoverM,by="Genome")

df12 <- Ganesh_SRR960580_CoverM %>% full_join(Ganesh_SRR961673_CoverM,by="Genome") %>% full_join(Ganesh_SRR961676_CoverM,by="Genome")%>% full_join(Ganesh_SRR961679_CoverM,by="Genome")
df13 <- Ganesh_SRR961671_CoverM %>% full_join(Ganesh_SRR961675_CoverM,by="Genome") %>% full_join(Ganesh_SRR961677_CoverM,by="Genome")%>% full_join(Ganesh_SRR961680_CoverM,by="Genome")
Ganesh_df <- df12 %>% full_join(df13,by="Genome")

#join all and write to a csv file
all_df <- AMAL_df %>% full_join(Fuchsman_df,by="Genome") %>% full_join(Glass_df,by="Genome") %>% full_join(Tsementzi_df,by="Genome") %>% full_join(Stewart_df,by="Genome") %>% full_join(Ganesh_df,by="Genome")
write.csv(all_df, file = "CoverM_GTDB.csv")

#work on getting GTDB taxonomy
#ETNP surface samples
AMAL1_table <- read_excel("all_MAG_GTDBtk.xlsx",
                          sheet="AMAL1-GTDBtk")
AMAL5_table <- read_excel("all_MAG_GTDBtk.xlsx",
                          sheet="AMAL5-GTDBtk")
AMAL13_table <- read_excel("all_MAG_GTDBtk.xlsx",
                           sheet="AMAL13-GTDBtk")
AMAL17_table <- read_excel("all_MAG_GTDBtk.xlsx",
                           sheet="AMAL17-GTDBtk")

#Arabian Sea 
AMAL9_table <- read_excel("all_MAG_GTDBtk.xlsx",
                          sheet="AMAL9-GTDBtk")
AMAL10_table <- read_excel("all_MAG_GTDBtk.xlsx",
                           sheet="AMAL10-GTDBtk")
AMAL11_table <- read_excel("all_MAG_GTDBtk.xlsx",
                           sheet="AMAL11-GTDBtk")
AMAL12_table <- read_excel("all_MAG_GTDBtk.xlsx",
                           sheet="AMAL12-GTDBtk")

#ETNP oxycline
AMAL18_table <- read_excel("all_MAG_GTDBtk.xlsx",
                           sheet="AMAL18-GTDBtk")
AMAL14_table <- read_excel("all_MAG_GTDBtk.xlsx",
                           sheet="AMAL14-GTDBtk")
AMAL2_table <- read_excel("all_MAG_GTDBtk.xlsx",
                          sheet="AMAL2-GTDBtk")
AMAL7_table <- read_excel("all_MAG_GTDBtk.xlsx",
                          sheet="AMAL7-GTDBtk")
AMAL8_table <- read_excel("all_MAG_GTDBtk.xlsx",
                          sheet="AMAL8-GTDBtk")

#ETNP ODZ core
AMAL3_table <- read_excel("all_MAG_GTDBtk.xlsx",
                          sheet="AMAL3-GTDBtk")
AMAL15_table <- read_excel("all_MAG_GTDBtk.xlsx",
                           sheet="AMAL15-GTDBtk")
AMAL20_table <- read_excel("all_MAG_GTDBtk.xlsx",
                           sheet="AMAL20-GTDBtk")

#Fuchsman all GTDBtk
Fuchsman_table <- read_excel("all_MAG_GTDBtk.xlsx",
                           sheet="Fuchsman-GTDBtk")

#Glass all GTDBtk
Glass_table <- read_excel("all_MAG_GTDBtk.xlsx",
                           sheet="Glass-GTDBtk")

#Tsementzi all GTDBtk
Tsementzi_table <- read_excel("all_MAG_GTDBtk.xlsx",
                           sheet="Tsementzi-GTDBtk")


#combine GTDBtk tables into one vertically with rbind
GTDB_AMAL_tables <- rbind(AMAL1_table, AMAL2_table, AMAL3_table, AMAL5_table, AMAL7_table,
                          AMAL8_table, AMAL9_table, AMAL10_table, AMAL11_table, AMAL12_table, 
                          AMAL13_table, AMAL14_table, AMAL15_table, AMAL17_table, AMAL18_table, 
                          AMAL20_table)
GTDB_other_tables <- rbind(Fuchsman_table, Glass_table, Tsementzi_table)
all_GTDB_tables <- rbind(GTDB_AMAL_tables, GTDB_other_tables)

#combine CoverM tables into one vertically with rbind
combined_taxa <- rbind(AMAL_taxa, Fuchsman_taxa, Glass_taxa, Tsementzi_taxa, Stewart_taxa,
                       Ganesh_taxa)
sorted_taxa <- combined_taxa[order(combined_taxa$Relative_Abundance, decreasing=TRUE),]
sorted_taxa[sorted_taxa==0] <- NA
filtered_taxa<-sorted_taxa[complete.cases(sorted_taxa),]

#left join with CoverM table as left table and GTDB table as right table
joined= filtered_taxa %>% left_join(all_GTDB_tables,by="Genome")

#write to a csv file
write.csv(joined, file = "CoverM_GTDB_joined.csv")

All_MAGs_tax_table <- read_excel("CoverM_GTDB_joined.xlsx",
                        sheet="tax_table")
All_MAGs_denit_table <- read_excel("CoverM_GTDB_joined.xlsx",
                                 sheet="denitrification")
joined_denite = All_MAGs_tax_table %>% left_join(All_MAGs_denit_table,by="MAG")

write.csv(joined_denite, file = "MAGs_with_denite_joined.csv")
