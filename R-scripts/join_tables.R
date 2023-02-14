library(dplyr)
library(readxl)
library(xlsx) 

setwd("/Users/irene/Desktop/Bess_Xin_Amal/PROKKA-MAGs") 

#read in all MAG tables

#ETNP surface samples
AMAL1_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                          sheet="AMAL1")
AMAL5_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                          sheet="AMAL5")
AMAL13_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                           sheet="AMAL13")
AMAL17_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                           sheet="AMAL17")

#Arabian Sea 
AMAL9_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                          sheet="AMAL9")
AMAL10_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                           sheet="AMAL10")
AMAL11_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                           sheet="AMAL11")
AMAL12_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                           sheet="AMAL12")

#ETNP oxycline
AMAL18_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                           sheet="AMAL18")
AMAL14_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                           sheet="AMAL14")
AMAL2_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                          sheet="AMAL2")
AMAL7_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                          sheet="AMAL7")
AMAL8_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                          sheet="AMAL8")

#ETNP ODZ core
AMAL3_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                          sheet="AMAL3")
AMAL15_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                           sheet="AMAL15")
AMAL20_table <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                           sheet="AMAL20")

#read in the gene presence absence table
Presence_absence_all <- read_excel("Bess-Xin-Amal-MAGs.xlsx",
                          sheet="Presence-absence-all")



AMAL2_join <- left_join(AMAL2_table, Presence_absence_all, by='PROKKAID')
AMAL3_join <- left_join(AMAL3_table, Presence_absence_all, by='PROKKAID')
AMAL7_join <- left_join(AMAL7_table, Presence_absence_all, by='PROKKAID')
AMAL8_join <- left_join(AMAL8_table, Presence_absence_all, by='PROKKAID')
AMAL9_join <- left_join(AMAL9_table, Presence_absence_all, by='PROKKAID')
AMAL10_join <- left_join(AMAL10_table, Presence_absence_all, by='PROKKAID')
AMAL11_join <- left_join(AMAL11_table, Presence_absence_all, by='PROKKAID')
AMAL12_join <- left_join(AMAL12_table, Presence_absence_all, by='PROKKAID')
AMAL14_join <- left_join(AMAL14_table, Presence_absence_all, by='PROKKAID')
AMAL15_join <- left_join(AMAL15_table, Presence_absence_all, by='PROKKAID')
AMAL18_join <- left_join(AMAL18_table, Presence_absence_all, by='PROKKAID')
AMAL20_join <- left_join(AMAL20_table, Presence_absence_all, by='PROKKAID')

#write full abundance table for all NTUs in all samples to csv file
write.xlsx(AMAL2_join, file = "AMAL_MAG_tables_combined.xlsx", sheetName="AMAL2", append=TRUE)
write.xlsx(AMAL3_join, file = "AMAL_MAG_tables_combined.xlsx", sheetName="AMAL3", append=TRUE)
write.xlsx(AMAL7_join, file = "AMAL_MAG_tables_combined.xlsx", sheetName="AMAL7", append=TRUE)
write.xlsx(AMAL8_join, file = "AMAL_MAG_tables_combined.xlsx", sheetName="AMAL8", append=TRUE)
write.xlsx(AMAL9_join, file = "AMAL_MAG_tables_combined.xlsx", sheetName="AMAL9", append=TRUE)
write.xlsx(AMAL10_join, file = "AMAL_MAG_tables_combined.xlsx", sheetName="AMAL10", append=TRUE)
write.xlsx(AMAL11_join, file = "AMAL_MAG_tables_combined.xlsx", sheetName="AMAL11", append=TRUE)
write.xlsx(AMAL12_join, file = "AMAL_MAG_tables_combined.xlsx", sheetName="AMAL12", append=TRUE)
write.xlsx(AMAL14_join, file = "AMAL_MAG_tables_combined.xlsx", sheetName="AMAL14", append=TRUE)
write.xlsx(AMAL15_join, file = "AMAL_MAG_tables_combined.xlsx", sheetName="AMAL15", append=TRUE)
write.xlsx(AMAL18_join, file = "AMAL_MAG_tables_combined.xlsx", sheetName="AMAL18", append=TRUE)
write.xlsx(AMAL20_join, file = "AMAL_MAG_tables_combined.xlsx", sheetName="AMAL20", append=TRUE)

