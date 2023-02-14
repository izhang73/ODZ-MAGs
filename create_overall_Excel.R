library(readxl)
library(dplyr)
library(xlsx)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/PROKKA-MAGs/PROKKA-MAGs-all/") 

#read in tab-delimited text files containing MAG name and PROKKA IDs and combine
Arabian_table <- read.delim("Arabian_PROKKA/Arabian_match_PROKKA_to_MAGs.txt", header = FALSE, sep = "\t")
ETNP_table <- read.delim("ETNP_PROKKA/ETNP_match_PROKKA_to_MAGs.txt", header = FALSE, sep = "\t")

#set column names as MAG and PROKKAID
colnames(Arabian_table)[1] <- "MAG"
colnames(Arabian_table)[2] <- "PROKKAID"
colnames(ETNP_table)[1] <- "MAG"
colnames(ETNP_table)[2] <- "PROKKAID"

all_table <- rbind(Arabian_table, ETNP_table)

#set a list of all files in the IDs_clean directory to then loop over
data_files <- list.files("new_HMM_run/all_IDs_clean/")

#loop over IDs_clean directory and read files listed, storing as hits
for(i in 1:length(data_files)) {         
  assign(paste0("hits", i),                                   
         read.delim(paste0("new_HMM_run/all_IDs_clean/", 
                          data_files[i]), header = FALSE, col.names = data_files[i]))
}

hitslist <- c(hits1, hits2, hits3, hits4, hits5, hits6, hits7, hits8, hits9, hits10, hits11)

add_cols <- mapply(cbind, hitslist, "hits"=1, SIMPLIFY=F)
add_cols_new <- mapply(unique, add_cols)

add_cols_new <-lapply(add_cols_new, function(x) {colnames(x)[1] <- "PROKKAID"; x})

merged_table <- merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(
  all_table, 
  add_cols_new$all_ODZ_amoA_AOA_IDs_clean.txt, by="PROKKAID", all = T), 
  add_cols_new$all_ODZ_amoA_AOB_IDs_clean.txt, by="PROKKAID", all = T),
  add_cols_new$all_ODZ_napA_IDs_clean.txt, by="PROKKAID", all = T),
  add_cols_new$all_ODZ_narG_IDs_clean.txt, by="PROKKAID", all = T),
  add_cols_new$all_ODZ_nirK_IDs_clean.txt, by="PROKKAID", all = T),
  add_cols_new$all_ODZ_nirS_IDs_clean.txt, by="PROKKAID", all = T),
  add_cols_new$all_ODZ_norB_IDs_clean.txt, by="PROKKAID", all = T),
  add_cols_new$all_ODZ_nosZ_IDs_clean.txt, by="PROKKAID", all = T),
  add_cols_new$all_ODZ_nifH_IDs_clean.txt, by="PROKKAID", all = T),
  add_cols_new$all_ODZ_nrfA_IDs_clean.txt, by="PROKKAID", all = T),
  add_cols_new$all_ODZ_nxrA_IDs_clean.txt, by="PROKKAID", all = T)

colnames(merged_table) <- c("PROKKAID","MAG","amoA_AOA","amoA_AOB","napA","narG","nirK","nirS","norB","nosZ","nifH","nrfA","nxrA")

merged_table<- merged_table %>% 
  mutate_at(vars(amoA_AOA, amoA_AOB, napA, narG, nirK, nirS, norB, nosZ, nifH, nrfA, nxrA), as.numeric)%>%
  mutate_all(~replace(., is.na(.), 0))


write.csv(merged_table, file = "all_ODZ_gene_presence_absence.csv")

#dropped .fa from MAG names in Excel and moved to new directory and replaced periods in MAG names with underscores

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets/") 

#read in new, edited files and join phylogeny and old HMMhits
merged_table<- read.csv("all_ODZ_gene_presence_absence.csv", row.names = NULL) 
merged_table <- merged_table[,-1]
old_HMM_table<- read.csv("MAGs_with_denite_joined.csv", row.names=NULL) 

merged_table_v1 <- merge(
  merged_table, 
  old_HMM_table, by="MAG", all = T)

write.csv(merged_table_v1, file = "merged_phylo_old_and_new_HMMhits.csv")

#edited file outside of R to combine hits from both hold and new HMM runs
#read in new, edited file
merged_table_v2 <- read_excel("merged_phylo_old_and_new_HMMhits.xlsx",
                          sheet="overall_phylo_and_HMMhits")

#read in file with completeness and contamination information
complete_contam <- read_excel("all-MAGs-info.xlsx",
                              sheet="total")

#merge completeness and contamination info with previous tables
merged_table_v3<- merge(
  merged_table_v2, 
  complete_contam, by="MAG", all = T)

#drop GTDBtk info because it only has the dereplicated genomes
merged_table_v3 <- merged_table_v3[,-3:-9]

#read in GTDBtk info for all genomes
all_phylo <- read_excel("all_MAG_GTDBtk.xlsx",
                              sheet="all-ODZ-GTDBtk")
#merge phylogeny with previous tables
merged_table_v4<- merge(
  merged_table_v3, 
  all_phylo, by="MAG", all = T)

#write an overall Excel table including: MAG name, PROKKA IDs, GTDBtk phylogeny, overall hits for all genes
  #influding denitrification genes and nxr, nrf, nif, amo
  #also includes completeness and contamination information for all genomes
write.csv(merged_table_v4, file = "ODZ_MAG_all_info.csv")

#read in new, overall file
overall_table <- read_excel("ODZ_MAG_all_info.xlsx")


setwd("/Users/irene/Desktop/Bess_Xin_Amal/PROKKA-MAGs/trimmed_alignments_trees/new_alignments_both_HMMruns") 

data_files_2 <- list.files("final_HMM_IDs_validated/all_IDs_clean/")

#loop over IDs_clean directory and read files listed, storing as hits
for(i in 1:length(data_files_2)) {         
  assign(paste0("hits", i),                                   
         read.delim(paste0("final_HMM_IDs_validated/all_IDs_clean/", 
                           data_files_2[i]), header = FALSE, col.names = data_files_2[i]))
}

hitslist_2 <- c(hits1, hits2, hits3, hits4, hits5, hits6)

add_cols <- mapply(cbind, hitslist_2, "hits"=1, SIMPLIFY=F)
add_cols_new <- mapply(unique, add_cols)

add_cols_new <-lapply(add_cols_new, function(x) {colnames(x)[1] <- "PROKKAID"; x})

merged_table_final <- merge(merge(merge(merge(merge(merge(
  overall_table,
  add_cols_new$final_napA_MAGs_IDs.txt_IDs_clean.txt, by="PROKKAID", all = T),
  add_cols_new$final_narG_MAGs_IDs.txt_IDs_clean.txt, by = "PROKKAID", all = T),
  add_cols_new$final_nirK_MAGs_IDs.txt_IDs_clean.txt, by = "PROKKAID", all = T),
  add_cols_new$final_nirS_MAGs_IDs.txt_IDs_clean.txt, by = "PROKKAID", all = T),
  add_cols_new$final_norB_MAGs_IDs.txt_IDs_clean.txt, by = "PROKKAID", all = T),
  add_cols_new$final_nosZ_MAGs_IDs.txt_IDs_clean.txt, by = "PROKKAID", all = T
)

names(merged_table_final)[29:34] <- c("napA_val", "narG_val","nirK_val","nirS_val","norB_val","nosZ_val")

merged_table_final<- merged_table_final %>% 
  mutate_at(vars(napA_val, narG_val, nirK_val, nirS_val, norB_val, nosZ_val), as.numeric)%>%
  mutate_all(~replace(., is.na(.), 0))

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets/") 
write.xlsx(merged_table_final, file = "ODZ_MAG_all_info.xlsx", sheetName="all-ODZ-info")

#add in hzo hit info
overall_table <- read.xlsx("ODZ_MAG_all_info.xlsx", 
                                  sheetIndex=1)

hzo_hits <- read.delim("/Users/irene/Desktop/Bess_Xin_Amal/PROKKA-MAGs/PROKKA-MAGs-all/new_HMM_run/all_ODZ_hzo_unique.hits", header=FALSE) 
names(hzo_hits) <- c("PROKKAID")
hzo_hits <- cbind(hzo_hits, hzo='1')

merged_table <- merge(
  overall_table,
  hzo_hits, by="PROKKAID", all = T)

write.xlsx(merged_table, file = "ODZ_MAG_all_info_v2.xlsx", sheetName="all-ODZ-info")

#add ODZ MAG all info hzo and amoA genes to CoverM_GTDB_joined
overall_table <- read_excel("ODZ_MAG_all_info.xlsx")
tax_mat <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=5)

simple_table <- overall_table[,c(2,24:26,30)]
merged_table_final <- merge(x = tax_mat, y = simple_table, by = "MAG", all.x = TRUE)
write.xlsx(merged_table_final, file = "CoverM_GTDB_joined_v2.xlsx")


#add in nod and new nors info
setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets/") 
overall_table <- read.xlsx("ODZ_MAG_all_info.xlsx", 
                           sheetIndex=1)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/PROKKA-MAGs/trimmed_alignments_trees/new_alignments_both_HMMruns/final_HMM_IDs_validated/") 
nod_hits <- read.delim("all_IDs_clean/final_nod_MAGs_IDs.txt_IDs_clean.txt", 
                       header=FALSE) 
names(nod_hits) <- c("PROKKAID")
nod_hits <- cbind(nod_hits, nod='1')

merged_table_nod <- merge(
  overall_table,
  nod_hits, by="PROKKAID", all = T)

newnor_hits <- read.delim("all_IDs_clean/final_allnors_MAGs_IDs.txt_IDs_clean.txt", 
                       header=FALSE) 
names(newnor_hits) <- c("PROKKAID")
newnor_hits <- cbind(newnor_hits, newnor='1')

merged_table_newnor <- merge(
  merged_table_nod,
  newnor_hits, by="PROKKAID", all = T)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets/") 
write.xlsx(merged_table_newnor, file = "ODZ_MAG_all_info_v3.xlsx", sheetName="all-ODZ-info")

#add ODZ MAG all info nod and new nor genes to CoverM_GTDB_joined
overall_table_newnor <- read_excel("ODZ_MAG_all_info_v3.xlsx")
tax_mat <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=5)

simple_table_newnor <- overall_table_newnor[,c(2,11:12)]
merged_table_newnor <- merge(x = tax_mat, y = simple_table_newnor, by = "MAG", all.x = TRUE)
write.xlsx(merged_table_newnor, file = "CoverM_GTDB_joined_v3.xlsx")


#separate nors in qnor/cnor and weird nors
setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets/") 
overall_table <- read.xlsx("ODZ_MAG_all_info_v3.xlsx", 
                           sheetIndex=1)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/PROKKA-MAGs/trimmed_alignments_trees/new_alignments_both_HMMruns/final_HMM_IDs_validated/") 

qnor_cnor_hits <- read.delim("all_IDs_clean/all_ODZ_cnor_qnor_IDs_clean.txt", 
                          header=FALSE) 
names(qnor_cnor_hits ) <- c("PROKKAID")
qnor_cnor_hits <- cbind(qnor_cnor_hits, nor_qc='1')

merged_table_norqc <- merge(
  overall_table,
  qnor_cnor_hits, by="PROKKAID", all = T)

weirdnor_hits <- read.delim("all_IDs_clean/all_ODZ_weirdnors_IDs_clean.txt", 
                             header=FALSE) 
names(weirdnor_hits ) <- c("PROKKAID")
weirdnor_hits <- cbind(weirdnor_hits, nor_weird='1')

merged_table_weird <- merge(
  merged_table_norqc,
  weirdnor_hits, by="PROKKAID", all = T)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets/") 
write.xlsx(merged_table_weird, file = "ODZ_MAG_all_info_v4.xlsx", sheetName="all-ODZ-info")

#add ODZ MAG all info nod and new nor genes to CoverM_GTDB_joined
overall_table_nor <- read_excel("ODZ_MAG_all_info_v4.xlsx")
tax_mat <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=7)

simple_table_weirdnor <- overall_table_nor[,c(2,36:37)]
merged_table_weirdnor <- merge(x = tax_mat, y = simple_table_weirdnor, by = "MAG", all.x = TRUE)
write.xlsx(merged_table_weirdnor, file = "CoverM_GTDB_joined_v3.xlsx")

#merge with chemot info
tax_mat_nors <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=7)
tax_mat_chemot <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=6)

simple_table_chemot <- tax_mat_chemot[,c(1,26:38)]
merged_table_nor_chemot <- merge(x = tax_mat_nors, y = simple_table_chemot, by = "MAG", all.x = TRUE)
write.xlsx(merged_table_nor_chemot, file = "CoverM_GTDB_joined_v4.xlsx")


