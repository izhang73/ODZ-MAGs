library(readxl)
library(dplyr)
library(xlsx)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets/") 
overall_table <- read_excel("ODZ_MAG_all_info.xlsx", 
                          sheet="all-ODZ-validated-info")

setwd("/Users/irene/Desktop/Bess_Xin_Amal/PROKKA-MAGs/PROKKA-MAGs-all/che_mot_genes_HMMrun") 

#set a list of all files in the IDs_clean directory to then loop over
data_files <- list.files("HMMout/IDs_clean_unique/")

#loop over IDs_clean directory and read files listed, storing as hits
for(i in 1:length(data_files)) {         
  assign(paste0("hits", i),                                   
         read.delim(paste0("HMMout/IDs_clean_unique/", 
                          data_files[i]), header = FALSE, col.names = data_files[i]))
}

hitslist <- c(hits1, hits2, hits3, hits4, hits5, hits6, hits7, hits8, hits9)

add_cols <- mapply(cbind, hitslist, "hits"=1, SIMPLIFY=F)
add_cols_new <- mapply(unique, add_cols)

add_cols_new <-lapply(add_cols_new, function(x) {colnames(x)[1] <- "PROKKAID"; x})

merged_table3 <- merge(
  overall_table, 
  add_cols_new$all_aer_IDs_clean.txt, by="PROKKAID", all = T) 
colnames(merged_table3)[33] = "aer"
merged_table3 <- merge(
  merged_table3, 
  add_cols_new$all_cheA_IDs_clean.txt, by="PROKKAID", all = T) 
colnames(merged_table3)[34] = "cheA"
merged_table3 <- merge(
  merged_table3, 
  add_cols_new$all_cheB_IDs_clean.txt, by="PROKKAID", all = T) 
colnames(merged_table3)[35] = "cheB"
merged_table3 <- merge(
  merged_table3, 
  add_cols_new$all_cheR_IDs_clean.txt, by="PROKKAID", all = T) 
colnames(merged_table3)[36] = "cheR"
merged_table3 <- merge(
  merged_table3, 
  add_cols_new$all_fliG_IDs_clean.txt, by="PROKKAID", all = T) 
colnames(merged_table3)[37] = "fliG"
merged_table3 <- merge(
  merged_table3, 
  add_cols_new$all_fliM_IDs_clean.txt, by="PROKKAID", all = T) 
colnames(merged_table3)[38] = "fliM"
merged_table3 <- merge(
  merged_table3, 
  add_cols_new$all_fliNY_IDs_clean.txt, by="PROKKAID", all = T) 
colnames(merged_table3)[39] = "fliNY"
merged_table3 <- merge(
  merged_table3, 
  add_cols_new$all_motAC_IDs_clean.txt, by="PROKKAID", all = T) 
colnames(merged_table3)[40] = "motAC"
merged_table3 <- merge(
  merged_table3, 
  add_cols_new$all_motBD_IDs_clean.txt, by="PROKKAID", all = T) 
colnames(merged_table3)[41] = "motBD"

write.csv(merged_table3, file = "all_ODZ_gene_presence_absence_chemot.csv")

#saved csv file under new name in directory
setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets/") 
write.csv(merged_table3, file = "ODZ_MAG_all_info_chemot.csv")

#add motility and chemotaxis info to 
tax_mat <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=5)
overall_table <- read_excel("ODZ_MAG_all_info_chemot.xlsx", 
                            sheet="ODZ_MAG_all_info_chemot")

simple_table <- overall_table[,c(2,33:44)]
merged_table_chemot <- merge(x = tax_mat, y = simple_table, by = "MAG", all.x = TRUE)
write.xlsx(merged_table_chemot, file = "CoverM_GTDB_joined_chemot.xlsx")

#add in archaeal motility 
setwd("/Users/irene/Desktop/Bess_Xin_Amal/PROKKA-MAGs/PROKKA-MAGs-all/che_mot_genes_HMMrun/archaeal/HMMhits/IDs") 
flaG_hits <- read.delim("all_flaG_IDs.txt", 
                             header=FALSE) 
names(flaG_hits) <- c("PROKKAID")
flaG_hits <- cbind(flaG_hits, flaG='1')

flaH_hits <- read.delim("all_flaH_IDs.txt", 
                        header=FALSE) 
names(flaH_hits) <- c("PROKKAID")
flaH_hits <- cbind(flaH_hits, flaH='1')

flaJ_hits <- read.delim("all_flaJ_IDs.txt", 
                        header=FALSE) 
names(flaJ_hits) <- c("PROKKAID")
flaJ_hits <- cbind(flaJ_hits, flaJ='1')

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets/") 
overall_table <- read_excel("ODZ_MAG_all_info_chemot.xlsx", 
                            sheet="ODZ_MAG_all_info_chemot")
merged_table_arch <- merge(
  overall_table, 
  flaG_hits, by="PROKKAID", all = T) 
merged_table_arch <- merge(
  merged_table_arch, 
  flaH_hits, by="PROKKAID", all = T) 
merged_table_arch <- merge(
  merged_table_arch, 
  flaJ_hits, by="PROKKAID", all = T) 

merged_table_arch <- merged_table_arch[,-2:-32]

overall_table_new <- read_excel("ODZ_MAG_all_info.xlsx", 
                            sheet="all_ODZ_validated_newnors")
merged_table_all <- merge(
  overall_table_new, 
  merged_table_arch, by="PROKKAID", all = T) 

write.xlsx(merged_table_all, file = "ODZ_MAG_all_info_chemot_v2.xlsx")

tax_mat <- read.xlsx2("CoverM_GTDB_joined.xlsx", sheetIndex=9)
overall_table <- read_excel("ODZ_MAG_all_info.xlsx", 
                            sheet="ODZ_validated_chemot")

simple_table <- overall_table[,c(2,38:51)]
merged_table_chemot <- merge(x = tax_mat, y = simple_table, by = "MAG", all.x = TRUE)
write.xlsx(merged_table_chemot, file = "CoverM_GTDB_joined_chemot.xlsx")
