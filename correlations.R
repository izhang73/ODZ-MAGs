library(tidyverse)
library(readxl)
library(xlsx)

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets") 

AMAL_profile_data <- read_excel("Data_for_R_profiles.xlsx", 
                                sheet="AMAL")

num_cols <- unlist(lapply(AMAL_profile_data, is.numeric))
num_cols
profile_numeric <- AMAL_profile_data[ , num_cols]
profile_numeric

# correlation for all variables
Pearson_correlations <- round(cor(profile_numeric,use="pairwise.complete.obs"),
      digits = 2 # rounded to 2 decimals
)

Spearman_correlations <- round(cor(profile_numeric,use="pairwise.complete.obs", method ="spearman"),
                              digits = 2 # rounded to 2 decimals
)

#write.csv(Pearson_correlations, file = "AMAL_Pearson_correlations_factors.csv")
#write.csv(Spearman_correlations, file = "AMAL_Spearman_correlations_factors.csv")

library(ggplot2)
library(reshape2)

melted_pearson <- melt(Pearson_correlations)
melted_spearman <- melt(Spearman_correlations)

melted_pearson_plot <- ggplot(data = melted_pearson, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation")+
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title="Correlations") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1))+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) 

melted_spearman_plot <- ggplot(data = melted_spearman, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation")+
  labs(x = NULL, y = NULL, fill = "Spearman\nCorrelation", title="Correlations") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1))+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0))


### Get lower triangle of the correlation matrix (Pearson only)
Pearson <- Pearson_correlations
Pearson[lower.tri(Pearson,diag=TRUE)] <- NA

melted_pearson_lower <- melt(Pearson, na.rm = TRUE)
melted_pearson_plot_lower <- ggplot(data = melted_pearson_lower, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation")+
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title="Correlations") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1))+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) 


#subset only a few columns and rows from correaltions matrices
Pearson_keep <- Pearson_correlations[c(1, 3:5, 33:43),c(1, 3:5, 33:43)]
#Spearman_keep <- Spearman_correlations[c(1, 3:5, 32:40),c(1, 3:5, 32:40)]

### Get lower triangle of the correlation matrix again for subsets
Pearson <- Pearson_keep
Pearson[lower.tri(Pearson,diag=TRUE)] <- NA
Pearson


melted_pearson_lower <- melt(Pearson, na.rm = TRUE)
melted_pearson_plot_lower <- ggplot(data = melted_pearson_lower, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation")+
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title="Correlations") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1))+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) 

melted_pearson_simple <- melt(Pearson_keep)


library(psych)

#do correlations tests and use Benjamini-Hochberg multiple hypothesis correction method
Pearson_p<- corr.test(Pearson_keep, adjust="BH")
melted_pearson_p <- melt(Pearson_p$p)
#Spearman_p<- corr.test(Spearman_keep, method="spearman",adjust="BH")

melted_pearson_simple$Pvalue <- melted_pearson_p$value
colnames(melted_pearson_simple) <- c("Var1","Var2","Correlations","Pvalue")
sig_p = ifelse(melted_pearson_simple$Pvalue < .05, T, F)
p_if_sig = ifelse(melted_pearson_simple$Pvalue <.05, melted_pearson_simple$Pvalue, NA)
Pearson_if_sig = ifelse(melted_pearson_simple$Pvalue <.05, melted_pearson_simple$Correlations, NA)
melted_pearson_simple$sig_p <- sig_p
melted_pearson_simple$p_if_sig <- p_if_sig
melted_pearson_simple$Pearson_if_sig <- Pearson_if_sig


#heatmap format, keeping only significant p values
melted_pearson_plot_p_1 <- melted_pearson_simple %>% 
  ggplot(aes(Var1, Var2, fill=Correlations, label=round(p_if_sig,3))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title="Pearson Correlations AMAL", subtitle="Only significant Pearson's correlation coefficients shown") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0))

setwd("/Users/irene/Desktop/Bess_Xin_Amal/Excel_sheets") 
write.xlsx(Pearson_p$r, file = "Pearson_pvalues_BH.xlsx", sheetName="Pearson_pvalues_r")
write.xlsx(Pearson_p$p, file = "Pearson_pvalues_BH.xlsx", sheetName="Pearson_pvalues_p",append=TRUE)
write.xlsx(Pearson_p$p.adj, file = "Pearson_pvalues_BH.xlsx", sheetName="Pearson_pvalues_padj",append=TRUE)
write.xlsx(Pearson_p$se, file = "Pearson_pvalues_BH.xlsx", sheetName="Pearson_pvalues_se",append=TRUE)
write.xlsx(Pearson_p$ci, file = "Pearson_pvalues_BH.xlsx", sheetName="Pearson_pvalues_ci",append=TRUE)
write.xlsx(Pearson_p$ci2, file = "Pearson_pvalues_BH.xlsx", sheetName="Pearson_pvalues_ci2",append=TRUE)
write.xlsx(Pearson_p$ci.adj, file = "Pearson_pvalues_BH.xlsx", sheetName="Pearson_pvalues_ciadj",append=TRUE)
write.xlsx(Pearson_p$stars, file = "Pearson_pvalues_BH.xlsx", sheetName="Pearson_pvalues_stars",append=TRUE)

#ggsave(filename="Pearson_heatmap.svg", 
 #      melted_pearson_plot_p_1,bg='transparent')

#need to add a dendrogram, using heatmaply package for ths function

library(heatmaply)

ordered_heatmap_Spearman <- heatmaply(
  Spearman_keep, xlab = NULL, ylab = NULL, k_col = 3, k_row = 3,
  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "#0C6291", 
    high = "#A63446", 
    midpoint = 0, 
    limits = c(-1, 1))
  )

ordered_heatmap_Spearman

ordered_heatmap_Pearson <- heatmaply(
  Pearson_keep, xlab = NULL, ylab = NULL, k_col = 3, k_row = 3,
  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "#0C6291", 
    high = "#A63446", 
    midpoint = 0, 
    limits = c(-1, 1)),
)
ordered_heatmap_Pearson
#commands to save files as various formats - do not run unless needed

setwd("/Users/irene/Desktop/Bess_Xin_Amal/figures/new_figures_v2")

heatmaply(
  Pearson_keep, xlab = NULL, ylab = NULL, k_col = 3, k_row = 3,
  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "#0C6291", 
    high = "#A63446", 
    midpoint = 0, 
    limits = c(-1, 1)),
  file = "ordered_heatmap_Pearson.pdf"
)

heatmaply(
  Spearman_keep, xlab = NULL, ylab = NULL, k_col = 3, k_row = 3,
  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "#0C6291", 
    high = "#A63446", 
    midpoint = 0, 
    limits = c(-1, 1)),
  file = "ordered_heatmap_Spearman.pdf"
)


#plots for tests of linearity - do not run unless needed
napA_nirK <- ggplot() +
  geom_point(data=profile_numeric, aes(x=napA, y=nirK), pch=18, size= 3,colour = "blue" )

narG_nirK <- ggplot() +
  geom_point(data=profile_numeric, aes(x=narG, y=nirK), pch=18, size= 3,colour = "blue" )

nirS_nirK <- ggplot() +
  geom_point(data=profile_numeric, aes(x=nirS, y=nirK), pch=18, size= 3,colour = "blue" ) 

nosZ_nirK <- ggplot() +
  geom_point(data=profile_numeric, aes(x=nosZ, y=nirK), pch=18, size= 3,colour = "blue" ) 


norB_nirS <- ggplot() +
  geom_point(data=profile_numeric, aes(x=norB, y=nirS), pch=18, size= 3,colour = "blue" )

nosZ_nirS <- ggplot() +
  geom_point(data=profile_numeric, aes(x=nosZ, y=nirS), pch=18, size= 3,colour = "blue" )

nitrite_nirS <- ggplot() +
  geom_point(data=profile_numeric, aes(x=Nitrite, y=nirS), pch=18, size= 3,colour = "red" )

nirK_nitrite <- ggplot() +
  geom_point(data=profile_numeric, aes(x=nirK, y=Nitrite), pch=18, size= 3,colour = "red" )


nosZ_norB <- ggplot() +
  geom_point(data=profile_numeric, aes(x=nosZ, y=norB), pch=18, size= 3,colour = "blue" )

nosZ_nitrate <- ggplot() +
  geom_point(data=profile_numeric, aes(x=nosZ, y=Nitrate), pch=18, size= 3,colour = "red" )

nosZ_nitrite <- ggplot() +
  geom_point(data=profile_numeric, aes(x=nosZ, y=Nitrite), pch=18, size= 3,colour = "red" )



ggsave(filename="AMAL Pearson correlations nirK.pdf", width=16, height=5, 
       useDingbats=FALSE,
       plot_grid(napA_nirK, 
                 narG_nirK, 
                 nirS_nirK, 
                 nosZ_nirK,
                 nirK_nitrite,
                 ncol = 3))

ggsave(filename="AMAL Pearson correlations nirS.pdf", width=10, height=3, 
       useDingbats=FALSE,
       plot_grid(norB_nirS, 
                 nirS_nirK,
                 nosZ_nirS, 
                 nitrite_nirS, 
                 ncol = 3))

ggsave(filename="AMAL Pearson correlations nosZ.pdf", width=12, height=3, 
       useDingbats=FALSE,
       plot_grid(nosZ_norB, 
                 nosZ_nirS,
                 nosZ_nirK, 
                 nosZ_nitrate,
                 nosZ_nitrite, 
                 ncol = 3))
