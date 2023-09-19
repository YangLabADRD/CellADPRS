## Cell_ADPRS Ver 7. Clean (Aug 2023)

library(plyr)
library(dplyr)
library(nlme)
library(data.table)
library(MASS)
library(ggplot2)
library(mediation)
require(utils)
library(UpSetR)
library(pwr)
library(ivmodel)
library(fmsb)
library(lavaan)
library(semPlot)
library(sensemakr)
library(fmsb)

source("./Functions_Cell_ADPRS_upload_18SEP2023.R")
source("./demosum_HSY_CellADPRS_upload_18SEP2023.R")

### Study Participants
# d0: ROSMAP phenotype data download 04/25/2022 (www.radc.rush.edu)
d2<-subset(d0, !is.na(gbatch) & !is.na(apoe4n) & !is.na(EV01) & !is.na (Mic_Sg_30kb) & !is.na(age_death) & (!is.na(amyloid)|!is.na(tangles)))
# gbatch: genotyping batch
# EV: genotype principal components

# mathys_meta<-read.table("..[Mathys H et al Nature 2019]/filtered_column_metadata.txt",header=TRUE)
mathys_ID<-unique(mathys_meta$projid)
sum(mathys_ID %in% d0$projid)
sum(mathys_ID %in% subset(d0,pathoAD==0)$projid)
mathys_ID_noAD<-mathys_ID[mathys_ID %in% subset(d0,pathoAD==0)$projid]
# n=24

d2<-subset(d0, !is.na(gbatch) & !is.na(apoe4n) & !is.na(EV01) & !is.na (Mic_Sg_100kb) & !is.na(age_death) & (!is.na(amyloid)|!is.na(tangles)))
d2<-subset(d2,!(projid %in% mathys_ID_noAD))
# Exclude participants overlapping with the Mathys snRNA dataset
# n=1457


### #1 "Study Participants" (Table 1)
demosum_ROSMAP(data=d2)
#        Variables Mean.Number SD.Percent    n
# 1   Age at Death    89.74455  6.5458797 1457
# 2     Female (%)   973.00000 66.7810570 1457
# 3 White Race (%)  1455.00000 99.8627316 1457
# 4 Education(yrs)    16.27404  3.5514538 1456
# 5 E4 carrier (%)   376.00000 25.8064516 1457
# 6        pathoAD   954.00000  0.6547701 1457
# 7      Last MMSE    24.00000 13.0000000 1456

summary.factor(d2$cogdx>=4 & d2$pathoAD==1) # AD with dementia
# FALSE  TRUE 
# 918   539
539/1457
# [1] 0.3699382

# Amyloid +/-
summary.factor(d2$ceradsc)
#   1   2   3   4 
# 485 523 119 330 
#1 (definite) or 2 (probable)

(485+523)/(485+523+119+330)
# [1] 0.6918325

a<-demosum_A4(a_d1);a
#        Variables Mean.Number  SD.Percent    n
# 1   Baseline Age   71.371955   4.7636526 2921
# 2     Female (%) 1740.000000  59.5686409 2921
# 3 White Race (%) 2921.000000 100.0000000 2921
# 4 Education(yrs)   16.746146   2.6948339 2919
# 5 E4 carrier (%) 1038.000000  35.5357754 2921
# 6    Baseline Ab    1.098138   0.1943858 2921
# 7  Baseline MMSE   29.000000   2.0000000 2921

# Amyloid-B +/-
table(a_d1$Ab_read)
# negative positive 
#     2030      890
890/(2030+890)*100
# [1] 30.47945


### #2 "Derivation of cell-type-specific ADPRS" (Figure 1, Supplementary Table 1)
# PRS-CS: see https://github.com/getian107/PRScs
# We used AD GWAS summary statistics from: PMID: 35379992
# We used bim files (SNP list) without APOE plus minus 1 MB
# reference panel: 1000G phase 3 (see https://github.com/getian107/PRScs)

# Derivation of cell-type-specific gene list:
# Cell type annotation defined as in Mathys H et al Nature 2019 
# avg_exp: data frame, average expression of each gene in each cell type (*n=24 no AD subjects only)
#       Seurat, AverageExpression command
# per_exp: data frame, percent of each cell type expressing each gene (*n=24 no AD subjects only)
#       Seurat, Percent_Expressing command (scCustomize)

ead(avg_exp);head(per_exp)
names(avg_exp)[1]<-"Genes";names(per_exp)[1]<-"Genes"
names(avg_exp)[2:7]<-paste0(names(avg_exp)[2:7],"_expr")
names(avg_exp)[2:7]<-substring(names(avg_exp)[2:7],5,)
names(per_exp)[2:7]<-paste0(names(per_exp)[2:7],"_percent")
head(avg_exp);head(per_exp)

expr<-merge(avg_exp,per_exp,by="Genes")

expr$All_expr<-rowSums(expr[,c(2:7)])

# Calculating Sg value, following Skene Nat Genet 2018 and Bryois Nat Genet 2020
for (i in names(expr)[2:7]){
  j<-substring(i,1,nchar(i)-5)
  j1<-paste0(j,"_Sg")
  expr[,j1]<-expr[,i]/expr[,14]
}

# Gene coordinates (hg38) - select autosomal genes not within 1 MB from APOE
# hg38: gene coordinates in BED format (downloaded from PLINK website 9/15/2022)
names(hg38)<-c("Chr","Pos1","Pos2","Gene")
# n=27033
sum(duplicated(hg38$Gene))
# n=1410 # to be removed
dup<-hg38$Gene[duplicated(hg38$Gene)]
dup<-unique(dup)
# n=614 genes
hg38<-subset(hg38,!(Gene %in% dup))
# n=25009

hg38<-subset(hg38,Chr %in% c(1:22))
# n=23995
hg38[hg38$Gene=="APOE",]
#     Chr     Pos1     Pos2 Gene
# 950  19 44905781 44909393 APOE

hg38<-subset(hg38,!(Chr==19 & ((43905781<Pos1)&(Pos2<45909393))))

# Select genes with at least 1% expression in at least one cell type
expr_1<-subset(expr,(Ex_percent>1|Oli_percent>1|In_percent>1|Mic_percent>1|Ast_percent>1|Opc_percent>1))
expr_1<-subset(expr_1,Genes %in% hg38$Gene)
# --> 13,438 autosomal genes expressed in >1% of cells from one or more cell types, excluding APOE region
# --> top 10%: n=1343

# Generate cell-type specific gene lists
l_genes<-list()
for (i in 1:6){
  a<-expr_1[,c("Genes",names(expr_1)[i+14],names(expr_1)[i+7])]
  a<-a[a[,3]>1,]
  a1<-a[order(a[,2],decreasing=TRUE),]$Genes
  b<-1343
  a2<-a1[c(1:b)]
  l_genes_auto[[(i)]]<-a2
}
names(l_genes)<-names(expr[15:20])

# UpSetR - Fig 1b
p<-upset(fromList(l_genes_auto),sets=c(cells),keep.order=TRUE,order.by="freq",
         mainbar.y.label="Intersections",sets.x.label="Cell-type-specific genes",
         text.scale=c(1.2,1.2,1.2,1.2,1.2,1.05),nintersects = 15)

# 30 kb, 0 kb, and 100 kb tracks
for (i in c(1:6)){
  a<-subset(hg38,Gene %in% l_genes[[i]])
  assign(paste0("hg38_",as.character(names(l_genes[i]))),a)
  a<-a[order(a$Chr,a$Pos1),]
  # write as a file, derive PRS with PLINK --extract range flag
  
  a_30<-a
  a_30$Pos1<-a_30$Pos1-30000
  a_30$Pos2<-a_30$Pos2+30000
  a_30$Pos1[a_30$Pos1<0]<-0
  assign(paste0("hg38_",as.character(names(l_genes[i])),"_30"),a_30)
  # write as a file, derive PRS with PLINK --extract range flag
  
  a_100<-a
  a_100$Pos1<-a_100$Pos1-100000
  a_100$Pos2<-a_100$Pos2+100000
  a_100$Pos1[a_100$Pos1<0]<-0
  assign(paste0("hg38_",as.character(names(l_genes[i])),"_100"),a_100)
  # write as a file, derive PRS with PLINK --extract range flag
}

# Use this list and BED file (BED tracks) to generate PRS: use PRS-CS posterior effect size + above BED tracks for cell-type-specific PRS calculation
# Also used above BED tracks to calculate cell-type-specific PRS, using PRSet (p=1 as recommended by Choi SW et al. PLOS Genetics 2023)

## Generating non-overlapping cell-type-specific gene tracks to assess Mic-ADPRS-independent effects
Mic_noOli<-l_genes$Mic_Sg[!(l_genes$Mic_Sg %in% l_genes$Oli_Sg)]
# n=1207
Oli_noMic<-l_genes$Oli_Sg[!(l_genes$Oli_Sg %in% l_genes$Mic_Sg)]
# n=1207
Ast_noMic<-l_genes$Ast_Sg[!(l_genes$Ast_Sg %in% l_genes$Mic_Sg)]
# n=1184
Ex_noMic<-l_genes$Ex_Sg[!(l_genes$Ex_Sg %in% l_genes$Mic_Sg)]
# n=1327

hg38_Mic_noOli<-subset(hg38,Gene %in% Mic_noOli)
# n=1207
hg38_Oli_noMic<-subset(hg38,Gene %in% Oli_noMic)
hg38_Ast_noMic<-subset(hg38,Gene %in% Ast_noMic)
hg38_Ex_noMic<-subset(hg38,Gene %in% Ex_noMic)

head(hg38_Mic_noOli)
hg38_Mic_noOli_30<-hg38_Mic_noOli
hg38_Mic_noOli_30$Pos1<-hg38_Mic_noOli_30$Pos1-30000
hg38_Mic_noOli_30$Pos2<-hg38_Mic_noOli_30$Pos2+30000
hg38_Mic_noOli_30$Pos1[hg38_Mic_noOli_30$Pos1<0]<-0
# write as a file, derive PRS with PLINK --extract range flag

head(hg38_Oli_noMic)
hg38_Oli_noMic_30<-hg38_Oli_noMic
hg38_Oli_noMic_30$Pos1<-hg38_Oli_noMic_30$Pos1-30000
hg38_Oli_noMic_30$Pos2<-hg38_Oli_noMic_30$Pos2+30000
hg38_Oli_noMic_30$Pos1[hg38_Oli_noMic_30$Pos1<0]<-0
# write as a file, derive PRS with PLINK --extract range flag

head(hg38_Ast_noMic)
hg38_Ast_noMic_30<-hg38_Ast_noMic
hg38_Ast_noMic_30$Pos1<-hg38_Ast_noMic_30$Pos1-30000
hg38_Ast_noMic_30$Pos2<-hg38_Ast_noMic_30$Pos2+30000
hg38_Ast_noMic_30$Pos1[hg38_Ast_noMic_30$Pos1<0]<-0
# write as a file, derive PRS with PLINK --extract range flag

head(hg38_Ex_noMic)
hg38_Ex_noMic_30<-hg38_Ex_noMic
hg38_Ex_noMic_30$Pos1<-hg38_Ex_noMic_30$Pos1-30000
hg38_Ex_noMic_30$Pos2<-hg38_Ex_noMic_30$Pos2+30000
hg38_Ex_noMic_30$Pos1[hg38_Ex_noMic_30$Pos1<0]<-0
# write as a file, derive PRS with PLINK --extract range flag

# PRSs from ROSMAP and A4 generated (PLINK --extract range flag; PRSet with default [no p-value thresholing] setting) 
# and merged to the phenotype data. 

# d2: ROSMAP data frame (n=1457), a_d1: A4 data frame (n=2921)

# All PRSs standardized before phenotype analyses.
snRNAname_30<-c("Ex_Sg_30kb","In_Sg_30kb","Ast_Sg_30kb","Mic_Sg_30kb","Oli_Sg_30kb","Opc_Sg_30kb")
snRNAname_all30<-c("All_Bellenguez",snRNAname_30)
for (i in snRNAname_all30){
  a<-d2[,i]
  a_std<-(a-mean(a))/sd(a)
  d2[,paste0(i,"_std")]<-a_std
}

for (i in snRNAname_all){
  a<-a_d1[,i]
  a_std<-(a-mean(a))/sd(a)
  a_d1[,paste0(i,"_std")]<-a_std
}
names(a_d1)

for (i in snRNAname_all30){
  a<-a_d1[,i]
  a_std<-(a-mean(a))/sd(a)
  a_d1[,paste0(i,"_std")]<-a_std
}
names(a_d1)

snRNAname_all30_std<-paste0(snRNAname_all30,"_std")

## Heatmap with 30 kb - Fig1c/d
# ROSMAP
PRSs<-snRNAname_all30_std

markers<-as.matrix(d2[,c(PRSs)]) 
row.names(markers)<-d2$projid

a1<-PRSs;a2<-matrix(NA,1,length(a1))

for (i in 1:length(a1)){
  a2[i]<-substr(a1[i],1,nchar(a1[i]))
}
a1;a2

colnames(markers)<-gsub("_.*","",a2)
# colnames(markers)[2]<-"Ctx"
markers<-markers[complete.cases(markers),]
# markers.mat<-round(cor(markers,method="pearson"),2) # All r's are positive
markers.mat<-(round(cor(markers,method="pearson"),2))^2 # Will use (R2)

# melted_markers.mat<-melt(markers.mat)

# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
# Reorder the correlation matrix
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# cormat <- reorder_cormat(markers.mat)
cormat<-markers.mat
cormat
#         Ex     In    Ast    Mic    Oli    Opc
# Ex  1.0000 0.0169 0.0036 0.0049 0.0036 0.0081
# In  0.0169 1.0000 0.0081 0.0064 0.0036 0.0289
# Ast 0.0036 0.0081 1.0000 0.0121 0.0100 0.0961
# Mic 0.0049 0.0064 0.0121 1.0000 0.3136 0.0256
# Oli 0.0036 0.0036 0.0100 0.3136 1.0000 0.0400
# Opc 0.0081 0.0289 0.0961 0.0256 0.0400 1.0000

cormat>0.1
# All False except Oli-Mic

# upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(cormat, na.rm = TRUE)
# melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat_edit<-melted_cormat #will censor rho<0.1
# for (i in 1:nrow(melted_cormat)){
#   if(abs(melted_cormat_edit[i,3])<0.1){
#     melted_cormat_edit[i,3]<-0
#  }
# }

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", 
                       limit = c(0,1), space = "Lab", 
                       name="R-squared") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap

## A4
markers<-as.matrix(a_d1[,c(PRSs)]) 
row.names(markers)<-a_d1$projid

a1<-PRSs;a2<-matrix(NA,1,length(a1))

for (i in 1:length(a1)){
  a2[i]<-substr(a1[i],1,nchar(a1[i]))
}
a1;a2
colnames(markers)<-gsub("_.*","",a2)
markers<-markers[complete.cases(markers),]
markers.mat<-(round(cor(markers,method="pearson"),2))^2 # Will use (R2)

cormat<-markers.mat
cormat
#         Ex     In    Ast    Mic    Oli    Opc
# Ex  1.0000 0.0225 0.0064 0.0121 0.0081 0.0121
# In  0.0225 1.0000 0.0100 0.0064 0.0025 0.0289
# Ast 0.0064 0.0100 1.0000 0.0100 0.0144 0.0841
# Mic 0.0121 0.0064 0.0100 1.0000 0.2601 0.0100
# Oli 0.0081 0.0025 0.0144 0.2601 1.0000 0.0256
# Opc 0.0121 0.0289 0.0841 0.0100 0.0256 1.0000

cormat>0.1
# All FALSE exept Mic-Oli

melted_cormat <- melt(cormat, na.rm = TRUE)

ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", 
                       limit = c(0,1), space = "Lab", 
                       name="R-squared") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap

## ROSMAP (Fig. 2, supplementary tables 2-13)
# Functions to calculate associations between ADPRS and traits are in Functions_Cell_ADPRS_ver6_3pc_08242023.R
# (For both ROSMAP and A4 analyses)

# d2_ADdem: subset for Fig2a (Control vs. AD with dementia)
d2$ADdem<-as.numeric((d2$cogdx>=4 & d2$pathoAD==1))
d2_ADdem<-subset(d2,(cogdx==1 & pathoAD==0)|((cogdx==4|cogdx==5|cogdx==6) & pathoAD==1))

## Revised main result code:
## ROSMAP (Fig. 2, supplementary tables 2-13)
# Functions to calculate associations between ADPRS and traits are in Functions_Cell_ADPRS_ver6_3pc_08242023.R
# (For both ROSMAP and A4 analyses)

# d2_ADdem: subset for Fig2a (Control vs. AD with dementia)
d2$ADdem<-as.numeric((d2$cogdx>=4 & d2$pathoAD==1))
d2_ADdem<-subset(d2,(cogdx==1 & pathoAD==0)|((cogdx==4|cogdx==5|cogdx==6) & pathoAD==1))

# Code for the main results (ROSMAP)
a1<-ADdem_PRScs_ROSMAP(PRS = snRNAname_all30_std, data=d2_ADdem);b<-boot_R2.AD(data=d2,predictors=snRNAname_all30_std)
a1<-merge(a1,b,by="PRS",all.x=TRUE);a1$PRS<-paste0(a1$PRS,"_AD")

a2<-amyloid_PRScs_ROSMAP(PRS = snRNAname_all30_std,data = d2);b<-boot_R2(data=d2,trait="amyloid_sqrt",predictors=snRNAname_all30_std)
a2<-merge(a2,b,by="PRS",all.x=TRUE);a2$PRS<-paste0(a2$PRS,"_A")

a3<-dp_PRScs_ROSMAP(PRS = snRNAname_all30_std,data = d2);b<-boot_R2(data=d2,trait="plaq_d_sqrt",predictors=snRNAname_all30_std)
a3<-merge(a3,b,by="PRS",all.x=TRUE);a3$PRS<-paste0(a3$PRS,"_dp")

a4<-np_PRScs_ROSMAP(PRS = snRNAname_all30_std,data = d2);b<-boot_R2(data=d2,trait="plaq_n_sqrt",predictors=snRNAname_all30_std)
a4<-merge(a4,b,by="PRS",all.x=TRUE);a4$PRS<-paste0(a4$PRS,"_np")

a5<-tangles_PRScs_ROSMAP(PRS = snRNAname_all30_std,data = d2);b<-boot_R2(data=d2,trait="tangles_sqrt",predictors=snRNAname_all30_std)
a5<-merge(a5,b,by="PRS",all.x=TRUE);a5$PRS<-paste0(a5$PRS,"_T")

a6<-nft_PRScs_ROSMAP(PRS = snRNAname_all30_std,data = d2);b<-boot_R2(data=d2,trait="nft_sqrt",predictors=snRNAname_all30_std)
a6<-merge(a6,b,by="PRS",all.x=TRUE);a6$PRS<-paste0(a6$PRS,"_nft")

a7<-cogdec_PRScs_ROSMAP(PRS = snRNAname_all30_std,data = d2);b<-boot_R2.cog(data=d2,predictors=snRNAname_all30_std)
a7<-merge(a7,b,by="PRS",all.x=TRUE);a7$PRS<-paste0(a7$PRS,"_C")

A_ROSMAP<-rbind(a1,a2);A_ROSMAP<-rbind(A_ROSMAP,a3);A_ROSMAP<-rbind(A_ROSMAP,a4)
A_ROSMAP<-rbind(A_ROSMAP,a5);A_ROSMAP<-rbind(A_ROSMAP,a6);A_ROSMAP<-rbind(A_ROSMAP,a7)
A_ROSMAP$FDR<-p.adjust(A_ROSMAP$p,method="fdr") #study-wise FDR (ROSMAP) - threshold 0.025

m1<-glm(factor(ADdem)~All_Bellenguez_std+apoe4n+apoe2n+age_death+msex+educ+gbatch+EV01+EV02+EV03,data=d2_ADdem,family="binomial");summary(m1)
m2<-glm(factor(ADdem)~All_Bellenguez_std+apoe2n+age_death+msex+educ+gbatch+EV01+EV02+EV03,data=d2_ADdem,family="binomial");summary(m2)
NagelkerkeR2(m1)$R2 - NagelkerkeR2(m2)$R2
# [1] 0.1167123
m2<-glm(factor(ADdem)~All_Bellenguez_std+apoe4n+age_death+msex+educ+gbatch+EV01+EV02+EV03,data=d2_ADdem,family="binomial");summary(m2)
NagelkerkeR2(m1)$R2 - NagelkerkeR2(m2)$R2
# [1] 0.02308067

m<-lm(amyloid_sqrt~All_Bellenguez_std+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.1327
m<-lm(amyloid_sqrt~All_Bellenguez_std+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Multiple R-squared:  0.06173,	Adjusted R-squared:  0.05557 
0.1327-0.05557 
# [1] 0.07713
m<-lm(amyloid_sqrt~All_Bellenguez_std+apoe4n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Multiple R-squared:  0.1249,	Adjusted R-squared:  0.1191 
0.1327-0.1191
# [1] 0.0136

m<-lm(plaq_d_sqrt~All_Bellenguez_std+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.09234
m<-lm(plaq_d_sqrt~All_Bellenguez_std+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.03086
0.09234-0.03086
# [1] 0.06148
m<-lm(plaq_d_sqrt~All_Bellenguez_std+apoe4n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.08142
0.09234-0.08142
# [1] 0.01092

m<-lm(plaq_n_sqrt~All_Bellenguez_std+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.147 
m<-lm(plaq_n_sqrt~All_Bellenguez_std+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.06551
0.147-0.06551
# [1] 0.08149
m<-lm(plaq_n_sqrt~All_Bellenguez_std+apoe4n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.1304
0.147-0.1304
# [1] 0.0166

m<-lm(tangles_sqrt~All_Bellenguez_std+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.1566
m<-lm(tangles_sqrt~All_Bellenguez_std+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.09487 
0.1566-0.09487
# [1] 0.06173
m<-lm(tangles_sqrt~All_Bellenguez_std+apoe4n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.1506 
0.1566-0.1506
# [1] 0.006

m<-lm(nft_sqrt~All_Bellenguez_std+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.1543
m<-lm(nft_sqrt~All_Bellenguez_std+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.08994 
0.1543-0.08994
# [1] 0.06436
m<-lm(nft_sqrt~All_Bellenguez_std+apoe4n+age_death+msex+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.1439
0.1543-0.1439
# [1] 0.0104

m<-lm(cogng_demog_slope~All_Bellenguez_std+apoe4n+apoe2n+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.09821 
m<-lm(cogng_demog_slope~All_Bellenguez_std+apoe2n+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.02842 
0.09821-0.02842
# [1] 0.06979
m<-lm(cogng_demog_slope~All_Bellenguez_std+apoe4n+gbatch+EV01+EV02+EV03,data=d2);summary(m)
# Adjusted R-squared:  0.09432
0.09821-0.09432
# [1] 0.00389

## ROSMAP (Fig. 3) - mediation model / structural equation modeling
# Fig 3a. Ast --> DP --> NP
# Note: ACME, average causal mediated effects. ADE, average direct effects. treat=independent variable.  
set.seed(1005)
m1<-lm(plaq_d_sqrt~Ast_Sg_30kb_std+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=subset(d2,!is.na(plaq_n_sqrt)));summary(m1)
m2<-lm(plaq_n_sqrt~plaq_d_sqrt+Ast_Sg_30kb_std+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=subset(d2));summary(m2)
m_med_AstDpNp<-mediate(m1,m2,sims=10000,boot=TRUE,treat="Ast_Sg_30kb_std",mediator="plaq_d_sqrt")
summary(m_med_AstDpNp)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            0.02343      0.00688         0.04  0.0054 ** 
# ADE             0.02794      0.00880         0.05  0.0040 ** 
# Total Effect    0.05137      0.02614         0.08  0.0002 ***
# Prop. Mediated  0.45617      0.18423         0.74  0.0056 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 1452 
# Simulations: 10000 

# Fig 3b. Ast --> NP --> NFT
# Note: ACME, average causal mediated effects. ADE, average direct effects. treat=independent variable.  
set.seed(1205)
m1<-lm(plaq_n_sqrt~Ast_Sg_30kb_std+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=subset(d2,!is.na(nft_sqrt)));summary(m1)
m2<-lm(nft_sqrt~plaq_n_sqrt+Ast_Sg_30kb_std+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=subset(d2));summary(m2)
m_med_AstNpNft<-mediate(m1,m2,sims=10000,boot=TRUE,treat="Ast_Sg_30kb_std",mediator="plaq_n_sqrt")
summary(m_med_AstNpNft)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            0.02262      0.01157         0.03  <2e-16 ***
# ADE             0.01220     -0.00429         0.03  0.1486    
# Total Effect    0.03482      0.01505         0.05  0.0002 ***
# Prop. Mediated  0.64967      0.38198         1.24  0.0002 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 1474 
# Simulations: 10000 

# Fig 3c. Mic --> NP --> NFT
# Note: ACME, average causal mediated effects. ADE, average direct effects. treat=independent variable.  
set.seed(1235)
m1<-lm(plaq_n_sqrt~Mic_Sg_30kb_std+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=subset(d2,!is.na(nft_sqrt)));summary(m1)
m2<-lm(nft_sqrt~plaq_n_sqrt+Mic_Sg_30kb_std+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03,data=subset(d2));summary(m2)
m_med_MicNpNft<-mediate(m1,m2,sims=10000,boot=TRUE,treat="Mic_Sg_30kb_std",mediator="plaq_n_sqrt")
summary(m_med_MicNpNft)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                Estimate 95% CI Lower 95% CI Upper p-value    
# ACME             0.0238       0.0127         0.03  <2e-16 ***
# ADE              0.0305       0.0136         0.05  <2e-16 ***
# Total Effect     0.0543       0.0341         0.07  <2e-16 ***
# Prop. Mediated   0.4379       0.2676         0.66  <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 1474
# Simulations: 10000 

# Fig 3d. Mic --> NFT --> CogDec (adj. NP)
# Note: ACME, average causal mediated effects. ADE, average direct effects. treat=independent variable.  
set.seed(1334)
m1<-lm(nft_sqrt~Mic_Sg_30kb_std+plaq_n_sqrt+apoe4n+apoe2n+gbatch+EV01+EV02+EV03,data=subset(d2,!is.na(cogng_demog_slope)));summary(m1)
m2<-lm(cogng_demog_slope~nft_sqrt+Mic_Sg_30kb_std+plaq_n_sqrt+apoe4n+apoe2n+gbatch+EV01+EV02+EV03,data=subset(d2));summary(m2)
m_med_MicNftCog<-mediate(m1,m2,sims=10000,boot=TRUE,treat="Mic_Sg_30kb_std",mediator="nft_sqrt")
summary(m_med_MicNftCog)
# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the Percentile Method
#                Estimate 95% CI Lower 95% CI Upper p-value    
# ACME           -0.00176     -0.00304         0.00  0.0058 ** 
# ADE            -0.00549     -0.00965         0.00  0.0060 ** 
# Total Effect   -0.00725     -0.01150         0.00  0.0002 ***
# Prop. Mediated  0.24276      0.07953         0.56  0.0060 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 1392 
# Simulations: 10000 

## SEM (Fig 3d)
a<-subset(d2,!is.na(nft) & !is.na(plaq_d) & !is.na(plaq_n) & !is.na(cogng_demog_slope))
a$plaq_d_resid<-resid(lm(plaq_d_sqrt~apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03+EV04+EV05+EV06+EV07+EV08+EV09+EV010,data=a))
a$plaq_n_resid<-resid(lm(plaq_n_sqrt~apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03+EV04+EV05+EV06+EV07+EV08+EV09+EV010,data=a))
a$nft_resid<-resid(lm(nft_sqrt~apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03+EV04+EV05+EV06+EV07+EV08+EV09+EV010,data=a))
a$cogdec_resid<-resid(lm(cogng_demog_slope~apoe4n+apoe2n+gbatch+EV01+EV02+EV03+EV04+EV05+EV06+EV07+EV08+EV09+EV010,data=a))

model2<-'
plaq_d_resid ~ Ast_Sg_30kb_std
plaq_n_resid ~ Mic_Sg_30kb_std + Ast_Sg_30kb_std + plaq_d_resid 
nft_resid ~ Mic_Sg_30kb_std + plaq_n_resid 
cogdec_resid ~ Mic_Sg_30kb_std + nft_resid + plaq_n_resid
'

set.seed(1778)
model2.fitted<-sem(model2, data=subset(a), se="boot", bootstrap=10000)
summary(model2.fitted, standardized=TRUE, rsquare=TRUE, fit.measures=TRUE)
# lavaan 0.6.13 ended normally after 11 iterations
# Estimator                                         ML
# Optimization method                           NLMINB
# Number of model parameters                        13


semPaths(model2.fitted,"std",edge.label.cex=1.5, label.font=1, curvePivot = TRUE, 
         nCharNodes=0, sideMan=50, Weighted=TRUE, esize=6,asize=6,
         node.width=3, sizeMan2=6, layout="tree2",residuals=FALSE,thresholds = TRUE,
         exoCov=FALSE, edge.label.position=0.65)

### A4 (Fig. 4)
## A4 Main (with delta.R2)
a1<-Ab_PRScs(PRS = snRNAname_all30_std,data = a_d1);a1$PRS<-paste0(a1$PRS,"_A")

m<-lm(Ab_Composite_Summary~All_Bellenguez_std+apoe4n+apoe2n+PTAGE+PTGENDER+pc1+pc2+pc3,data=a_d1);summary(m)
# Multiple R-squared:  0.1877,	Adjusted R-squared:  0.1855 
m<-lm(Ab_Composite_Summary~All_Bellenguez_std+apoe2n+PTAGE+PTGENDER+pc1+pc2+pc3,data=a_d1);summary(m)
# Multiple R-squared:  0.04149,	Adjusted R-squared:  0.03918 
0.1855 - 0.03918
# [1] 0.14632
m<-lm(Ab_Composite_Summary~All_Bellenguez_std+apoe4n+PTAGE+PTGENDER+pc1+pc2+pc3,data=a_d1);summary(m)
# Multiple R-squared:  0.1852,	Adjusted R-squared:  0.1833 
0.1855 - 0.1833 
# [1] 0.0022

a2<-FTP_composite_PRScs(PRS = snRNAname_all30_std,data = a_d1);a2$PRS<-paste0(a2$PRS,"_T")

m<-lm(FTP_Composite_Unweighted~All_Bellenguez_std+apoe4n+apoe2n+PTAGE+PTGENDER+pc1+pc2+pc3,data=a_d1);summary(m)
# Multiple R-squared:  0.1162,	Adjusted R-squared:  0.09207 
m<-lm(FTP_Composite_Unweighted~All_Bellenguez_std+apoe2n+PTAGE+PTGENDER+pc1+pc2+pc3,data=a_d1);summary(m)
# Multiple R-squared:  0.08871,	Adjusted R-squared:  0.06701 
0.09207  - 0.06701
# [1] 0.02506
m<-lm(FTP_Composite_Unweighted~All_Bellenguez_std+apoe4n+PTAGE+PTGENDER+pc1+pc2+pc3,data=a_d1);summary(m)
# Multiple R-squared:  0.1023,	Adjusted R-squared:  0.08091 
0.09207  - 0.08091 
# [1] 0.01116

a3<-HV_PRScs(PRS = snRNAname_all30_std,data = a_d1);a3$PRS<-paste0(a3$PRS,"_HV")

m<-lm(HippVol~All_Bellenguez_std+apoe4n+apoe2n+PTAGE+PTGENDER+ICV+pc1+pc2+pc3,data=a_d1);summary(m)
# Multiple R-squared:  0.3317,	Adjusted R-squared:  0.3269 
m<-lm(HippVol~All_Bellenguez_std+apoe2n+PTAGE+PTGENDER+ICV+pc1+pc2+pc3,data=a_d1);summary(m)
# Multiple R-squared:  0.3258,	Adjusted R-squared:  0.3215 
0.3269  - 0.3215 
# [1] 0.0054
m<-lm(HippVol~All_Bellenguez_std+apoe4n+PTAGE+PTGENDER+ICV+pc1+pc2+pc3,data=a_d1);summary(m)
# Multiple R-squared:  0.3313,	Adjusted R-squared:  0.3271 
0.3269  - 0.3271 
# [1] -2e-04

a4<-PACC_PRScs(PRS = snRNAname_all30_std,data = a_d1);a4$PRS<-paste0(a4$PRS,"_PACC")

m<-lm(PACC_RN~All_Bellenguez_std+apoe4n+apoe2n+PTAGE+PTGENDER+PTEDUCAT+pc1+pc2+pc3,data=a_d1);summary(m)
# Adjusted R-squared:  0.191
m<-lm(PACC_RN~All_Bellenguez_std+apoe2n+PTAGE+PTGENDER+PTEDUCAT+pc1+pc2+pc3,data=a_d1);summary(m)
# Multiple R-squared:  0.1902,	Adjusted R-squared:  0.188 
0.191  - 0.188
# [1] 0.003
m<-lm(PACC_RN~All_Bellenguez_std+apoe4n+PTAGE+PTGENDER+PTEDUCAT+pc1+pc2+pc3,data=a_d1);summary(m)
# Multiple R-squared:  0.1935,	Adjusted R-squared:  0.1913 
# NA (p>0.05)

A_A4<-rbind(a1,a2)
A_A4<-rbind(A_A4,a3)
A_A4<-rbind(A_A4,a4)
A_A4$FDR<-p.adjust(A_A4$p,method="fdr")
View(A_A4)

# Last edit: 9/18/2023 HSY











