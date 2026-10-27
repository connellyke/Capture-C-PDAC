##Cleaned up script for RNA-seq differential gene expression analysis using RSEM count files 

##Load libraries

##BiocManager::install("DESeq2")
library(DESeq2)

##BiocManager::install("tximport")
library(tximport)
##BiocManager::install("tximportData")
library(tximportData)

##setwd("path/to/RSEM/files")

###This is where my final analysis code starts
##HPDE samples
samples_hpde<-read.table("HPDE_samples.txt", header= TRUE, sep="\t")
samples_hpde$Condition<-factor(rep(c("Control", "Knockdown"), each =3))
row.names(samples_hpde)<-samples_hpde$sample

dir<-("path/to/RSEM/files")
files<-file.path(dir, paste0(samples_hpde$SampleID, "_rsem.genes.results"))
names(files)<-paste("sample", 1:6)
txi_hpde<-tximport(files, type = "rsem", txOut = FALSE)
txi_hpde$length[txi_hpde$length == 0] <- 1
ddsTxi_hpde <- DESeqDataSetFromTximport(txi_hpde,
                                   colData = samples_hpde,
                                   design = ~ Condition)
##Sets group size
smallestGroupSize <- 3
##Filter out genes with a count less than 10 in more than 3 samples
keep <- rowSums(counts(ddsTxi_hpde) >= 10) >= smallestGroupSize
ddsTxi_hpde <- ddsTxi_hpde[keep,]
#Set condition as the factor with control as the reference sample
ddsTxi_hpde$Condition <- factor(ddsTxi_hpde$Condition, levels = c("Control","Knockdown"))
ddsTxi_hpde <- DESeq(ddsTxi_hpde)
res_hpde<- results(ddsTxi_hpde)
res_hpde
summary(res_hpde)
plotMA(res_hpde, ylim=c(-2,2))
plotCounts(ddsTxi_hpde, gene = which.min(res_hpde$padj), intgroup="Condition") # plot counts for gene with the lowest pvalue adj.
write.csv(as.data.frame(res_hpde), 
          file="hpde_inhba_knockdown_DEG.csv") ##Supplementary Table 14 data

##subset HPDE results by adj pvalue
sig_HPDE <- subset(res_hpde, padj < 0.05) ##4646 DEG 
sig_HPDE$genes <- row.names(sig_HPDE)
sig_HPDE_FC <- subset(sig_HPDE, log2FoldChange >= 1 | log2FoldChange <= -1)


##SU.86.86 samples
samples_SU<-read.table("SU86_sample.txt", header= TRUE, sep="\t")
samples_SU$Condition<-factor(rep(c("Control", "Knockdown"), each =3))
row.names(samples_SU)<-samples_SU$sample

dir<-("C:/Users/connellyke/OneDrive - National Institutes of Health/CaptureC Manuscript/RSEM_count/RSEM")
files_su<-file.path(dir, paste0(samples_SU$SampleID, "_rsem.genes.results"))
names(files_su)<-paste("sample", 1:6)
txi_su<-tximport(files_su, type = "rsem", txOut = FALSE)
txi_su$length[txi_su$length == 0] <- 1
ddsTxi_su <- DESeqDataSetFromTximport(txi_su,
                                        colData = samples_SU,
                                        design = ~ Condition)
smallestGroupSize <- 3
keep <- rowSums(counts(ddsTxi_su) >= 10) >= smallestGroupSize
ddsTxi_su <- ddsTxi_su[keep,]
ddsTxi_su$Condition <- factor(ddsTxi_su$Condition, levels = c("Control","Knockdown"))
ddsTxi_su <- DESeq(ddsTxi_su)
res_su<- results(ddsTxi_su)
res_su
summary(res_su)
plotMA(res_su, ylim=c(-2,2))
plotCounts(ddsTxi_su, gene = which.min(res_su$padj), intgroup="Condition") # plot counts for gene with the lowest pvalue adj.
write.csv(as.data.frame(res_su), 
          file="su_inhba_knockdown_DEG.csv")##Supplementary Table 14 data

sig_SU <- subset(res_su, padj < 0.05) ##1767 DEG
sig_SU$genes <- row.names(sig_SU) 
sig_SU_FC <- subset(sig_SU, log2FoldChange >= 1 | log2FoldChange <= -1)


##Examining overlap between the two cell lines
library(VennDiagram)
library(ggplot2)
library(ggVennDiagram)
set1a<-sig_HPDE$genes
set2a<-sig_SU$genes
overlap_list <- list(Set1= set1a, Set2= set2a)
ggVennDiagram(overlap_list, category.names=c("HPDE", "SU.86.86"), label_color = "white", label_alpha = 0, label= "count") ##1041 overlapping genes

set1 <-sig_HPDE_FC$genes
set2 <- sig_SU_FC$genes
overlapFC<-list(Set1=set1, Set2=set2)
ggVennDiagram(overlapFC, category.names=c("HPDE", "SU.86.86"), label_color = "white", label_alpha = 0, label= "count") #97 overlapping genes



#BiocManager::install("fgsea")
#BiocManager::install("msigdbr")

## GSEA on the differentially expressed genes, ranked by the test statistic 
##Create .rnk file for GSEA (file uploaded is generated in Excel from the DEG results file)
###RNK file for GSEA
HPDE_GSEA<-read.table("HPDE_Stat_FC_GSEA_060226.txt", header = TRUE, sep="\t")
#set up rank for HPDE GSEA
rank_HPDE<-HPDE_GSEA$stat
names(rank_HPDE)<-HPDE_GSEA$GeneID
rank_HPDE<-sort(rank_HPDE, decreasing=TRUE)

##Create .rnk file for GSEA (file uploaded is generated in Excel from the DEG results file)
SU_GSEA<- read.table("SU_Stat_FC_GSEA_060226.txt", header = TRUE, sep="\t")
#set up rank for SU.86.86 GSEA
rank_SU<-SU_GSEA$stat
names(rank_SU)<-SU_GSEA$GeneID
rank_SU<-sort(rank_SU, decreasing=TRUE)

BiocManager::install("clusterProfiler", version = "3.22")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

gse <- gseGO(geneList=rank_HPDE, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "fdr",
             seed=TRUE)
gse_res<-gse@result
hpde_simplified<-simplify(gse)
hpde_simplified_res<-as.data.frame(hpde_simplified@result)
dotplot(hpde_simplified, x= "GeneRatio", showCategory= 8, split=".sign", color="p.adj") + facet_grid(.~.sign)+ scale_size(range=c(5, 20))
write.table(hpde_simplified_res, file = "HPDE_GSEA_GO_simplified.txt") ##Supplementary Table 15


gse_SU <- gseGO(geneList=rank_SU, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "fdr",
             seed=TRUE)
dotplot(su_simplified, x= "GeneRatio", showCategory= 8, split=".sign", color="p.adj") + facet_grid(.~.sign)+ scale_size(range=c(5, 20))
SU_res<-as.data.frame(gse_SU@result)
su_simplified<-simplify(gse_SU)
su_simplified_res<-as.data.frame(su_simplified@result)
write.table(su_simplified_res, file = "SU_GSEA_GO_simplified.txt")  ##Supplementary Table 15

