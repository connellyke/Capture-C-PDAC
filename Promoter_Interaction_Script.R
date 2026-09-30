## ========================================================
## Cleaned Capture-C promoter interaction analyses
## for GWAS Capture-C/7p14.1 manuscript (Connelly et al. 2026)

## This script reproduces the analyses used in:
## Figures 2b-d, f
## Supplementary Tables 7 and 8

## Libraries 

library("GenomicRanges")
library("rtracklayer")

## Paths and inputs
## setwd("path/to/CaptureC_interaction_files_and_annotation_file")

## Load in the Gencode annotation file and extract the transcript data
gtf <- import("gencode.v49lift37.annotation.gtf.gz")
tx <- gtf[gtf$type == "transcript"]

#Extract coordinates for promoters, defined +/- 2.5 kb of the TSS
proms=promoters(tx, upstream = 2500, downstream = 2500)

## Create a GRanges object keeping a subset of the data
proms_reduced <- GRanges(
  seqnames = seqnames(proms),
  ranges = ranges(proms),
  strand = strand(proms),
  gene_id = tx$gene_id,
  tx_id = tx$transcript_id,
  hgnc_name = tx$gene_name,
  biotype = tx$transcript_type,
  biotype_gene=tx$gene_type
)

## write.table(proms_reduced, file="Gencode_Promoter_coordinates_2.5kb.txt")

## All-Cell-Line (ACL) interaction data reformatted from the WashU file format and interchromosomal interactions removed

ACL_CaptureC<-read.table("All_6cultures_4F_Granges_track.txt", header = TRUE, sep = "\t")

## Transform the Capture-C data from a data frame into a GRanges object, keeping the extra dataframe columns as metadata
ACL_CaptureC_GR<-makeGRangesFromDataFrame(ACL_CaptureC,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=TRUE,
                                      seqinfo=NULL,
                                      seqnames.field=c("seqnames", "seqname",
                                                       "chromosome", "chrom",
                                                       "chr.1","chromosome_name",
                                                       "seqid"),
                                      start.field="start",
                                      end.field=c("end", "stop"),
                                      strand.field="strand",
                                      starts.in.df.are.0based=FALSE
)

#Overlaps between the CaptureC ranges and the promoter coordinates
ACL_overlaps<-findOverlaps(ACL_CaptureC_GR, proms_reduced, type = "any", select="all", ignore.strand=TRUE)
ACL_CaptureC_GR_overlaps <- ACL_CaptureC_GR[queryHits(ACL_overlaps)]
Promoters_ACL_overlaps <- proms_reduced[subjectHits(ACL_overlaps)]

## Create data frame of the promoter interactions including interacting fragment coordinates, promoter coordinates involved in interaction, transcript ID,
## ENSEMBL ID, gene ID, transcript and gene biotypes, CHiCAGO score, and Interaction Value. Output is Supplementary Table 7
Overlap_df<-data.frame(
  CaptureC_GR_chr = seqnames(ACL_CaptureC_GR_overlaps),
  CaptureC_GR_start = start(ACL_CaptureC_GR_overlaps),
  CaptureC_GR_end = end(ACL_CaptureC_GR_overlaps),
  proms_reduced_chr = seqnames(Promoters_ACL_overlaps),
  proms_reduced_start = start(ACL_Promoters_overlaps),
  proms_reduced_end = end(Promoters_ACL_overlaps),
  gene_id= Promoters_ACL_overlaps$gene_id,
  tx_id=Promoters_ACL_overlaps$tx_id,
  biotype=Promoters_ACL_overlaps$biotype,
  biotype_gene=Promoters_ACL_overlaps$biotype_gene,
  ensembl_id=Promoters_ACL_overlaps$gene_id,
  gene_id=Promoters_overlaps$hgnc_name,
  Bait_chr=CaptureC_GR_overlaps$chr,
  Bait_start=CaptureC_GR_overlaps$bait.start,
  Bait_end = CaptureC_GR_overlaps$bait.end,
  Score=CaptureC_GR_overlaps$score,
  InteractionValue = CaptureC_GR_overlaps$NA.
)
write.table(Overlap_df, file="Supplementary_Table_7.txt",sep="\t")

##Extract details on interactions to functional transcript promoters, such as # of transcripts/gene, interactions/gene, etc. Data used for Supplementary Figure 2
#subset to only functional transcripts based on listed biotypes
functional_df<- Overlap_df[Overlap_df$biotype %in% c("lncRNA","miRNA", "protein_coding", "antisense", "snoRNA", "snRNA"),]
length(unique(functional_df$tx_id))

##Find how many unique transcripts there are per gene
transcript_counts <- aggregate(tx_id ~ gene_id, data = functional_df, function(x) length(unique(x)))
median(transcript_counts$tx_id) #3
mean(transcript_counts$tx_id) #7.72
max(transcript_counts$tx_id) #114
min(transcript_counts$tx_id) #1
write.table(transcript_counts, file="Functional_Transcript_counts_per_gene.txt", sep="\t") 

##Find number of unique interactions per gene
interaction_counts <- aggregate(InteractionValue ~ gene_id, data = functional_df, function(x) length(unique(x)))
median(interaction_counts$InteractionValue) #26
mean(interaction_counts$InteractionValue) #70.28
max(interaction_counts$InteractionValue) #978
min(interaction_counts$InteractionValue) #1
write.table(interaction_counts, file="Gene~Interaction_counts_per_gene_All_061726.txt", sep="\t")


##Overlapping ChromHMM active elements (Active Element, Active TSS, Bivalent TSS) with promoters to identify epigenetically active promoters in pancreas.
#ChromHMM elements were extracted from the 8-cell line ChromHMM bed files published in Zhong et al. Gut 2025

Promoter_HMM<- read.table("Path/to/ChromHMM/elements/from/cell/lines.txt", header=TRUE, sep = "\t")
Promoter_HMM_GR<-makeGRangesFromDataFrame(Promoter_HMM,
                                          keep.extra.columns=TRUE,
                                          ignore.strand=TRUE,
                                          seqinfo=NULL,
                                          seqnames.field=c("seqnames", "seqname",
                                                           "chromosome", "chrom",
                                                           "chr","chromosome_name",
                                                           "seqid"),
                                          start.field="start",
                                          end.field="end",
                                          strand.field="strand",
                                          starts.in.df.are.0based=FALSE
)

#Find overlaps between defined ChromHMM elements and +/-2.5kb of gene transcripts listed in proms_reduced 
overlaps_HMM_promoters<-findOverlaps(Promoter_HMM_GR, proms_reduced, type = "any", select="all", ignore.strand=TRUE)
#labels the query hits with the active HMM information
HMM_overlaps <- Promoter_HMM_GR[queryHits(overlaps_HMM_promoters)]
#labels the subject hits with the promoter information
promoter_overlaps <- proms_reduced[subjectHits(overlaps_HMM_promoters)]

##Generate a dataframe for further analysis
Promoter_HMM_df<-data.frame(
  Promoter_chr = seqnames(promoter_overlaps),
  Promoter_start = start(promoter_overlaps),
  Promoter_end = end(promoter_overlaps),
  Promoter_gene=promoter_overlaps$gene_id,
  Promoter_transcript=promoter_overlaps$tx_id,
  ChromHMM_chr = seqnames(HMM_overlaps),
  ChromHMM_start = start(HMM_overlaps),
  ChromHMM_end = end(HMM_overlaps),
  Cell_line= HMM_overlaps$cell_line,
  Element_type=HMM_overlaps$ChromHMM,
  biotype = promoter_overlaps$biotype,
  biotype_gene=promoter_overlaps$biotype_gene,
  gene_name = promoter_overlaps$gene_name,
  hgnc_name =promoter_overlaps$hgnc_name
)
#Remove duplicate transcripts as an element or multiple elements may overlap with one promoter. Duplicate promoters aren't needed.  Generates Supplementary Table 10
active_unique <- Promoter_HMM_df[!duplicated(Promoter_HMM_df$Promoter_transcript), ]
write.table(active_unique, file="Supplementary_Table_10.txt", sep="\t")

#Remove duplicate tx and genes from ChromHMM overlap analysis in GRanges format for 
dup_idx <- duplicated(mcols(promoter_overlaps)$tx_id)
gr_nodup <- promoter_overlaps[!dup_idx]


#Identify overlaps between the Capture-C interactions and the epigenetically active promoter coordinates
overlaps_ChromHMM<-findOverlaps(CaptureC_GR, gr_nodup, type = "any", select="all", ignore.strand=TRUE)
CaptureC_HMM_overlaps <- CaptureC_GR[queryHits(overlaps_ChromHMM)]
HMM_overlaps <- gr_nodup[subjectHits(overlaps_ChromHMM)]

Overlap_HMM_df<-data.frame(
  CaptureC_GR_chr = seqnames(CaptureC_HMM_overlaps),
  CaptureC_GR_start = start(CaptureC_HMM_overlaps),
  CaptureC_GR_end = end(CaptureC_HMM_overlaps),
  promoter_chr = seqnames(HMM_overlaps),
  promoter_start = start(HMM_overlaps),
  promoter_end = end(HMM_overlaps),
  gene_id= HMM_overlaps$gene_id,
  tx_id=HMM_overlaps$tx_id,
  biotype=HMM_overlaps$biotype,
  biotype_gene=HMM_overlaps$biotype_gene,
  ensembl_id=HMM_overlaps$gene_id,
  hgnc_id=HMM_overlaps$hgnc_name,
  Bait_chr=CaptureC_HMM_overlaps$chr,
  Bait_start=CaptureC_HMM_overlaps$bait.start,
  Bait_end = CaptureC_HMM_overlaps$bait.end,
  Score=CaptureC_HMM_overlaps$score,
  InteractionValue = CaptureC_HMM_overlaps$NA.
)

write.table(Overlap_HMM_df, file="Supplementary_Table_11.txt",sep="\t")

##For individual cell lines
#PANC-1
CaptureC_PANC<-read.table("Path/to/PANC-1/4F/Track/file", header = TRUE, sep = "\t")
#MIA PaCa-2
CaptureC_MP<-read.table("Path/to/MIA PaCa-2/4F/Track/file", header = TRUE, sep = "\t")
#SU.86.86
CaptureC_SU<-read.table("Path/to/SU8686/4F/Track/file", header = TRUE, sep = "\t")
#Hs766T
CaptureC_HS<-read.table("Path/to/Hs766t/4F/Track/file", header = TRUE, sep = "\t")
#HPDE
CaptureC_HPDE<-read.table("Path/to/HPDE/4F/Track/file", header = TRUE, sep = "\t")
#AsPC1
CaptureC_ASPC1<-read.table("Path/to/AsPC1/4F/Track/file", header = TRUE, sep = "\t")

#Convert CaptureC data into a GRanges format using the interaction ends and keeping bait info in the metadata along with scores
CaptureC_PANC_GR<-makeGRangesFromDataFrame(CaptureC_PANC,
                                           keep.extra.columns=TRUE,
                                           ignore.strand=TRUE,
                                           seqinfo=NULL,
                                           seqnames.field=c("seqnames", "seqname",
                                                            "chromosome", "chrom",
                                                            "chr.1","chromosome_name",
                                                            "seqid"),
                                           start.field="start",
                                           end.field=c("end", "stop"),
                                           strand.field="strand",
                                           starts.in.df.are.0based=FALSE
)

#Overlaps between the CaptureC ranges and the promoter coordinates
overlaps_PANC<-findOverlaps(CaptureC_PANC_GR, proms_reduced, type = "any", select="all", ignore.strand=TRUE)
CaptureC_PANC_GR_overlaps <- CaptureC_PANC_GR[queryHits(overlaps_PANC)]
Promoters_overlaps_PANC <- proms_reduced[subjectHits(overlaps_PANC)]
Overlap_PANC_df<-data.frame(
  CaptureC_PANC_GR_chr = seqnames(CaptureC_PANC_GR_overlaps),
  CaptureC_PANC_GR_start = start(CaptureC_PANC_GR_overlaps),
  CaptureC_PANC_GR_end = end(CaptureC_PANC_GR_overlaps),
  proms_reduced_chr = seqnames(Promoters_overlaps_PANC),
  proms_reduced_start = start(Promoters_overlaps_PANC),
  proms_reduced_end = end(Promoters_overlaps_PANC),
  gene_id= Promoters_overlaps_PANC$gene_id,
  tx_id=Promoters_overlaps_PANC$tx_id,
  biotype=Promoters_overlaps_PANC$biotype,
  biotype_gene=Promoters_overlaps_PANC$biotype_gene,
  ensembl_id=Promoters_overlaps_PANC$gene_id,
  gene_id=Promoters_overlaps_PANC$hgnc_name,
  Bait_start=CaptureC_PANC_GR_overlaps$bait.start,
  Bait_end = CaptureC_PANC_GR_overlaps$bait.end,
  Score=CaptureC_PANC_GR_overlaps$score,
  InteractionValue = CaptureC_PANC_GR_overlaps$NA.
)

write.table(Overlap_PANC_df, file="PANC1.txt",sep="\t") ## Generates PANC-1 tab of Supplementary Table 8

#### MIAPACA2

CaptureC_MP_GR<-makeGRangesFromDataFrame(CaptureC_MP,
                                         keep.extra.columns=TRUE,
                                         ignore.strand=TRUE,
                                         seqinfo=NULL,
                                         seqnames.field=c("seqnames", "seqname",
                                                          "chromosome", "chrom",
                                                          "chr.1","chromosome_name",
                                                          "seqid"),
                                         start.field="start",
                                         end.field=c("end", "stop"),
                                         strand.field="strand",
                                         starts.in.df.are.0based=FALSE
)

#Overlaps between the CaptureC ranges and the promoter coordinates
overlaps_MP<-findOverlaps(CaptureC_MP_GR, proms_reduced, type = "any", select="all", ignore.strand=TRUE)
CaptureC_MP_GR_overlaps <- CaptureC_MP_GR[queryHits(overlaps_MP)]
Promoters_overlaps_MP <- proms_reduced[subjectHits(overlaps_MP)]
Overlap_MP_df<-data.frame(
  CaptureC_MP_GR_chr = seqnames(CaptureC_MP_GR_overlaps),
  CaptureC_MP_GR_start = start(CaptureC_MP_GR_overlaps),
  CaptureC_MP_GR_end = end(CaptureC_MP_GR_overlaps),
  proms_reduced_chr = seqnames(Promoters_overlaps_MP),
  proms_reduced_start = start(Promoters_overlaps_MP),
  proms_reduced_end = end(Promoters_overlaps_MP),
  gene_id= Promoters_overlaps_MP$gene_id,
  tx_id=Promoters_overlaps_MP$tx_id,
  biotype=Promoters_overlaps_MP$biotype,
  biotype_gene=Promoters_overlaps_MP$biotype_gene,
  ensembl_id=Promoters_overlaps_MP$gene_id,
  gene_id=Promoters_overlaps_MP$hgnc_name,
  Bait_start=CaptureC_MP_GR_overlaps$bait.start,
  Bait_end = CaptureC_MP_GR_overlaps$bait.end,
  Score=CaptureC_MP_GR_overlaps$score,
  InteractionValue = CaptureC_MP_GR_overlaps$NA.
)

write.table(Overlap_MP_df, file="MP2.txt",sep="\t")## Generates MIA PaCa-2 tab of Supplementary Table 8

##SU8686
CaptureC_SU_GR<-makeGRangesFromDataFrame(CaptureC_SU,
                                         keep.extra.columns=TRUE,
                                         ignore.strand=TRUE,
                                         seqinfo=NULL,
                                         seqnames.field=c("seqnames", "seqname",
                                                          "chromosome", "chrom",
                                                          "chr.1","chromosome_name",
                                                          "seqid"),
                                         start.field="start",
                                         end.field=c("end", "stop"),
                                         strand.field="strand",
                                         starts.in.df.are.0based=FALSE
)

#Overlaps between the CaptureC ranges and the promoter coordinates
overlaps_SU<-findOverlaps(CaptureC_SU_GR, proms_reduced, type = "any", select="all", ignore.strand=TRUE)
CaptureC_SU_GR_overlaps <- CaptureC_SU_GR[queryHits(overlaps_SU)]
Promoters_overlaps_SU <- proms_reduced[subjectHits(overlaps_SU)]
Overlap_SU_df<-data.frame(
  CaptureC_SU_GR_chr = seqnames(CaptureC_SU_GR_overlaps),
  CaptureC_SU_GR_start = start(CaptureC_SU_GR_overlaps),
  CaptureC_SU_GR_end = end(CaptureC_SU_GR_overlaps),
  proms_reduced_chr = seqnames(Promoters_overlaps_SU),
  proms_reduced_start = start(Promoters_overlaps_SU),
  proms_reduced_end = end(Promoters_overlaps_SU),
  gene_id= Promoters_overlaps_SU$gene_id,
  tx_id=Promoters_overlaps_SU$tx_id,
  biotype=Promoters_overlaps_SU$biotype,
  biotype_gene=Promoters_overlaps_SU$biotype_gene,
  ensembl_id=Promoters_overlaps_SU$gene_id,
  gene_id=Promoters_overlaps_SU$hgnc_name,
  Bait_start=CaptureC_SU_GR_overlaps$bait.start,
  Bait_end = CaptureC_SU_GR_overlaps$bait.end,
  Score=CaptureC_SU_GR_overlaps$score,
  InteractionValue = CaptureC_SU_GR_overlaps$NA.
)

write.table(Overlap_SU_df, file="SU8686.txt",sep="\t") ## Generates SU.86.86 tab of Supplementary Table 8

##Hs766T

CaptureC_HS_GR<-makeGRangesFromDataFrame(CaptureC_HS,
                                         keep.extra.columns=TRUE,
                                         na.rm=TRUE,
                                         ignore.strand=TRUE,
                                         seqinfo=NULL,
                                         seqnames.field=c("seqnames", "seqname",
                                                          "chromosome", "chrom",
                                                          "chr.1","chromosome_name",
                                                          "seqid"),
                                         start.field="start",
                                         end.field=c("end", "stop"),
                                         strand.field="strand",
                                         starts.in.df.are.0based=FALSE
)

#Overlaps between the CaptureC ranges and the promoter coordinates
overlaps_HS<-findOverlaps(CaptureC_HS_GR, proms_reduced, type = "any", select="all", ignore.strand=TRUE)
CaptureC_HS_GR_overlaps <- CaptureC_HS_GR[queryHits(overlaps_HS)]
Promoters_overlaps_HS <- proms_reduced[subjectHits(overlaps_HS)]
Overlap_HS_df<-data.frame(
  CaptureC_HS_GR_chr = seqnames(CaptureC_HS_GR_overlaps),
  CaptureC_HS_GR_start = start(CaptureC_HS_GR_overlaps),
  CaptureC_HS_GR_end = end(CaptureC_HS_GR_overlaps),
  proms_reduced_chr = seqnames(Promoters_overlaps_HS),
  proms_reduced_start = start(Promoters_overlaps_HS),
  proms_reduced_end = end(Promoters_overlaps_HS),
  gene_id= Promoters_overlaps_HS$gene_id,
  tx_id=Promoters_overlaps_HS$tx_id,
  biotype=Promoters_overlaps_HS$biotype,
  biotype_gene=Promoters_overlaps_HS$biotype_gene,
  ensembl_id=Promoters_overlaps_HS$gene_id,
  gene_id=Promoters_overlaps_HS$hgnc_name,
  Bait_start=CaptureC_HS_GR_overlaps$bait.start,
  Bait_end = CaptureC_HS_GR_overlaps$bait.end,
  Score=CaptureC_HS_GR_overlaps$score,
  InteractionValue = CaptureC_HS_GR_overlaps$NA.
)

write.table(Overlap_HS_df, file="HS766T.txt",sep="\t") ## Generates Hs766T tab of Supplementary Table 8

##HPDE

CaptureC_HPDE_GR<-makeGRangesFromDataFrame(CaptureC_HPDE,
                                           keep.extra.columns=TRUE,
                                           ignore.strand=TRUE,
                                           seqinfo=NULL,
                                           seqnames.field=c("seqnames", "seqname",
                                                            "chromosome", "chrom",
                                                            "chr.1","chromosome_name",
                                                            "seqid"),
                                           start.field="start",
                                           end.field=c("end", "stop"),
                                           strand.field="strand",
                                           starts.in.df.are.0based=FALSE
)

#Overlaps between the CaptureC ranges and the promoter coordinates
overlaps_HPDE<-findOverlaps(CaptureC_HPDE_GR, proms_reduced, type = "any", select="all", ignore.strand=TRUE)
CaptureC_HPDE_GR_overlaps <- CaptureC_HPDE_GR[queryHits(overlaps_HPDE)]
Promoters_overlaps_HPDE <- proms_reduced[subjectHits(overlaps_HPDE)]
Overlap_HPDE_df<-data.frame(
  CaptureC_HPDE_GR_chr = seqnames(CaptureC_HPDE_GR_overlaps),
  CaptureC_HPDE_GR_start = start(CaptureC_HPDE_GR_overlaps),
  CaptureC_HPDE_GR_end = end(CaptureC_HPDE_GR_overlaps),
  proms_reduced_chr = seqnames(Promoters_overlaps_HPDE),
  proms_reduced_start = start(Promoters_overlaps_HPDE),
  proms_reduced_end = end(Promoters_overlaps_HPDE),
  gene_id= Promoters_overlaps_HPDE$gene_id,
  tx_id=Promoters_overlaps_HPDE$tx_id,
  biotype=Promoters_overlaps_HPDE$biotype,
  biotype_gene=Promoters_overlaps_HPDE$biotype_gene,
  ensembl_id=Promoters_overlaps_HPDE$gene_id,
  gene_id=Promoters_overlaps_HPDE$hgnc_name,
  Bait_start=CaptureC_HPDE_GR_overlaps$bait.start,
  Bait_end = CaptureC_HPDE_GR_overlaps$bait.end,
  Score=CaptureC_HPDE_GR_overlaps$score,
  InteractionValue = CaptureC_HPDE_GR_overlaps$NA.
)

write.table(Overlap_HPDE_df, file="HPDE.txt",sep="\t")## Generates HPDE tab of Supplementary Table 8

##AsPC1

CaptureC_ASPC1_GR<-makeGRangesFromDataFrame(CaptureC_ASPC1,
                                            keep.extra.columns=TRUE,
                                            ignore.strand=TRUE,
                                            seqinfo=NULL,
                                            seqnames.field=c("seqnames", "seqname",
                                                             "chromosome", "chrom",
                                                             "chr.1","chromosome_name",
                                                             "seqid"),
                                            start.field="start",
                                            end.field=c("end", "stop"),
                                            strand.field="strand",
                                            starts.in.df.are.0based=FALSE
)

#Overlaps between the CaptureC ranges and the promoter coordinates
overlaps_ASPC1<-findOverlaps(CaptureC_ASPC1_GR, proms_reduced, type = "any", select="all", ignore.strand=TRUE)
CaptureC_ASPC1_GR_overlaps <- CaptureC_ASPC1_GR[queryHits(overlaps_ASPC1)]
Promoters_overlaps_ASPC1 <- proms_reduced[subjectHits(overlaps_ASPC1)]
Overlap_ASPC1_df<-data.frame(
  CaptureC_ASPC1_GR_chr = seqnames(CaptureC_ASPC1_GR_overlaps),
  CaptureC_ASPC1_GR_start = start(CaptureC_ASPC1_GR_overlaps),
  CaptureC_ASPC1_GR_end = end(CaptureC_ASPC1_GR_overlaps),
  proms_reduced_chr = seqnames(Promoters_overlaps_ASPC1),
  proms_reduced_start = start(Promoters_overlaps_ASPC1),
  proms_reduced_end = end(Promoters_overlaps_ASPC1),
  gene_id= Promoters_overlaps_ASPC1$gene_id,
  tx_id=Promoters_overlaps_ASPC1$tx_id,
  biotype=Promoters_overlaps_ASPC1$biotype,
  biotype_gene=Promoters_overlaps_ASPC1$biotype_gene,
  ensembl_id=Promoters_overlaps_ASPC1$gene_id,
  gene_id=Promoters_overlaps_ASPC1$hgnc_name,
  Bait_start=CaptureC_ASPC1_GR_overlaps$bait.start,
  Bait_end = CaptureC_ASPC1_GR_overlaps$bait.end,
  Score=CaptureC_ASPC1_GR_overlaps$score,
  InteractionValue = CaptureC_ASPC1_GR_overlaps$NA.
)

write.table(Overlap_ASPC1_df, file="AsPC1.txt",sep="\t")## Generates AsPC1 tab of Supplementary Table 8


