##Cleaned up script for identifying variants in promoters, promoter-interacting fragments, ATAC peaks, or ChromHMM active elements
##Used in Connelly et al. 2026


library("GenomicRanges")

##For promoters use Gencode transcripts
## setwd("path/to/finemapped variants_and_annotation_file")

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

## SNPs located in GENCODE PROMOTERS
SNPs<-read.table("path/to/fine-mapped/variants/coordinates", header = TRUE, sep="\t")
snp_gr<-makeGRangesFromDataFrame(SNPs,
                                 keep.extra.columns=TRUE,
                                 ignore.strand=TRUE,
                                 seqinfo=NULL,
                                 seqnames.field=c("seqnames", "seqname",
                                                  "chromosome", "chrom",
                                                  "chr","chromosome_name",
                                                  "seqid"),
                                 start.field="start",
                                 end.field=c("end", "stop"),
                                 starts.in.df.are.0based=FALSE
)                

overlaps_snps<-findOverlaps(snp_gr, proms_reduced, type = "within", select="all", ignore.strand=TRUE)
snp_promoter_overlaps <- snp_gr[queryHits(overlaps_snps)]
prom_overlaps <- proms_reduced[subjectHits(overlaps_snps)]

##create dataframe
SNP_promoter_df<-data.frame(
  SNP_chr = seqnames(snp_promoter_overlaps),
  SNP_start = start(snp_promoter_overlaps),
  SNP_end = end(snp_promoter_overlaps),
  promoter_chr = seqnames(prom_overlaps),
  promoter_start = start(prom_overlaps),
  promoter_end = end(prom_overlaps),
  tx_id = prom_overlaps$tx_id,
  gene_id = prom_overlaps$gene_id,
  biotype = prom_overlaps$biotype,
  rsID = snp_promoter_overlaps$SNP
)
write.table(SNP_promoter_df, file="SNPs in Gencode promoters.txt", sep="\t") ##Tab of Supplementary Table 9, used in Figure 

##Find SNP overlaps with PANCREAS-active promoters
Pancreas_promoters<-read.table("path/to/pancreas/active/promoter/coordinates.txt", header = TRUE, sep="\t") ## Coordinates data file generated in Promoter_Interaction_Script.R
PP_gr<-makeGRangesFromDataFrame(Pancreas_promoters,
                                keep.extra.columns=TRUE,
                                ignore.strand=TRUE,
                                seqinfo=NULL,
                                seqnames.field=c("Promoter_chr"),
                                start.field="Promoter_start",
                                end.field="Promoter_end",
                                starts.in.df.are.0based=FALSE
)

PP_snps<-findOverlaps(snp_gr, PP_gr, type = "within", select="all", ignore.strand=TRUE)
snp_PP_overlaps <- snp_gr[queryHits(PP_snps)]
pp_overlaps <- PP_gr[subjectHits(PP_snps)]

SNP_activepromoter_df<-data.frame(
  SNP_chr = seqnames(snp_PP_overlaps),
  SNP_start = start(snp_PP_overlaps),
  SNP_end=end(snp_PP_overlaps),
  rsID= snp_PP_overlaps$SNP,
  promoter_chr = seqnames(pp_overlaps),
  promoter_start = start(pp_overlaps),
  promoter_end = end(pp_overlaps),
  geneID= pp_overlaps$Promoter_gene,
  tx_ID=pp_overlaps$Promoter_transcript,
  biotype=pp_overlaps$Biotype
)
write.table(SNP_activepromoter_df, file="SNPs in pancreas-active promoters.txt", sep="\t") ##Data not presented in paper directly, but used for Supplementary Figure 5



##Identifies SNPs that lie in an interacting fragment
# First generate the All-Cell-Line (ACL) GRanges object as done in Promoter_Interaction_Script.R

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


allinteractions_snps<-findOverlaps(snp_gr, ACL_CaptureC_GR, type = "within", select="all", ignore.strand=TRUE)
snp_interaction_overlaps <- snp_gr[queryHits(allinteractions_snps)]
interaction_overlaps <- CaptureC_GR[subjectHits(allinteractions_snps)]

SNP_interaction_df<-data.frame(
  SNP_chr = seqnames(snp_interaction_overlaps),
  SNP_start = start(snp_interaction_overlaps),
  SNP_end=end(snp_interaction_overlaps),
  rsID= snp_interaction_overlaps$SNP,
  interaction_chr = seqnames(interaction_overlaps),
  interaction_start = start(interaction_overlaps),
  interaction_end = end(interaction_overlaps),
  bait_chr=interaction_overlaps$chr,
  bait_start= interaction_overlaps$bait.start,
  bait_end=interaction_overlaps$bait.end,
  InteractionValue=interaction_overlaps$NA.
)

write.table(SNP_interaction_df, file="SNPs in all interactions.txt", sep="\t") ###1278 SNPs in interactions ##not used in manuscript

## Identify variants in promoter-interacting fragments
Prom_inter<-Overlap_df #from Promoter_Interaction_Script.R, using this narrows down to the interactions only involving transcript promoters for analysis
Prom_inter_gr<-makeGRangesFromDataFrame(Prom_inter,
                                        keep.extra.columns=TRUE,
                                        ignore.strand=TRUE,
                                        seqinfo=NULL,
                                        seqnames.field=c("Bait_chr"),
                                        start.field="Bait_start",
                                        end.field=c("Bait_end"),
                                        starts.in.df.are.0based=FALSE
)         

bait_int_snps<-findOverlaps(snp_gr, Prom_inter_gr, type = "within", select="all", ignore.strand=TRUE)
snp_intbait_overlaps <- snp_gr[queryHits(bait_int_snps)]
bait_int_overlaps <- Prom_inter_gr[subjectHits(bait_int_snps)]

##create dataframe
SNP_int_bait_df<-data.frame(
  SNP_chr = as.character(seqnames(snp_intbait_overlaps)),
  SNP_start = as.integer(start(snp_intbait_overlaps)),
  SNP_end = as.integer(end(snp_intbait_overlaps)),
  Bait_chr= seqnames(bait_int_overlaps),
  Bait_start=start(bait_int_overlaps),
  Bait_end = end(bait_int_overlaps),
  rsID= snp_intbait_overlaps$SNP,
  promoter_chr=bait_int_overlaps$proms_reduced_chr,
  promoter_start=bait_int_overlaps$proms_reduced_start,
  promoter_end=bait_int_overlaps$proms_reduced_end,
  tx_id=bait_int_overlaps$tx_id,
  gene_id=bait_int_overlaps$gene_id,
  biotype=bait_int_overlaps$biotype,
  Capture_chr= bait_int_overlaps$ACL_CaptureC_GR_chr,
  Capture_start=bait_int_overlaps$ACL_CaptureC_GR_start,
  Capture_end=bait_int_overlaps$ACL_CaptureC_GR_end,
  InteractionValue=bait_int_overlaps$InteractionValue
)

write.table(SNP_int_bait_df, file="SNPs in promoter interacting fragments.txt", sep="\t") # Generates Supplementary Table 9

##SNPs in ChromHMM active promoter regions (Active Element, Active TSS, Bivalent TSS)
##This is strictly just defining elements based on Zhong et al. Gut 2025, without considering promoters. This is separate from the Active Promoters, which were overlapped with the +/-2.5 kb coordinates.
Promoter_HMM<- read.table("path/to/active/promoter/ChromHMM/coordinates", header=TRUE, sep = "\t")  
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

overlaps_snps_HMM<-findOverlaps(snp_gr, Promoter_HMM_GR, type = "within", select="all", ignore.strand=TRUE)
snp_HMM_overlaps <- snp_gr[queryHits(overlaps_snps_HMM)]
HMM_overlaps<-Promoter_HMM_GR[subjectHits(overlaps_snps_HMM)]

SNP_HMM_df<-data.frame(
  SNP_chr = seqnames(snp_HMM_overlaps),
  SNP_start = start(snp_HMM_overlaps),
  SNP_end = end(snp_HMM_overlaps),
  HMM_chr = seqnames(HMM_overlaps),
  HMM_start = start(HMM_overlaps),
  HMM_end = end(HMM_overlaps)
)
SNP_HMM_df$ChromHMM <- mcols(HMM_overlaps)$ChromHMM
SNP_HMM_df$SNP <- mcols(snp_HMM_overlaps)$SNP
SNP_HMM_df$cell_line<-mcols(HMM_overlaps)$cell_line
write.table(SNP_HMM_df, file="SNPs in ChromHMM.txt", sep="\t") ##This does not include active enhancer elements, combined with the analysis below to identify all variants in ChromHMM active elements used in Supplementary Figure 6 and Supplementary Table 12

## SNPs in ENHANCERS based on ChromHMM (Active Enhancer, Weak Enhancer 1 and Weak Enhancer 2)
Pancreas_enhancers<-read.table("path/to/pancreas/active/enhancers", header = FALSE, sep="\t") ##generated from the ChromHMM bed files from Zhong et al. Gut 2025
Enhancer_gr<-makeGRangesFromDataFrame(Pancreas_enhancers,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=TRUE,
                                      seqinfo=NULL,
                                      seqnames.field= 'V1',
                                      start.field= 'V2',
                                      end.field= 'V3',
                                      starts.in.df.are.0based=FALSE
)

enhancer_snps<-findOverlaps(snp_gr, Enhancer_gr, type = "within", select="all", ignore.strand=TRUE)
snp_enhancer_overlaps <- snp_gr[queryHits(enhancer_snps)]
enhancer_overlaps <- Enhancer_gr[subjectHits(enhancer_snps)]

SNP_enhancer_df<-data.frame(
  SNP_chr = seqnames(snp_enhancer_overlaps),
  SNP_start = start(snp_enhancer_overlaps),
  SNP_end=end(snp_enhancer_overlaps),
  rsID= snp_enhancer_overlaps$SNP,
  enhancer_chr = seqnames(enhancer_overlaps),
  enhancer_start = start(enhancer_overlaps),
  enhancer_end = end(enhancer_overlaps),
  cell_line = enhancer_overlaps$V10
)
write.table(SNP_enhancer_df, file="SNPs_within_Enhancers.txt", sep="\t") ##This was combined with pancreas-active promoter variants to identify all variants within active ChromHMM elements. Duplicate variants removed in Excel. 
##Variants used in Supplementary Figure 6 and Supplementary Table 12


##Identify interactions containing fine-mapped variants that interact with a pancreas-active promoter
Prom_inter_HMM<-Overlap_HMM_df #from extracting promoter coordinates script
Prom_inter_HMM_gr<-makeGRangesFromDataFrame(Prom_inter_HMM,
                                            keep.extra.columns=TRUE,
                                            ignore.strand=TRUE,
                                            seqinfo=NULL,
                                            seqnames.field=c("seqnames", "seqname",
                                                             "chromosome", "chrom",
                                                             "Bait_chr","chromosome_name",
                                                             "seqid"),
                                            start.field="Bait_Start",
                                            end.field=c("Bait_End"),
                                            starts.in.df.are.0based=FALSE
)         


bait_HMMint_snps<-findOverlaps(snp_gr, Prom_inter_HMM_gr, type = "within", select="all", ignore.strand=TRUE)
snp_intHMMbait_overlaps <- snp_gr[queryHits(bait_HMMint_snps)]
bait_intHMM_overlaps <- Prom_inter_HMM_gr[subjectHits(bait_HMMint_snps)]

##create dataframe
SNP_intHMM_bait_df<-data.frame(
  SNP_chr = as.character(seqnames(snp_intHMMbait_overlaps)),
  SNP_start = as.integer(start(snp_intHMMbait_overlaps)),
  SNP_end = as.integer(end(snp_intHMMbait_overlaps)),
  Bait_chr= seqnames(bait_intHMM_overlaps),
  Bait_start=start(bait_intHMM_overlaps),
  Bait_end = end(bait_intHMM_overlaps),
  rsID= snp_intHMMbait_overlaps$SNP,
  promoter_chr=bait_intHMM_overlaps$promoter_chr,
  promoter_start=bait_intHMM_overlaps$promoter_start,
  promoter_end=bait_intHMM_overlaps$promoter_end,
  tx_id=bait_intHMM_overlaps$tx_id,
  gene_id=bait_intHMM_overlaps$gene_id,
  biotype=bait_intHMM_overlaps$biotype,
  Capture_chr= bait_intHMM_overlaps$CaptureC_GR_chr,
  Capture_start=bait_intHMM_overlaps$CaptureC_GR_start,
  Capture_end=bait_intHMM_overlaps$CaptureC_GR_end,
  InteractionValue=bait_intHMM_overlaps$InteractionValue
)

write.table(SNP_intHMM_bait_df, file="SNPs in interacting HMM bait frag gencode 080526.txt", sep="\t") #Data used to generate Supplementary Table 12 and Supplementary Figures 5,6


## SNPS in ATAC peaks, both normal-derived immortalizated cell lines and PDAC-derived cell lines. Data from Zhong et al. Gut 2025
ATAC_normal <-read.table("path/to/normal/ATAC/bed/file", header =TRUE, sep = "\t")
Normal_gr<-makeGRangesFromDataFrame(ATAC_normal,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=TRUE,
                                    seqinfo=NULL,
                                    seqnames.field= 'X.cmd',
                                    start.field= 'X.',
                                    end.field= 'mergePeaks',
                                    starts.in.df.are.0based=FALSE
)

ATAC_normal_snps<-findOverlaps(snp_gr, Normal_gr, type = "within", select="all", ignore.strand=TRUE)
snp_N_ATAC_overlaps <- snp_gr[queryHits(ATAC_normal_snps)]
N_ATAC_overlaps <- Normal_gr[subjectHits(ATAC_normal_snps)]

SNP_Normal_ATAC_df<-data.frame(
  SNP_chr = seqnames(snp_N_ATAC_overlaps),
  SNP_start = start(snp_N_ATAC_overlaps),
  SNP_end=end(snp_N_ATAC_overlaps),
  rsID= snp_N_ATAC_overlaps$SNP,
  ATAC_chr = seqnames(N_ATAC_overlaps),
  ATAC_start = start(N_ATAC_overlaps),
  ATAC_end = end(N_ATAC_overlaps),
  cell_line = N_ATAC_overlaps$chr
)
write.table(SNP_Normal_ATAC_df, file="SNPs_within_Normal_ATAC.txt", sep="\t") ##Data used to generate Figure 3, Supplementary Figures 5,6 Supplementary Table 12

ATAC_PDAC <-read.table("T:/Amundadottir/ATAC/FinalPeaks_ATAC_2020/All_ATAC_PDAC.txt", header =TRUE, sep = "\t")
PDAC_gr<-makeGRangesFromDataFrame(ATAC_PDAC,
                                  seqinfo=NULL,
                                  seqnames.field= 'chr',
                                  start.field= 'start',
                                  end.field= 'end',
                                  strand.field="strand",
                                  keep.extra.columns=TRUE,
                                  starts.in.df.are.0based=FALSE
)

ATAC_PDAC_snps<-findOverlaps(snp_gr, PDAC_gr, type = "within", select="all", ignore.strand=TRUE)
snp_PDAC_ATAC_overlaps <- snp_gr[queryHits(ATAC_PDAC_snps)]
PDAC_ATAC_overlaps <- PDAC_gr[subjectHits(ATAC_PDAC_snps)]

SNP_PDAC_ATAC_df<-data.frame(
  SNP_chr = seqnames(snp_PDAC_ATAC_overlaps),
  SNP_start = start(snp_PDAC_ATAC_overlaps),
  SNP_end=end(snp_PDAC_ATAC_overlaps),
  rsID= snp_PDAC_ATAC_overlaps$SNP,
  ATAC_chr = seqnames(PDAC_ATAC_overlaps),
  ATAC_start = start(PDAC_ATAC_overlaps),
  ATAC_end = end(PDAC_ATAC_overlaps)
)
write.table(SNP_PDAC_ATAC_df, file="SNPs_within_PDAC_ATAC.txt", sep="\t") ##Data used to generate Figure 3, Supplementary Figures 5,6 Supplementary Table 12. Variants from Normal and PDAC were combined for a general ATAC-peak residing variant list.


## To quantitate the number of variants within the original baitmap designs, 1F and 4F for original comparison to verify which captured more fine-mapped variants
bait_map<-read.table("path/to/4F/bait/map", header =TRUE, sep = "\t")
baitmapGR<- makeGRangesFromDataFrame(bait_map,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo=NULL,
                                     seqnames.field=c("seqnames", "seqname",
                                                      "chromosome", "chrom",
                                                      "chr","chromosome_name",
                                                      "seqid"),
                                     start.field="start",
                                     end.field=c("end", "stop"),
                                     starts.in.df.are.0based=FALSE
)        

bait_map_1F<-read.table("path/to/1F/bait/map", header =TRUE, sep = "\t")
baitmapGR_1F<- makeGRangesFromDataFrame(bait_map_1F,
                                        keep.extra.columns=TRUE,
                                        ignore.strand=TRUE,
                                        seqinfo=NULL,
                                        seqnames.field=c("seqnames", "seqname",
                                                         "chromosome", "chrom",
                                                         "chr","chromosome_name",
                                                         "seqid"),
                                        start.field="start",
                                        end.field=c("end", "stop"),
                                        starts.in.df.are.0based=FALSE
)        

overlaps_baits_snps_4F<-findOverlaps(snp_gr, baitmapGR, type = "within", select="all", ignore.strand=TRUE)
snp_bait_overlaps_4F <- snp_gr[queryHits(overlaps_baits_snps_4F)]
bait_overlaps_4F <- baitmapGR[subjectHits(overlaps_baits_snps_4F)]

SNP_bait_4F_df<-data.frame(
  SNP_chr = seqnames(snp_bait_overlaps_4F),
  SNP_start = start(snp_bait_overlaps_4F),
  SNP_end = end(snp_bait_overlaps_4F),
  bait_chr = seqnames(bait_overlaps_4F),
  bait_start = start(bait_overlaps_4F),
  bait_end = end(bait_overlaps_4F),
  rsID=snp_bait_overlaps_4F$SNP
)

~~~~
  
overlaps_baits_snps_1F<-findOverlaps(snp_gr, baitmapGR_1F, type = "within", select="all", ignore.strand=TRUE)
snp_bait_overlaps_1F <- snp_gr[queryHits(overlaps_baits_snps_1F)]
bait_overlaps_1F <- baitmapGR_1F[subjectHits(overlaps_baits_snps_1F)]

SNP_bait_1F_df<-data.frame(
  SNP_chr = seqnames(snp_bait_overlaps_1F),
  SNP_start = start(snp_bait_overlaps_1F),
  SNP_end = end(snp_bait_overlaps_1F),
  bait_chr = seqnames(bait_overlaps_1F),
  bait_start = start(bait_overlaps_1F),
  bait_end = end(bait_overlaps_1F)
)
SNP_bait_1F_df$SNP<-mcols(snp_bait_overlaps_1F)$SNP
write.table(SNP_bait_1F_df, file="SNPs in 1F.txt", sep="\t")  