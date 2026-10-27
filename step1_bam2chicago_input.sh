#!/bin/bash
#SBATCH --mem=240G
#SBATCH --cpus-per-task=2
#SBATCH --time=36:00:00
#SBATCH --mail-type=FAIL,END,TIME_LIMIT_80
module load R/3.6
module load bedtools/2.25.0

for file in /data/Brown_lab/Laufey_lab_capture_C/data/merged_bam_files/*bam
do
var=$(echo "$file"| cut -d/ -f7)
var1=$(echo "$file")
chmod 777 /scratch/thakurr2/region_capture_c/oneFrag/bam2chicago/"$var"
/data/Brown_lab/chicago/region_capture_c/Chicago_package_scriptFiles/bam2chicago.sh "$var1" /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/design_file_4F/GWAS_Mbo_4F.baitmap /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/design_file_4F/MboIhg19_4fragment.rmap /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/"$var"
done
