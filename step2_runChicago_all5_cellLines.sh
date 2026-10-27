#!/bin/bash
#SBATCH --mem=240G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mail-type=FAIL,END,TIME_LIMIT_80

set -e
module load bedtools/2.25.0
module load perl/5.24
module load R/4.0.5

cd /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis
Rscript /data/Brown_lab/chicago/region_capture_c/Chicago_package_scriptFiles/runChicago.R --settings-file /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/script/settingsFile --export-format seqMonk,interBed,washU_text,washU_track --design-dir /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/design_file_4F /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/AsPc1.bam/AsPc1.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/Hs766T.bam/Hs766T.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/MIA.bam/MIA.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/Panc.bam/Panc.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/SU8686.bam/SU8686.bam.chinput All_5cultures_4F_bam_merged

