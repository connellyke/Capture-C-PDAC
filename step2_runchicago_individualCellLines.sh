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

Rscript /data/Brown_lab/chicago/region_capture_c/Chicago_package_scriptFiles/runChicago.R --settings-file /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/script/settingsFile --export-format seqMonk,interBed,washU_text,washU_track --design-dir /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/design_file_4F /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/SU8686_Rep1_R1_2_merge.hicup.sam.bam/SU8686_Rep1_R1_2_merge.hicup.sam.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/SU8686_Rep2_R1_2_merge.hicup.sam.bam/SU8686_Rep2_R1_2_merge.hicup.sam.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/SU8686_Rep3_S1_L006_R1_2_001.hicup.sam.bam/SU8686_Rep3_S1_L006_R1_2_001.hicup.sam.bam.chinput SU8686_4F

Rscript /data/Brown_lab/chicago/region_capture_c/Chicago_package_scriptFiles/runChicago.R --settings-file /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/script/settingsFile --export-format seqMonk,interBed,washU_text,washU_track --design-dir /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/design_file_4F /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/Panc_T1_R1_2_merge.hicup.sam.bam/Panc_T1_R1_2_merge.hicup.sam.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/Panc_T2_R1_2_merge.hicup.sam.bam/Panc_T2_R1_2_merge.hicup.sam.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/Panc_T3_R1_2_merge.hicup.sam.bam/Panc_T3_R1_2_merge.hicup.sam.bam.chinput PANC_4F

Rscript /data/Brown_lab/chicago/region_capture_c/Chicago_package_scriptFiles/runChicago.R --settings-file /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/script/settingsFile --export-format seqMonk,interBed,washU_text,washU_track --design-dir /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/design_file_4F /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/MIA_T1_R1_2_merge.hicup.sam.bam/MIA_T1_R1_2_merge.hicup.sam.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/MIA_T2_R1_2_merge.hicup.sam.bam/MIA_T2_R1_2_merge.hicup.sam.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/MIA_T3_R1_2_merge.hicup.sam.bam/MIA_T3_R1_2_merge.hicup.sam.bam.chinput MIA_4F

Rscript /data/Brown_lab/chicago/region_capture_c/Chicago_package_scriptFiles/runChicago.R --settings-file /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/script/settingsFile --export-format seqMonk,interBed,washU_text,washU_track --design-dir /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/design_file_4F /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/Hs766T_Rep1_R1_2_merge.hicup.sam.bam/Hs766T_Rep1_R1_2_merge.hicup.sam.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/Hs766T_Rep2_R1_2_merge.hicup.sam.bam/Hs766T_Rep2_R1_2_merge.hicup.sam.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/Hs766T_Rep3_S1_L008_R1_2_001.hicup.sam.bam/Hs766T_Rep3_S1_L008_R1_2_001.hicup.sam.bam.chinput Hs766T_4F

Rscript /data/Brown_lab/chicago/region_capture_c/Chicago_package_scriptFiles/runChicago.R --settings-file /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/script/settingsFile --export-format seqMonk,interBed,washU_text,washU_track --design-dir /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/design_file_4F /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/HPDE_Rep1_R1_2_merge.hicup.sam.bam/HPDE_Rep1_R1_2_merge.hicup.sam.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/HPDE_Rep2_R1_2_merge.hicup.sam.bam/HPDE_Rep2_R1_2_merge.hicup.sam.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/HPDE_Rep3_S1_L007_R1_2_001.hicup.sam.bam/HPDE_Rep3_S1_L007_R1_2_001.hicup.sam.bam.chinput HPDE_4F

Rscript /data/Brown_lab/chicago/region_capture_c/Chicago_package_scriptFiles/runChicago.R --settings-file /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/script/settingsFile --export-format seqMonk,interBed,washU_text,washU_track --design-dir /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/design_file_4F /data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/AsPc1_Rep1_R1_2_merge.hicup.sam.bam/AsPc1_Rep1_R1_2_merge.hicup.sam.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/AsPc1_Rep2_R1_2_merge.hicup.sam.bam/AsPc1_Rep2_R1_2_merge.hicup.sam.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/four_frag_analysis/bam2chicago/AsPc1_Rep3_S1_L005_R1_2_001.hicup.sam.bam/AsPc1_Rep3_S1_L005_R1_2_001.hicup.sam.bam.chinput AsPc1_4F

