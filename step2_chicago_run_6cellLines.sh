#!/bin/bash
#SBATCH --mem=240G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mail-type=FAIL,END,TIME_LIMIT_80

set -e
module load bedtools/2.25
module load perl/5.24
module load R/4.0.5

cd /data/Brown_lab/Laufey_lab_capture_C/one_frag_analysis
Rscript /data/Brown_lab/chicago/region_capture_c/Chicago_package_scriptFiles/runChicago.R --settings-file /data/Brown_lab/Laufey_lab_capture_C/one_frag_analysis/script/settingsFile --export-format seqMonk,interBed,washU_text,washU_track --design-dir /data/Brown_lab/Laufey_lab_capture_C/one_frag_analysis/design_file_1F /data/Brown_lab/Laufey_lab_capture_C/one_frag_analysis/bam2chicago_input/AsPc1.bam/AsPc1.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/one_frag_analysis/bam2chicago_input/Hs766T.bam/Hs766T.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/one_frag_analysis/bam2chicago_input/MIA.bam/MIA.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/one_frag_analysis/bam2chicago_input/Panc.bam/Panc.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/one_frag_analysis/bam2chicago_input/SU8686.bam/SU8686.bam.chinput,/data/Brown_lab/Laufey_lab_capture_C/one_frag_analysis/bam2chicago_input/HPDE.bam/HPDE.bam.chinput All_6_cell_lines_1F
