#!/bin/bash -l
#SBATCH --job-name=annotate_fdr
#SBATCH --output=logs/annotate_fdr_%j
#SBATCH --partition gpu,cpu
#SBATCH --mem=128GB
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=1-0

module load r

echo "Filtering meQTLs for 450K legacy probes ---------------------------------------------"
Rscript --vanilla annotate_fdr_probes_450K_legacy.R

echo "Filtering meQTLs for EPIC only probes -----------------------------------------------"
Rscript --vanilla annotate_fdr_probes_epic_only.R

