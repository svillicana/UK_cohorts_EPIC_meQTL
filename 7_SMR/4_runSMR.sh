#!/bin/sh -l
#SBATCH --job-name=runSMR
#SBATCH --output=logs/runSMR_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=40GB
#SBATCH --ntasks=10
#SBATCH --time=1:0:0

cd /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR

### use R script to create files

module load use.own
module load smr

### use SMR software to process before running analysis per GWAS trait
### converts mQTL results into required format
smr_Linux --eqtl-flist Blood.flist --make-besd --out Blood 

# Especify reference genome
ref_genome=/users/k19063501/scratch/k19063501/Reference/1000G/P3/EUR_1000G_LiftedOver_hg19

# Run analysis for each trait
# Replication was tested with
# smr_Linux --bfile ${ref_genome} --gwas-summary reference_files/jointGwasMc_LDL_SMRFormat.txt --beqtl-summary FormattedSMRFiles/Blood --out Results/SMROUT/LDL_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/ALS/2019/meta.all.saige.stats.SMR.txt --beqtl-summary Blood --out SMROUT/ALS_2019 --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/Schizophrenia/GWAS/PGC2+CLOZUK/CLOZUK_PGC2noclo.wCHRX.w1000Gfrq.wNtot.METAL.assoc.dosage_SMRFormat.txt  --beqtl-summary Blood --out SMROUT/SCZ_CLOZUK_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/BMI/GIANT_BMI_Locke_et_al_European_Ancestry_SMR.txt --beqtl-summary Blood --out SMROUT/BMI_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/Height/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq_SMR.txt --beqtl-summary Blood --out SMROUT/Height_GWAS --thread-num 10 --diff-freq-prop 0.1;

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/WHRadjBMI/GIANT_2015_WHRadjBMI_COMBINED_EUR_SMR.txt --beqtl-summary Blood --out SMROUT/WHRadjBMI_GWAS --thread-num 10 --diff-freq-prop 0.1;

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/RA/RA_GWASmeta_European_v2_SMRFormat.txt --beqtl-summary Blood --out SMROUT/RA_GWAS --thread-num 10 --diff-freq-prop 0.1;

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/Alzheimers/IGAP_stage_1_SMRFormat.txt --beqtl-summary Blood --out SMROUT/AD_Stage1_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/Cardiogram/CARDIoGRAM_GWAS_RESULTS_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Cardiogram_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/Cardiogram/cad.add.160614.website_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Cardiogram_2_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/EGG_BW2_DISCOVERY_SMRFormat.txt --beqtl-summary Blood --out SMROUT/EGG_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/T2Diabetes/GWAS/DIAGRAM_Gaulton_2015_SMRFormat.txt --beqtl-summary Blood --out SMROUT/T2Diabetes_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/T2Diabetes/GWAS/DIAGRAM.website.GWAS.metabochip_SMRFormat.txt --beqtl-summary Blood --out SMROUT/T2Diabetes_Metabochip_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/T2Diabetes/GWAS/DIAGRAMv3.2012DEC17_SMRFormat.txt --beqtl-summary Blood --out SMROUT/T2Diabetes_2012DEC17_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/IBD/iibdgc-trans-ancestry-summary-stats/EUR.CD.gwas.assoc_SMRFormat.txt --beqtl-summary Blood --out SMROUT/IBD_CD_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/IBD/iibdgc-trans-ancestry-summary-stats/EUR.UC.gwas.assoc_SMRFormat.txt --beqtl-summary Blood --out SMROUT/IBD_UC_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/IBD/iibdgc-trans-ancestry-summary-stats/EUR.IBD.gwas.assoc_SMRFormat.txt --beqtl-summary Blood --out SMROUT/IBD_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/Migraine/any_mig.gwama_.out_.isq75.nstud12.clean_.p1e-5_3_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Migraine_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/T2Diabetes/GWAS/DIAGRAM_Gaulton_2015_SMRFormat.txt --beqtl-summary Blood --out SMROUT/T2Diabetes_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/TAG/tag.evrsmk.tbl_SMRFormat.txt --beqtl-summary Blood --out SMROUT/EverSmoked_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/TAG/tag.cpd.tbl_SMRFormat.txt --beqtl-summary Blood --out SMROUT/CPD_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/BipolarDisorder/pgc.bip.full.2012-04_SMRFormat.txt --beqtl-summary Blood --out SMROUT/BPD_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/MDD/pgc.mdd.full.2012-04_SMRFormat.txt --beqtl-summary Blood --out SMROUT/MDD_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/SSGAC/SSGAC_College_Rietveld2013_publicrelease_SMRFormat.txt --beqtl-summary Blood --out SMROUT/CollegeCompletion_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/SSGAC/SSGAC_EduYears_Rietveld2013_publicrelease_SMRFormat.txt --beqtl-summary Blood --out SMROUT/EduYears_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/SSGAC/Neuroticism_Full_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Neuroticism_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/SSGAC/DS_Full_SMRFormat.txt --beqtl-summary Blood --out SMROUT/DS_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/SSGAC/SWB_Full_SMRFormat.txt --beqtl-summary Blood --out SMROUT/SWB_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/Reprogen/Menopause_HapMap2_DayNG2015_18112015_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Menopause_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GlobalLipidsConsortium/jointGwasMc_HDL_SMRFormat.txt --beqtl-summary Blood --out SMROUT/HDL_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GlobalLipidsConsortium/jointGwasMc_LDL_SMRFormat.txt --beqtl-summary Blood --out SMROUT/LDL_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GlobalLipidsConsortium/jointGwasMc_TC_SMRFormat.txt --beqtl-summary Blood --out SMROUT/TC_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GlobalLipidsConsortium/jointGwasMc_TG_SMRFormat.txt --beqtl-summary Blood --out SMROUT/TG_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/Pubertal_growth_10F_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Pubertal_growth_10F_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/Pubertal_growth_12M_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Pubertal_growth_12M_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/Pubertal_growth_10F_12M_combined_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Pubertal_growth_10F_12M_combined_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/Pubertal_growth_PGF_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Pubertal_growth_PGF_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/Pubertal_growth_PGM_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Pubertal_growth_PGM_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/Pubertal_growth_PGF_PGM_combined_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Pubertal_growth_PGF_PGM_combined_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/Pubertal_growth_PTF_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Pubertal_growth_PTF_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/Pubertal_growth_PTM_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Pubertal_growth_PTM_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/Pubertal_growth_PTF_PTM_combined_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Pubertal_growth_PTF_PTM_combined_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/EGG-GWAS-BL_SMRFormat.txt --beqtl-summary Blood --out SMROUT/EGG_BL_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/EGG_TANNER_females_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Tanner_females_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/EGG_TANNER_males_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Tanner_males_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/EGG/EGG_TANNER_males_females_combined_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Tanner_male_females_combined_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/GIANT_OVERWEIGHT_Stage1_Berndt2013_publicrelease_HapMapCeuFreq_SMRFormat.txt --beqtl-summary Blood --out SMROUT/Overweight_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/GIANT_OBESITY_CLASS1_Stage1_Berndt2013_publicrelease_HapMapCeuFreq_SMRFormat.txt --beqtl-summary Blood --out SMROUT/ObesityClass1_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/GIANT_OBESITY_CLASS2_Stage1_Berndt2013_publicrelease_HapMapCeuFreq_SMRFormat.txt --beqtl-summary Blood --out SMROUT/ObesityClass2_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/GIANT_OBESITY_CLASS3_Stage1_Berndt2013_publicrelease_HapMapCeuFreq_SMRFormat.txt --beqtl-summary Blood --out SMROUT/ObesityClass3_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/GIANT_EXTREME_BMI_Stage1_Berndt2013_publicrelease_HapMapCeuFreq_SMRFormat.txt --beqtl-summary Blood --out SMROUT/EXTREME_BMI_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/GIANT_EXTREME_HEIGHT_Stage1_Berndt2013_publicrelease_HapMapCeuFreq_SMRFormat.txt --beqtl-summary Blood --out SMROUT/EXTREME_HEIGHT_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/GIANT_EXTREME_WHR_Stage1_Berndt2013_publicrelease_HapMapCeuFreq_SMRFormat.txt --beqtl-summary Blood --out SMROUT/EXTREME_WHR_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/GIANT_2015_HIP_COMBINED_EUR_SMRFormat.txt --beqtl-summary Blood --out SMROUT/HIP_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/GIANT_2015_WC_COMBINED_EUR_SMRFormat.txt --beqtl-summary Blood --out SMROUT/WC_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/GIANT_2015_WCadjBMI_COMBINED_EUR_SMRFormat.txt --beqtl-summary Blood --out SMROUT/WCadjBMI_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/GIANT/GIANT_2015_WHR_COMBINED_EUR_SMRFormat.txt --beqtl-summary Blood --out SMROUT/WHR_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/SSGAC/EduYears_Main_SMRFormat.txt --beqtl-summary Blood --out SMROUT/EduYearsMain_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/Autism/daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3_SMR.txt --beqtl-summary Blood --out SMROUT/Autism_lmm_GWAS --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary /mnt/data1/reference_files/ADHD/adhd_jul2017_SMRFormat.txt --beqtl-summary Blood --out SMROUT/ADHD --thread-num 10 --diff-freq-prop 0.1;
