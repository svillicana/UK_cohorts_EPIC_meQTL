library(data.table)
library(dplyr)
library(GenomicRanges)

# Load annotation data --------------------------------------------------------

path <- "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

domains <- fread(paste0(path,"3D_map_Rao_et_al/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt"),
                 select = 1:3, header = TRUE, data.table = FALSE)

subcompartments <- fread(paste0(path,"3D_map_Rao_et_al/GSE63525_GM12878_subcompartments.bed"),
                         select = 1:4, header = FALSE, data.table = FALSE)

loops <- fread(paste0(path,"3D_map_Rao_et_al/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt"),
               select = c(1:3,5:6), header = TRUE, data.table = FALSE)

# Prepare annotation data -----------------------------------------------------

# Contact domains
domains_gr <- makeGRangesFromDataFrame(domains, seqnames.field = "chr1", start.field = "x1",
                                       end.field = "x2", ignore.strand = TRUE)

# Subcompartments 
subcompartments <- subcompartments %>%
  transmute(chr = gsub("chr", "", V1), start = V2,
            end = V3, subcompartment = V4)

subcompartments_gr <- makeGRangesFromDataFrame(subcompartments, ignore.strand = TRUE,
                                               keep.extra.columns = TRUE, starts.in.df.are.0based = T)
# Loops
loops_1_gr <- makeGRangesListFromDataFrame(loops, seqnames.field = "chr1", start.field = "x1",
                                           end.field = "x2", ignore.strand = TRUE)
loops_2_gr <- makeGRangesListFromDataFrame(loops, seqnames.field = "chr1", start.field = "y1",
                                           end.field = "y2", ignore.strand = TRUE)
loops_grl <- pc(loops_1_gr, loops_2_gr)

# Overlaps for cis  -----------------------------------------------------------

# Read data
cis <- fread(paste0(path,"meQTL_annotation/cis_full_location_MAF.txt"),
             select = c("CpG", "SNP", "chr_cpg", "pos_cpg", "chr_snp", "pos_snp"),
             header = TRUE, data.table = FALSE)

# Prepare data
cis <- cis %>%
  filter(abs(pos_cpg-pos_snp) >= 30000) # Filter associations with distance > 30 kbp (min distance between loops)

cis_cpg_gr <- makeGRangesListFromDataFrame(select(cis, chr_cpg, pos_cpg, CpG),
                                           seqnames.field = "chr_cpg", start.field = "pos_cpg", end.field = "pos_cpg",
                                           ignore.strand = TRUE, keep.extra.columns = TRUE)
cis_snp_gr <- makeGRangesListFromDataFrame(select(cis, chr_snp, pos_snp, SNP),
                                           seqnames.field = "chr_snp", start.field = "pos_snp", end.field = "pos_snp",
                                           ignore.strand = TRUE, keep.extra.columns = TRUE)
cis_grl <- pc(cis_cpg_gr, cis_snp_gr)

# Overlaps
domains_overlap_count <- countOverlaps(cis_grl, domains_gr, type = "within")
subcompartments_overlap_count <- countOverlaps(cis_grl, subcompartments_gr, type = "within")
loops_overlap_count <- countOverlaps(cis_grl, loops_grl, type = "within")

# Save
cis <- cis %>% 
  bind_cols(contact_domains = domains_overlap_count,
            subcompartments = subcompartments_overlap_count,
            loops = loops_overlap_count)

fwrite(cis, paste0(path,"meQTL_annotation/cis_30kbp_3D_map.txt"), sep = "\t")

rm(cis, cis_cpg_gr, cis_snp_gr, cis_grl)
gc()

# Overlaps for trans  ---------------------------------------------------------

# Read data
trans <- fread(paste0(path,"meQTL_annotation/trans_full_location_MAF.txt"),
               select = c("CpG", "SNP", "chr_cpg", "pos_cpg", "chr_snp", "pos_snp"),
               header = TRUE, data.table = FALSE)

# Prepare data
trans <- trans %>%
  filter(chr_cpg == chr_snp) # Filter intrachromosome associations

trans_cpg_gr <- makeGRangesListFromDataFrame(select(trans, chr_cpg, pos_cpg, CpG),
                                             seqnames.field = "chr_cpg", start.field = "pos_cpg", end.field = "pos_cpg",
                                             ignore.strand = TRUE, keep.extra.columns = TRUE)
trans_snp_gr <- makeGRangesListFromDataFrame(select(trans, chr_snp, pos_snp, SNP),
                                             seqnames.field = "chr_snp", start.field = "pos_snp", end.field = "pos_snp",
                                             ignore.strand = TRUE, keep.extra.columns = TRUE)
trans_grl <- pc(trans_cpg_gr, trans_snp_gr)

# Overlaps
domains_overlap_count <- countOverlaps(trans_grl, domains_gr, type = "within")
subcompartments_overlap_count <- countOverlaps(trans_grl, subcompartments_gr, type = "within")
loops_overlap_count <- countOverlaps(trans_grl, loops_grl, type = "within")

# Save
trans <- trans %>% 
  bind_cols(contact_domains = domains_overlap_count,
            subcompartments = subcompartments_overlap_count,
            loops = loops_overlap_count)

fwrite(trans, paste0(path,"meQTL_annotation/trans_intra_3D_map.txt"), sep = "\t")

# Session info ----------------------------------------------------------------
sessionInfo()

