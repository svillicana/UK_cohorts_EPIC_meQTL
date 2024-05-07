library(data.table)
library(dplyr)
library(GenomicRanges)

# Load annotation data --------------------------------------------------------

path <- "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

TADs_files <- list.files(paste0(path,"TADs/hg19/"), full.names = TRUE)

TADs_l <- lapply(TADs_files, fread)
TADs_Rao <- TADs_l[[2]]
TADs_all <- rbindlist(TADs_l)

# Prepare annotation data -----------------------------------------------------

TADs_Rao[,V1:=gsub('chr','',V1)]
TADs_all[,V1:=gsub('chr','',V1)]

TADs_Rao_gr <- makeGRangesFromDataFrame(TADs_Rao, seqnames.field = "V1",
                                        start.field = "V2", end.field = "V3",
                                        ignore.strand = TRUE, starts.in.df.are.0based = TRUE)
TADs_all_gr <- makeGRangesFromDataFrame(TADs_all, seqnames.field = "V1",
                                        start.field = "V2", end.field = "V3",
                                        ignore.strand = TRUE, starts.in.df.are.0based = TRUE)

# Overlaps for cis  -----------------------------------------------------------

# Read data
cis <- fread(paste0(path,"meQTL_annotation/cis_full_location_MAF.txt"),
             select = c("CpG", "SNP", "chr_cpg", "pos_cpg", "chr_snp", "pos_snp"),
             header = TRUE)

# Prepare data
cis <- cis[abs(pos_cpg-pos_snp) >= 100000,] # Filter associations with distance > 100 kbp (min distance between loops)


cis_cpg_gr <- makeGRangesListFromDataFrame(select(cis, chr_cpg, pos_cpg, CpG),
                                           seqnames.field = "chr_cpg", start.field = "pos_cpg", end.field = "pos_cpg",
                                           ignore.strand = TRUE, keep.extra.columns = TRUE)
cis_snp_gr <- makeGRangesListFromDataFrame(select(cis, chr_snp, pos_snp, SNP),
                                           seqnames.field = "chr_snp", start.field = "pos_snp", end.field = "pos_snp",
                                           ignore.strand = TRUE, keep.extra.columns = TRUE)
cis_grl <- pc(cis_cpg_gr, cis_snp_gr)

# Overlaps
TADs_Rao_count <- countOverlaps(cis_grl, TADs_Rao_gr, type = "within")
TADs_all_count <- countOverlaps(cis_grl, TADs_all_gr, type = "within")

TADs_Rao_i <- findOverlaps(cis_grl, TADs_Rao_gr, type = "within")
TADs_all_i <- findOverlaps(cis_grl, TADs_all_gr, type = "within")

TADs_Rao_cis <- TADs_Rao[unique(subjectHits(TADs_Rao_i)),]
TADs_all_cis <- TADs_all[unique(subjectHits(TADs_all_i)),]

# Save
cis <- cis[,c("TADs_Rao","TADs_all"):=list(TADs_Rao_count,TADs_all_count)]
cis <- cis[TADs_Rao_count+TADs_all_count>0,]

fwrite(cis, paste0(path,"meQTL_annotation/cis_TADs.txt"), sep = "\t")

# Overlaps for trans  ---------------------------------------------------------

# Read data
trans <- fread(paste0(path,"meQTL_annotation/trans_full_location_MAF.txt"),
               select = c("CpG", "SNP", "chr_cpg", "pos_cpg", "chr_snp", "pos_snp"),
               header = TRUE)

# Prepare data
trans <- trans[chr_cpg==chr_snp,] # Filter intra associations


trans_cpg_gr <- makeGRangesListFromDataFrame(select(trans, chr_cpg, pos_cpg, CpG),
                                             seqnames.field = "chr_cpg", start.field = "pos_cpg", end.field = "pos_cpg",
                                             ignore.strand = TRUE, keep.extra.columns = TRUE)
trans_snp_gr <- makeGRangesListFromDataFrame(select(trans, chr_snp, pos_snp, SNP),
                                             seqnames.field = "chr_snp", start.field = "pos_snp", end.field = "pos_snp",
                                             ignore.strand = TRUE, keep.extra.columns = TRUE)
trans_grl <- pc(trans_cpg_gr, trans_snp_gr)

# Overlaps
TADs_Rao_count <- countOverlaps(trans_grl, TADs_Rao_gr, type = "within")
TADs_all_count <- countOverlaps(trans_grl, TADs_all_gr, type = "within")

TADs_Rao_i <- findOverlaps(trans_grl, TADs_Rao_gr, type = "within")
TADs_all_i <- findOverlaps(trans_grl, TADs_all_gr, type = "within")

TADs_Rao_trans <- TADs_Rao[unique(subjectHits(TADs_Rao_i)),]
TADs_all_trans <- TADs_all[unique(subjectHits(TADs_all_i)),]

# Save
trans <- trans[,c("TADs_Rao","TADs_all"):=list(TADs_Rao_count,TADs_all_count)]
trans <- trans[TADs_Rao_count+TADs_all_count>0,]

fwrite(trans, paste0(path,"meQTL_annotation/trans_TADs.txt"), sep = "\t")

# T-test ----------------------------------------------------------------------

# Size of TADs
TADs_all_cis <- TADs_all_cis[,size:=V3-V2]
TADs_all_trans <- TADs_all_trans[,size:=V3-V2]
TADs_Rao_cis <- TADs_Rao_cis[,size:=V3-V2]
TADs_Rao_trans <- TADs_Rao_trans[,size:=V3-V2]

# Test
t.test(TADs_Rao_cis[,size], TADs_Rao_trans[,size], var.equal=FALSE)
t.test(TADs_all_cis[,size], TADs_all_trans[,size], var.equal=FALSE)


# Session info ----------------------------------------------------------------
sessionInfo()

