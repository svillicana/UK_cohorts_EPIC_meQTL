## read table

setwd("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR")
dat=data.table::fread("trans_meQTL_full.txt",header=T,data.table=F)

## create SNP location file

snps.info<-unique(dat[,c("SNP", "REF", "ALT", "ALT_AF")])
snps<-unlist(lapply(strsplit(snps.info$SNP, "_"), head, n = 1))
chr<-unlist(lapply(strsplit(snps, ":"), head, n = 1))
bp<-unlist(lapply(strsplit(snps, ":"), tail, n = 1))

geno.map<-data.frame("Chr" = chr, "id"= snps, "cm" = 0, "bp" = bp, "A1" = snps.info$ALT, "A2" = snps.info$REF, "Freq" = snps.info$ALT_AF, stringsAsFactors = FALSE)

geno.map<-geno.map[order(as.numeric(geno.map[,1]), as.numeric(geno.map[,4])),]

## create list of duplicate SNP names i.e those at same position to exclude

dups<-unique(snps[duplicated(snps)])
geno.map<-geno.map[!geno.map[,2] %in% dups,]

write.table(geno.map, "FormattedSMRFiles/Blood_trans.esi", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


## create DNA methylation site file
annoObj <-  minfi::getAnnotationObject("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
all <- minfi:::.availableAnnotation(annoObj)$defaults
newfData <- do.call(cbind, lapply(all, function(wh) {
        minfi:::.annoGet(wh, envir = annoObj@data)
}))

newfData<-newfData[unique(dat$cpg),]


probes.tmp<-data.frame("Chr" = newfData$chr, "name" = newfData$Name,"cm" = 0, "pos" = newfData$pos, "Gene" = newfData$UCSC_RefGene_Name, "strand" = "+", stringsAsFactors = FALSE)
probes.tmp[,1]<-as.numeric(gsub("chr", "", probes.tmp[,1]))

#probes.tmp[,5]<-as.character(probes.tmp[,5])
probes.tmp[which(probes.tmp[,5] == ""),5]<-"NA"
probes.tmp<-probes.tmp[order(probes.tmp[,1], probes.tmp[,4]),]
write.table(probes.tmp, "FormattedSMRFiles/Blood_trans.epi", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

### create a file for every single probe and a reference file
for(each in unique(dat$cpg)){
	sub<-dat[which(dat$cpg == each),c("SNP", "ALT", "REF", "ALT_AF", "beta", "se", "p-value")]
	snps<-unlist(lapply(strsplit(sub$SNP, "_"), head, n = 1))
	chr<-unlist(lapply(strsplit(snps, ":"), head, n = 1))
	bp<-unlist(lapply(strsplit(snps, ":"), tail, n = 1))
	res<-cbind(chr, snps, bp, sub[,-1])
	colnames(res)<-c("Chr","SNP","Bp","A1","A2","Freq","Beta","se","p")
	## exclude duplicate SNPs IDs
	res<-res[!res$SNP %in% dups,]
	write.table(res, paste("mQTLResults_trans/", each, ".esd", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
}

probes<-unique(dat$cpg)
fileRef<-cbind(probes.tmp[match(probes, probes.tmp$name),], paste("mQTLResults_trans/", probes, ".esd", sep = ""))
colnames(fileRef)<-c("Chr","ProbeID","GeneticDistance","ProbeBp","Gene","Orientation","PathOfEsd")
fileRef$Gene<-as.character(fileRef$Gene)
fileRef$Gene[which(fileRef$Gene == "")]<-"NA"
write.table(fileRef, "FormattedSMRFiles/Blood_trans.flist", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)	
