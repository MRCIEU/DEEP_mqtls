library(data.table)

# Thanks to EasyQC harmonization
# This script is modified to harmonize indel alleles for afreq

args <- (commandArgs(TRUE));
afreq_file <- as.character(args[1]);
out_file <- as.character(args[2]);

afreq <- as.data.frame(fread(afreq_file))
afreq <- as.data.frame(fread("/user/work/er20212/bib_QC_phase1/processed_data/genetic_data/data.afreq"))

afreq[,2]<-as.character(afreq[,2])
afreq[,3]<-as.character(afreq[,3])
afreq[,4]<-as.character(afreq[,4])
SNPfail <-data.frame()
harmonize.alleles <- function(afreq=afreq,SNPfail=SNPfail) {

	message("Checking allele coding")
	a1<-afreq[,3]
	a2<-afreq[,4]
	SNPfail<-data.frame()
	### A1=NA & A2=NA
	isBothMissing <- which(is.na(a1) & is.na(a2))
	if(length(isBothMissing)>0) {SNPfail<-rbind(SNPfail,afreq[isBothMissing,2])}

	message("Allele harmonization:",length(isBothMissing)," alleles with NA are going to be removed")

	# Recode Sequence coding to D/I [noted]
	a1<-afreq[,3]
	a2<-afreq[,4]
	isRecode_seq1 <- which(nchar(a1)>nchar(a2))
	afreq[isRecode_seq1,3] <- "I"
	afreq[isRecode_seq1,4] <- "D"
	message("Allele harmonization:",length(isRecode_seq1)," alleles with sequence coding are recoded to REF: I and ALT: D")
	
	a1<-afreq[,3]
	a2<-afreq[,4]
	isRecode_seq2 <- which(nchar(a1)<nchar(a2))
	afreq[isRecode_seq2,3] <- "D"
	afreq[isRecode_seq2,4] <- "I"
	message("Allele harmonization:",length(isRecode_seq2)," alleles with sequence coding are recoded to REF: D and ALT: I")

	a1<-afreq[,3]
	a2<-afreq[,4]
	isRecode_seq3 <- which(nchar(a1)>1&nchar(a1)==nchar(a2))
	if(length(isRecode_seq3)>0) {SNPfail<-rbind(SNPfail,afreq[isRecode_seq3 ,2])}
	message("Allele harmonization:",length(isRecode_seq3)," alleles with sequence coding (REF and ALT same length) are removed")

	a1<-afreq[,3]
    a2<-afreq[,4]
    isInvalid <- !(a1%in%c("A","C","G","T","I","D")&a2%in%c("A","C","G","T","I","D")&a1!=a2)
	if(length(which(isInvalid))>0) {
	SNPfail<-rbind(SNPfail,afreq[which(isInvalid),2])}
    message("Allele harmonization:",length(which(isInvalid))," alleles with coding other than A,C,T,G,I,D are going to be removed")
	
	rm(a1,a2)
	
	recoded.afreq <- list(afreq, SNPfail)
	
	return(recoded.afreq)
}

recoded.afreq<-harmonize.alleles(afreq,SNPfail)

if(length(recoded.afreq[[2]])>0){
SNPfailures<-as.character(t(recoded.afreq[[2]]))}

write.table(recoded.afreq[[1]],out_file,sep="\t",quote=F,col.names=F,row.names=F)