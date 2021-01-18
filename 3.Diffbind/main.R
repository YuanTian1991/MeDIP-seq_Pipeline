# A script to use Diffbind to get Differentially 5hmC Enriched Genes
# Below code is designed for DiffBind 3.0.3
# Author: Tian

library("DiffBind")
message("Prepare a Sample Sheet for loading all file")
Samples <- read.csv("../2.PeakCalling/SampleSheet.csv", header=T, as.is=T)
Samples <- Samples[!is.na(Samples$InputToUse),]
Samples$Factor <- substr(Samples$Factor, 14,16)

SampleList <- data.frame(SampleID=paste(Samples$SampleName, Samples$Factor, sep="_"),
                         Tissue=Samples$Tissue,
                         Factor="",
                         Condition=Samples$Factor,
                         Treatment="5hmC",
                         Replicate=c(1,2,1,2,3,3,4,5,4,5),
                         bamReads=paste0("../1.Preprocess/myBigWig/", Samples$SampleName, ".grey_filtered.sorted.bam"),
                         ControlID="",
                         bamControl=paste0("../1.Preprocess/myBigWig/", Samples$InputToUse, ".grey_filtered.sorted.bam"),
                         Peaks=paste0("../2.PeakCalling/macs2/", Samples$SampleName , "_peaks.narrowPeak"),
                         PeakCaller="narrow")

# SampleList$Replicate <- c(1:5, 1:4, 1:5, 1:4)

write.csv(SampleList, file="SampleList.csv",quote=F, row.names=F)

message("Differential Analysis")
myDBA <- dba(sampleSheet="./SampleList.csv", dir="./")
myDBA <- dba.count(myDBA, summits=250, bParallel=40, minOverlap=0.33, bRemoveDuplicates=TRUE)

myDBA <- dba.normalize(myDBA, normalize=DBA_NORM_NATIVE) # This is the final choosen normalisation method. Below one is also not bad.
# myDBA <- dba.normalize(myDBA,  normalize=DBA_NORM_LIB, library=DBA_LIBSIZE_PEAKREADS, background=FALSE)

myDBA <- dba.contrast(myDBA, categories=DBA_CONDITION)
myDBA <- dba.analyze(myDBA, bParallel=40, bBlacklist=FALSE ,bGreylist=FALSE)
rownames(myDBA$binding) <- 1:nrow(myDBA$binding)

myDBA.peaks <- dba.report(myDBA, bCalled=TRUE, th=1)
BindingAffinityMatrix <- dba.peakset(myDBA, bRetrieve=TRUE)
dba.plotMA(myDBA, method=DBA_DESEQ2, sub="log consentration")

save(myDBA, myDBA.peaks, file="ResultDBA.rda")
save(myDBA.peaks, file="peaks_lib.rda")
save(BindingAffinityMatrix, file="BindingAffinityMatrix.rda")

