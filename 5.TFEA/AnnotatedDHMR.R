library(GenomicRanges)

load("../../3.Diffbind/peaks_lib.rda")
myPeaks <- myDBA.peaks
threshold <- 0.01

message("Annotate Regions")
library(org.Mm.eg.db)
library(ChIPpeakAnno)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
annoData <- toGRanges(TxDb.Mmusculus.UCSC.mm10.knownGene, feature="gene")
annoPromoter <- GenomicRanges::promoters(annoData)
myPeaks <- annotatePeakInBatch(myPeaks, AnnotationData=annoPromoter)
myPeaks <- addGeneIDs(myPeaks, "org.Mm.eg.db", IDs2Add = c("entrez_id",'ensembl_gene_id', 'symbol'), feature_id_type="entrez_id")

df <- as.data.frame(myPeaks)
df$status <- "No-Sig"
df$status[df$FDR <= threshold & df$Called2 >= 3 & df$Called1 >= 3 & df$Fold > 0] <- "Mid(Stem)_Higher"
df$status[df$FDR <= threshold & df$Called2 >= 3 & df$Called1 >= 3 & df$Fold < 0] <- "Neg(Diff)_Higher"
df$status[df$FDR <= threshold & df$Called2 >=3 & df$Called1 < 3 ] <- "Mig(Stem)_Unique"
df$status[df$FDR <= threshold & df$Called2 < 3 & df$Called1 >= 3 ] <- "Neg(Diff)_Unique"
