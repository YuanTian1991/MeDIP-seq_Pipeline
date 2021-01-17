# This script is written to do visialisation on DiffBind result. I am planning to show: Volcano Plot, barplot, karyoplote, top gene Heatmap.
# Author: Tian

load("../3.Diffbind/peaks_lib.rda")

myPeaks <- myDBA.peaks
threshold <- 0.001

message("Annotate Regions")
library(org.Mm.eg.db)
library(ChIPpeakAnno)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
annoData <- toGRanges(TxDb.Mmusculus.UCSC.mm10.knownGene, feature="gene")
myPeaks <- annotatePeakInBatch(myPeaks, AnnotationData=annoData)
myPeaks <- addGeneIDs(myPeaks, "org.Mm.eg.db", IDs2Add = c("entrez_id",'ensembl_gene_id', 'symbol'), feature_id_type="entrez_id")

message("Plotting VolcanoPlot")

library("ggplot2")
library("ggrepel")
    
df <- as.data.frame(myPeaks)

df$status <- "No-Sig"
df$status[df$FDR <= threshold & df$Called2 >= 3 & df$Called1 >= 3 & df$Fold > 0] <- "Mid(Stem)_Higher"
df$status[df$FDR <= threshold & df$Called2 >= 3 & df$Called1 >= 3 & df$Fold < 0] <- "Neg(Diff)_Higher"
df$status[df$FDR <= threshold & df$Called2 >=3 & df$Called1 < 3 ] <- "Mig(Stem)_Unique"
df$status[df$FDR <= threshold & df$Called2 < 3 & df$Called1 >= 3 ] <- "Neg(Diff)_Unique"

showGene <- c(which(df$FDR <= threshold & df$Called2 >= 3 & df$Called1 >= 3 & df$Fold > 0)[1:20],
              which(df$FDR <= threshold & df$Called2 >= 3 & df$Called1 >= 3 & df$Fold < 0)[1:20],
              which(df$FDR <= threshold & df$Called2 >=3 & df$Called1 < 3)[1:20],
              which(df$FDR <= threshold & df$Called2 < 3 & df$Called1 >= 3)[1:20])

peakList <- GRangesList(Stem_Higher=myPeaks[df$status == 'Mid(Stem)_Higher'],
                 Diff_Higher=myPeaks[df$status == 'Neg(Diff)_Higher'],
                 Stem_Unique=myPeaks[df$status == 'Mig(Stem)_Unique'],
                 Diff_Unique=myPeaks[df$status == 'Neg(Diff)_Unique'])

genomicElementDistribution(peakList, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, promoterRegion=c(upstream=2000, downstream=500),geneDownstream=c(upstream=0, downstream=2000))

    
p <- ggplot(df, aes(x = Fold, y = -log10(FDR), col = status)) + 
        geom_point(size = 0.5) + 
        scale_color_manual(values = c("#1565c0", '#bbdefb', "#fb8c00", '#ffe0b2', "#999999")) + 
        geom_hline(yintercept = c(-log10(0.001)), linetype = "dashed", color = "black") + 
        theme_minimal(base_size = 14) + 
        geom_text_repel(data = df[showGene, ], aes(label = symbol), 
            size = 3.5, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, 
                "lines"))

print(p)
#if (!file.exists("./Figure")) 
#        dir.create("./Figure")
#    pdf("./Figure/DMPVolcanoPlot.pdf", width = 10, height = 10)
#    print(p)
#    dev.off()
