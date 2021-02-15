# This is a script I wrote to draw TFEA Enrichment Heatmap
# Author: Tian

load("./MouseDHMR_PromoterAnnotated.rda")
load("../3.Diffbind/BindingAffinityMatrix.rda")
source("./champ.getTFPeaks.R")
source("./champ.TFEA.R")

BindingAffinityMatrix <- as.data.frame(BindingAffinityMatrix)
BindingAffinityMatrix <- BindingAffinityMatrix[match(df$start, BindingAffinityMatrix$start),]
beta <- BindingAffinityMatrix[, 6:15]

# myRandomRegion <- df[df$insideFeature %in% c("inside", "overlapStart", "overlapEnd"),]
message("User all 5hmC peaks (note not all are significant here) as random peaks.")
myRandomRegion <- df

# index <- which(df$status %in% c("Mid(Stem)_Higher", "Mig(Stem)_Unique") & df$insideFeature == "inside") # Get Higher/Unique Promoter-Related for StemCell
message("Prepare top 500 Stem_Higher DhMR and their BindingAffinityMatrix")
Stem.Index <- which(df$status %in% c("Mid(Stem)_Higher") & df$insideFeature == "inside")[1:500] # Get Higher/Unique Promoter-Related for StemCell
StemPromoterROI <- df[Stem.Index, ]
betaStemMatrix <- beta[Stem.Index, ]
betaStemMatrix[betaStemMatrix == 0] <- min(betaStemMatrix[betaStemMatrix != 0])

message("Prepare top 500 Diff_Higher DhMR and their BindingAffinityMatrix")
Diff.Index <- which(df$status %in% c("Neg(Diff)_Higher") & df$insideFeature == "inside")[1:500] # Get Higher/Unique Promoter-Related for DiffCell
DiffPromoterROI <- df[Diff.Index, ]
betaDiffMatrix <- beta[Diff.Index, ]
betaDiffMatrix[betaDiffMatrix == 0] <- min(betaDiffMatrix[betaDiffMatrix != 0])

# possible params: 2 beta matrix, 2 peak (CpG) list, 1 peakList, top n peak to plot.

StemTFPeaks <- champ.getTFPeaks("mouse", "stem_cell")
set.seed(12345)
StemTFEAonStemROI <- champ.TFEA(StemPromoterROI, myRandomRegion, StemTFPeaks)
StemTFEAonDiffROI <- champ.TFEA(DiffPromoterROI, myRandomRegion, StemTFPeaks)

DiffTFPeaks <- champ.getTFPeaks("mouse", "intestine")
set.seed(12345)
DiffTFEAonStemROI <- champ.TFEA(StemPromoterROI, myRandomRegion, DiffTFPeaks)
DiffTFEAonDiffROI <- champ.TFEA(DiffPromoterROI, myRandomRegion, DiffTFPeaks)

keyTFsID <- StemTFEAonStemROI$TFList[1:10, "ID"]
keyDiffTFsID <- DiffTFEAonDiffROI$TFList[1:2, "ID"]

StemTF <- sapply(keyTFsID, function(x) 1:nrow(StemPromoterROI) %in% StemTFEAonStemROI$PeakOV[[x]]$ovROI$ROIHits)
DiffTF <- sapply(keyTFsID, function(x) 1:nrow(DiffPromoterROI) %in% StemTFEAonDiffROI$PeakOV[[x]]$ovROI$ROIHits)
StemTF[StemTF == "FALSE"] <- DiffTF[DiffTF == "FALSE"] <- 0
StemTF[StemTF == "TRUE"] <- DiffTF[DiffTF == "TRUE"] <- 1

StemDiffTF <- sapply(keyDiffTFsID, function(x) 1:nrow(StemPromoterROI) %in% DiffTFEAonStemROI$PeakOV[[x]]$ovROI$ROIHits)
DiffDiffTF <- sapply(keyDiffTFsID, function(x) 1:nrow(DiffPromoterROI) %in% DiffTFEAonDiffROI$PeakOV[[x]]$ovROI$ROIHits)
StemDiffTF[StemDiffTF == "FALSE"] <- DiffDiffTF[DiffDiffTF == "FALSE"] <- 0
StemDiffTF[StemDiffTF == "TRUE"] <- DiffDiffTF[DiffDiffTF == "TRUE"] <- 1

library("pheatmap")
A <- pheatmap::pheatmap(betaStemMatrix)
StemOrder <- A$tree_row$order

B <- pheatmap::pheatmap(betaDiffMatrix)
DiffOrder <- B$tree_row$order

C <- pheatmap::pheatmap(rbind(betaStemMatrix, betaDiffMatrix))
ColOrder <- C$tree_col$order

MergeBeta <- log(rbind(betaStemMatrix[StemOrder, ColOrder], betaDiffMatrix[DiffOrder, ColOrder]))
MergeTF <- cbind(rbind(StemTF[StemOrder, ], DiffTF[DiffOrder, ]),
                 rbind(StemDiffTF[StemOrder, ], DiffDiffTF[DiffOrder, ]))

#pheatmap(MergeBeta, cluster_rows=F, cluster_cols=F, show_rownames=F, gaps_row=c(500))
colnames(mat) = paste0("column", seq_len(nc))
LeftAnno <- c("Stem_Higher_Promoter DhmR", "Diff_Higher_Promoter DhmR") 
library("ComplexHeatmap")
hit1 <- Heatmap(MergeBeta, 
        name = "MergeBeta", 
        cluster_rows = FALSE, 
        cluster_columns=FALSE,
        show_row_names=FALSE,
        row_split = factor(rep(LeftAnno, each=500), levels=LeftAnno))

hitTF <- list()
hit_list <- hit1
for(i in 1:10)
{
    hitTF[[i]] <- Heatmap(MergeTF[,i], 
                    name = substr(colnames(MergeTF)[i],20,100),
                    col =  colorRamp2(c(0, 1), c("white", "#1565c0")),
                    cluster_rows = FALSE, 
                    cluster_columns=FALSE,
                    show_heatmap_legend=FALSE,
                    show_row_names=FALSE)
    hit_list <- hit_list + hitTF[[i]]

}

for(i in 11:12)
{
    hitTF[[i]] <- Heatmap(MergeTF[,i], 
                    name = substr(colnames(MergeTF)[i],20,100),
                    col =  colorRamp2(c(0, 1), c("white", "#fb8c00")),
                    cluster_rows = FALSE, 
                    cluster_columns=FALSE,
                    show_heatmap_legend=FALSE,
                    show_row_names=FALSE)
    hit_list <- hit_list + hitTF[[i]]

}
