# This is a script I wrote to do TFEA between TFregulomeR and DhMR from Mouse Intestine
# Author: Tian

load("./MouseDHMR_PromoterAnnotated.rda")
source("./champ.getTFPeaks.R")
source("./champ.TFEA.R")

# myRandomRegion <- df[df$insideFeature %in% c("inside", "overlapStart", "overlapEnd"),]
message("User all 5hmC peaks (note not all are significant here) as random peaks.")
myRandomRegion <- df

# index <- which(df$status %in% c("Mid(Stem)_Higher", "Mig(Stem)_Unique") & df$insideFeature == "inside") # Get Higher/Unique Promoter-Related for StemCell
message("Enricher Top 1000 Stem Higher/Unique DhMR with Stem Cell Peaks.")
index <- which(df$status %in% c("Mid(Stem)_Higher", "Mig(Stem)_Unique")) # Get Higher/Unique Promoter-Related for StemCell
StemROI <- df[index[1:1000], ]
myTFPeaks <- champ.getTFPeaks("mouse", "stem_cell")
set.seed(12345)
StemTFEA <- champ.TFEA(StemROI, myRandomRegion, myTFPeaks)

message("Enricher Top 1000 Diff Higher/Unique DhMR with Intestine Cell Peaks.")
index <- which(df$status %in% c("Neg(Diff)_Higher", "Neg(Diff)_Unique")) # Get Higher/Unique Promoter-Related for DiffCell
DiffROI <- df[index[1:1000], ]
myTFPeaks <- champ.getTFPeaks("mouse", "intestine")
set.seed(12345)
DiffTFEA <- champ.TFEA(DiffROI, myRandomRegion, myTFPeaks)


