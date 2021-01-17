# This is script to do MACS2 calling
# Author: Tian

bams <- read.csv("./SampleSheet.csv", header=T)
bams <- bams[!is.na(bams$InputToUse),]

if (!file.exists("./macs2")) dir.create("./macs2")

commands <- list()
for(i in 1:nrow(bams)) {
    
    command <- paste0("macs2 callpeak -t ../1.Preprocess/myGreyList/", bams[i,'SampleName'], ".grey_filtered.bam -c ../1.Preprocess/myGreyList/", 
                      bams[i,'InputToUse'],".grey_filtered.bam -f BAM -g hs --outdir macs2 -n ", bams[i,'SampleName'], " -B -q 0.1 2> macs2/", bams[i,'SampleName'], "-macs2.log")
    message(command)
    commands <- c(commands, command)
}

runMacs <- function(i)
{
    system(i)
}

library(doParallel)
detectCores()
cl <- makeCluster(length(commands))
registerDoParallel(cl)
getDoParWorkers()

library(foreach)
x <- foreach(i = commands) %dopar% runMacs(i)

registerDoSEQ()
on.exit(stopCluster(cl))
