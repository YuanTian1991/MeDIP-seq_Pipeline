# This is a script I wrote to do preprocess on this medip-seq data. More specificlly, here I plan to do fastqc checking, fastp triming, and bowtie2 alignment, BlackList region     remove, and MACS2 Peak calling.
# Auhtor: Tian

directory <- "../Data/"
reference <- "./Genome/mm10"
threads <- 30


if (!file.exists("./myLog")) dir.create("./myLog")

message("\n[ Section 1: fastQC ] (require fastqc installed on server)")
if (!file.exists("./myFastQC")) dir.create("./myFastQC")
cmd <- paste0("fastqc --threads " , threads, "  --outdir ./myFastQC " , directory ,"* &> ./myLog/myFastQC.log")
message(cmd, '\n')
system(cmd)


message("\n[ Section 2: fastp ] (require parallel, fastp installed on server)")
if (!file.exists("./myFastp")) dir.create("./myFastp")
cmd <- paste0("parallel --plus 'fastp -h ./myFastp/{/..}.html -j ./myFastp/{/..}.json -i {} -o ./myFastp/{/..}.fastp.fq' ::: ", directory, "* &> ./myLog/myFastp.log")
message(cmd, '\n')
system(cmd)


message("\n[ Section 3: bowtie2 alignment ] (require bowtie2 installed on server)")
message("!!! Prepare genome from bowtie2 into a folder called Genome in this folder yourself, unzip it.")
if (!file.exists("./myAlignment")) dir.create("./myAlignment")

for(i in dir("./myFastp/", pattern="*.fastp.fq"))
{
    name <- strsplit(i, split="[.]")[[1]][1]
    cmd <- paste0("bowtie2 -p ",threads," -q --local -x ", reference ," -U  ./myFastp/", i , " | samtools view -bS - > ./myAlignment/", name,".bam")
    message(cmd)
    system(cmd)
}


message("\n[ Section 4: Remove Blacklist Regions ] (require GreylistChIP installed on server)")
message("In this step, I use GreyListChIP to generate Greylist from all Input sample, them merge them into one big file, then removed these regions across all samples.")
if (!file.exists("./myGreyList")) dir.create("./myGreyList")

library("GreyListChIP")
library("BSgenome.Mmusculus.UCSC.mm10")

files <- dir("./myAlignment")
files <- unique(sapply(files, function(x) strsplit(x, split="[.]")[[1]][1]))
Inps <- files[grep("Inp",files)]

for(i in Inps)
{
    message(i)
    gl <- greyListBS(BSgenome.Mmusculus.UCSC.mm10, paste0("./myAlignment/",i,".bam"))
    export(gl,con=paste0("./myGreyList/",i,"GreyList.bed"))
}

cmd <- "cat ./myGreyList/* > ./MergedGreyList.bed"
system(cmd)


cmd <- "parallel --plus 'bedtools intersect -v -abam {} -b ./MergedGreyList.bed > ./myGreyList/{/.}.grey_filtered.bam' ::: ./myAlignment/*.bam"
message(cmd)
system(cmd)

message("\n[ Section 5: Generate Bigwig file ] (require parallel, bamCoverage installed on server)")
if (!file.exists("./myBigWig")) dir.create("./myBigWig")
message("Sorting blacklist filtered bam files...")
 cmd <- "parallel --plus 'samtools sort {} -o ./myBigWig/{/.}.sorted.bam' ::: ./myGreyList/*.grey_filtered.bam"
# cmd <- "parallel --plus 'samtools sort {} -o ./myBigWig/{/.}.sorted.bam' ::: ./myAlignment/*.bam"
system(cmd)
cmd <- "parallel --plus 'samtools index {} ./myBigWig/{/..}.sorted.bam.bai' ::: ./myBigWig/*.bam"
system(cmd)
cmd <- "parallel --plus 'bamCoverage -p 5 -b {} -o ./myBigWig/{/.}.bw' ::: ./myBigWig/*.sorted.bam"
system(cmd)

