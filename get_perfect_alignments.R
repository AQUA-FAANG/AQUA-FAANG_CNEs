#!/usr/bin/Rscript
#SBATCH -c 1
#SBATCH --mem 16G
#SBATCH -p daniocode
#SBATCH --array 1-4011%10

# adapted script for Nash's kurtosis method

library(tidyverse)
library(rtracklayer)
library(Hmisc)
library(GenomicRanges)
library(data.table)
library(stringr)
library(Biostrings)
library(CNEr)
library(plyranges)

###################

getLengthsOfIdenticalSeqs<-function(x){
  #x must be an Axt object imported by CNEr
  qr<-DNAStringSet(apply(as.data.frame(querySeqs(x)), 2, function(z){
    tmp<-paste0("M", z)
    out<-paste0(tmp, "N")
  }))
  tr<-DNAStringSet(apply(as.data.frame(targetSeqs(x)), 2, function(z){
    tmp<-paste0("K", z)
    out<-paste0(tmp, "V")
  }))
  
  identical_bases <- mapply(x=qr, y=tr, function(x,y){
    out<-which(as.raw(x)==as.raw(y))
  })
  
  relative_ranges <- lapply(identical_bases, function(x){
    IRanges::reduce(IRanges(x-2, x-2))
  })
  
  absolute_ranges <- mapply(z=relative_ranges, y=as.list(x), function(z,y){
    if(length(z) == 0){
      out<-list(data.frame())
    } else {
      out<-list(data.frame("first.seqnames"=seqnames(CNEr::first(y)),
                           "first.start"=(start(CNEr::first(y)) + start(z)),
                           "first.end"=(start(CNEr::first(y)) + end(z)),
                           "first.width"=(end(z) - start(z) + 1),
                           "first.strand"="*",
                           "second.seqnames"=seqnames(CNEr::second(y)),
                           "second.start"=(start(CNEr::second(y)) + start(z)),
                           "second.end"=(start(CNEr::second(y)) + end(z)),
                           "second.width"=(end(z) - start(z) + 1),
                           "second.strand"="*"))
    }
  })
  
  return(absolute_ranges)
}

idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))


####################################Start of Script #####################

idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# args<-commandArgs(T)
# 
species1<- "Salmon"   #args[1]
species2<-  "Trout"   #args[2]
# chrom<-args[3]

root_folders <- list.files(path = "/mnt/orca/projects/AQUA-FAANG/species/", 
                           recursive = T, pattern = "^GCA_\\d+\\.\\d$",
                           include.dirs = T, full.names = T)
names(root_folders) <-  c("Carp", "Zebrafish", "Seabass", "Pike", "Trout", 
                          "Salmon", "Turbot", "Seabream")


speciesDir<- root_folders[species1]
outDir<-paste0(speciesDir, "/vs", species2, "/identical_seqs/")

if(!dir.exists(file.path(speciesDir))){
  dir.create(file.path(speciesDir))
}


if(!dir.exists(file.path(outDir))){
  dir.create(file.path(outDir))
}

###

axts<-lapply(species2, function(x) {
  lel <- list.files(paste0(root_folders[species1], "/vs", species2),
             pattern = "net\\.axt", full.names = T)
  lel<-lel[!grepl("Exon", lel)]
  lel<-lel[!grepl("broken", lel)]
  tfn<-list.files(file.path(root_folders[species1],
                            "genome"),
                  pattern = "2bit",
                  full.names = T)
  qfn<-list.files(file.path(root_folders[species2],
                            "genome"),
                  pattern = "2bit",
                  full.names = T)
  out<-readAxt(lel, tAssemblyFn=tfn, qAssemblyFn=qfn)
})

names(axts)<-species2

print("Axts Read")

args <- list(length = 4)
args[[4]] <- T
chrom <- "1"

if(args[4] == "TRUE"){
  chr_axts<-axts[[1]][which(seqnames(first(axts[[1]])) == chrom)]
} else {
  chr_axts<-axts[[1]]
}
rm(axts)

axtSubIndex<-cut2(1:length(chr_axts), m=5000)
print(paste0("number of cuts: ", length(levels(axtSubIndex))))

chr_axts_split<-split(chr_axts, axtSubIndex)

rm(chr_axts)


for(i in 1:length(chr_axts_split)){
  seqs<-getLengthsOfIdenticalSeqs(chr_axts_split[[i]])
  for(j in 1:length(seqs)){
    if(j==1){
      write.table(seqs[[j]], paste0(outDir, species1, "_", species2, "_", chrom, "_", i, "_identical_seqs.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
    } else {
      write.table(seqs[[j]], paste0(outDir, species1, "_", species2, "_", chrom, "_", i, "_identical_seqs.tsv"), sep="\t", col.names=F, row.names=F, append=T, quote=F)
    }
  }
  rm(seqs)
  gc()
  print(i)
}

print("Script Complete")