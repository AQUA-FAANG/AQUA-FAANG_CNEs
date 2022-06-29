library(rtracklayer)

sumitFiles <- list.files("/mnt/orca/projects/AQUA-FAANG/ATAC/Salmon/ATAC_Ss_fromDiego", 
                         pattern = "summits", full.names = T)

sampleNames <- sub("_summits.bed", "", basename(sumitFiles))
summits <- lapply(sumitFiles, import)

names(summits) <- sampleNames

lapply(sampleNames, function(x){
  export(summits[[x]] + 750, 
         file.path("/mnt/orca/projects/AQUA-FAANG/ATAC/Salmon/ATAC_Ss_fromDiego",
               paste0(x, "_750flank.bed")))
})
