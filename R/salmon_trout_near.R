#!/usr/bin/Rscript
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -p daniocode

load("data/alignment_env.RData")
library(CNEr)
library(tidyverse)
# if (!file.exists(paste0(root_folders[subject], "/lastdb"))){
#   dir.create(paste0(root_folders[subject], "/lastdb"))
#   system2(command = "lastdb", args = c("-c", 
#           paste0(root_folders[subject], "/lastdb/",subject), 
#           list.files(root_folders[subject], pattern = "softmasked\\.fa$", 
#                      full.names = T, recursive = T)
#           ))
# }

#idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
subject <- "Salmon"
query <- "Trout"

aln_folder <- paste0(root_folders[subject], "/vs", query)
dir.create(aln_folder)
maf <- paste0(aln_folder, "/", 
              paste(subject, query, "maf", sep = "."))
lastal(db=paste0(root_folders[subject], "/lastdb/",subject),
       queryFn=list.files(root_folders[query], pattern = "softmasked\\.fa$", 
                          full.names = T, recursive = T),
       outputFn=maf,
       distance="near", binary="lastal", mc.cores=8L)

psls <- paste0(aln_folder, "/",  paste(subject, query, "psl", sep = "."))
system2(command="maf-convert", args=c("psl", 
                                      maf,
                                      ">", psls))

## Join close alignments
assemblyTarget <- list.files(root_folders[subject], pattern = "\\.2bit", 
                             full.names = T, recursive = T)
assemblyQuery <- list.files(root_folders[query], pattern = "\\.2bit", 
                            full.names = T, recursive = T)
chains <- axtChain(psls, assemblyTarget=assemblyTarget,
                   assemblyQuery=assemblyQuery, 
                   distance="near",
                   removePsl=FALSE, binary="axtChain")

## Sort and combine
allChain <- chainMergeSort(chains, assemblyTarget, assemblyQuery,
                           allChain= paste0(aln_folder, "/",  
                                            paste(subject, query,
                                                  ".all.chain", sep = ".")), 
                           removeChains=FALSE, binary="chainMergeSort")

## Filtering out chains
allPreChain <- chainPreNet(allChain, assemblyTarget, assemblyQuery,
                           allPreChain=paste0(aln_folder,  "/", 
                                              paste(subject, query,
                                                    ".all.pre.chain", 
                                                    sep = ".")),
                           removeAllChain=FALSE, binary="chainPreNet")

## Keep the best chain and add synteny information
netSyntenicFile <- chainNetSyntenic(allPreChain, assemblyTarget, assemblyQuery,
                                    netSyntenicFile=paste0(aln_folder, "/",  
                                                           paste(subject, query,
                                                                 ".noClass.net", sep = ".")),
                                    binaryChainNet="chainNet", 
                                    binaryNetSyntenic="netSyntenic")

netToAxt(netSyntenicFile, allPreChain, assemblyTarget, assemblyQuery,
         axtFile=paste0(aln_folder, "/",  
                        paste(subject, query, 
                              ".net.axt", sep = ".")),
         removeFiles=FALSE,
         binaryNetToAxt="netToAxt", binaryAxtSort="axtSort")