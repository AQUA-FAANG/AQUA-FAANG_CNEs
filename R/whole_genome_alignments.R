#!/usr/bin/Rscript

# genome alignments for the AQUA-FAANG species

library(CNEr)
library(tidyverse)

species <- c("Salmon", "Trout", "Carp", "Turbot", "Seabass", "Seabream", "Zebrafish")

root_folders <- list.files(path = "/mnt/orca/projects/AQUA-FAANG/species/", 
                          recursive = T, pattern = "^GCA_\\d+\\.\\d$",
                          include.dirs = T, full.names = T)
names(root_folders) <- species[c(3, 7, 5, 2, 1, 4, 6)]

species_pairs <- combn(species, 2) %>% t() %>% as_tibble() %>% 
  dplyr::rename(Subject = V1, Query = V2) %>% 
  { rbind(., dplyr::rename(., Subject = Query, Query = Subject)) } %>% 
  dplyr::arrange(Subject, Query)

distance_matrix <- matrix(data = "far", nrow = length(species),
                          ncol = length(species))
rownames(distance_matrix) <- species
colnames(distance_matrix) <- species

medium_dist <- list(c(1, 2), c(2, 1), c(4, 5), c(4, 6), 
                    c(5, 4), c(5, 6), c(6, 4), c(6, 5),
                    c(3, 7), c(7, 3))

distance_matrix <- purrr::reduce(medium_dist, function(x, y) {
  x[y[1], y[2]] <- "medium"
  x
}, .init = distance_matrix)

library(BiocParallel)

params <- BatchtoolsParam(workers = 16, cluster = "slurm", 
                          resources = list(ncpus = 8, walltime = 4320,
                                           memory = 16384, partition = "low"))

bplapply(species, function(subject){
  load("data/alignment_env.RData")
  library(CNEr)
  library(tidyverse)
  dir.create(paste0(root_folders[subject], "/lastdb"))
  system2(command = "lastdb", args = c("-c", "-R10",
                                       paste0(root_folders[subject], 
                                              "/lastdb/",subject), 
                                       list.files(root_folders[subject], 
                                                  pattern = "softmasked\\.fa$", 
                                                  full.names = T, recursive = T)
  ))
}, BPPARAM = params)

bpmapply(function(subject, query){
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
  aln_folder <- paste0(root_folders[subject], "/vs", query)
  dir.create(aln_folder)
  maf <- paste0(aln_folder, "/", 
                paste(subject, query, "maf", sep = "."))
  lastal(db=paste0(root_folders[subject], "/lastdb/",subject),
         queryFn=list.files(root_folders[query], pattern = "softmasked\\.fa$", 
                            full.names = T, recursive = T),
         outputFn=maf,
         distance=distance_matrix[subject, query], binary="lastal", mc.cores=8L)
  
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
                     distance=distance_matrix[subject, query],
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
},
species_pairs$Subject, species_pairs$Query, BPPARAM = params
)

mapply(function(subject, query){
  vs_dir <- paste0(root_folders[subject], "/vs", query)
  dir.create(vs_dir)
  system2("ln", args = c("-s", list.files(paste0(root_folders[query], "/genome"),
                                         pattern = "softmasked\\.fa$",
                                         full.names = T, recursive = T),
                         paste0(vs_dir, "/")))
},
species_pairs$Subject, species_pairs$Query
)

# tsv ln create

mapply(function(subject, query){

  query_chromSizes <- list.files(file.path(root_folders[query], "genome"),
                                 pattern = "tsv$", full.names = T)
  system2("ln", args = c("-s", query_chromSizes,
                         paste0(root_folders[subject], "/vs", query)))

},
species_pairs$Subject, species_pairs$Query
)