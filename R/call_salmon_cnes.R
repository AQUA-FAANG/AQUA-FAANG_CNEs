#!/usr/bin/Rscript
#SBATCH -c 1
#SBATCH --mem 16G
#SBATCH -p daniocode
#SBATCH --array 1-6


library(CNEr)
library(rtracklayer)
library(plyranges)
library(dplyr)

idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# this will be a test run for Salmon - Zebrafish

root_folders <- list.files(path = "/mnt/orca/projects/AQUA-FAANG/species/", 
                           recursive = T, pattern = "^GCA_\\d+\\.\\d$",
                           include.dirs = T, full.names = T)
names(root_folders) <-  c("Carp", "Zebrafish", "Seabass", "Pike", "Trout", 
                          "Salmon", "Turbot", "Seabream")

vs_species <- c("Carp", "Seabass", "Pike", "Trout", "Turbot", "Seabream")
querySpecies <- vs_species[idx]

cutoffs <- c(16L, 8L, 8L, 16L, 8L, 8L)
names(cutoffs) <- vs_species

axtSalmonFile <- list.files(paste0(root_folders["Salmon"], "/vs", querySpecies),
                        pattern = "net\\.axt", full.names = T)
axtZebrafishFile <- list.files(file.path(root_folders[querySpecies], "vsSalmon"),
                           pattern = "net\\.axt", full.names = T)
salmonFastaFile <- list.files(file.path(root_folders["Salmon"],
                                        "genome"),
                              pattern = "softmasked.fa$",
                              full.names = T)
zebrafishFastaFile <- list.files(file.path(root_folders[querySpecies],
                                           "genome"),
                                 pattern = "softmasked.fa$",
                                 full.names = T)

# axtSalmon <- readAxt(axtSalmonFile, tAssemblyFn = salmonFastaFile,
#                      qAssemblyFn = zebrafishFastaFile)
# axtZebrafish <- readAxt(axtZebrafishFile, tAssemblyFn = zebrafishFastaFile,
#                         qAssemblyFn = salmonFastaFile)

salmonRepeats <- read.rmskFasta(salmonFastaFile)
zebrafishRepeats <- read.rmskFasta(zebrafishFastaFile)

salmonExons <- import.gff3(list.files(file.path(root_folders["Salmon"], 
                                                "geneset"), 
                                      pattern = "gff", 
                                      recursive = T, full.names = T)) %>%
  filter(type == "exon")

zebrafishExons <- import.gff3(list.files(file.path(root_folders[querySpecies], 
                                                "geneset"), 
                                      pattern = "gff", 
                                      recursive = T, full.names = T)) %>%
  filter(type == "exon")

salmonFilter <- c(salmonExons, salmonRepeats)
zebrafishFilter <- c(zebrafishExons, zebrafishRepeats)

salmon2bit <- list.files(file.path(root_folders["Salmon"],
                                        "genome"),
                              pattern = "2bit",
                              full.names = T)
zebrafish2bit <- list.files(file.path(root_folders[querySpecies],
                                           "genome"),
                                 pattern = "2bit",
                                 full.names = T)

cneSalSalVsDanRer <- CNE(assembly1Fn = salmon2bit,
                         assembly2Fn = zebrafish2bit, 
                         axt12Fn = axtSalmonFile,
                         axt21Fn = axtZebrafishFile,
                         cutoffs1 = 16L, cutoffs2 = cutoffs[querySpecies])

if (querySpecies == "Trout"){
  pecIdentity <- c(100L, 100L, 100L, 98L)
  window <- c(100L, 75L, 50L, 50L)
  identity <- as.integer(window * (pecIdentity / 100))
}else if (querySpecies == "Pike"){
  pecIdentity <- c(100L, 100L, 100L, 100L, 98L, 96L, 96L, 90L)
  window <- c(100L, 75L, 50L, 30L, 50L, 50L, 30L, 50L)
  identity <- as.integer(window * (pecIdentity / 100))
}else{
  pecIdentity <- rep(c(100L, 98L, 96L, 90L, 80L, 70L), each = 2)
  window <- rep(c(50L, 30L), times = 6)
  identity <- as.integer(window * (pecIdentity / 100))
}

cneSalSalVsDanRerList <- ceScan(x=cneSalSalVsDanRer, tFilter=salmonFilter,
                                qFilter=zebrafishFilter,
                                window=window, identity=identity)

cneSalSalVsDanRerList <- ceScan(x=cneSalSalVsDanRer, tFilter=salmonFilter,
                                qFilter=zebrafishFilter,
                                window=window, identity=identity)

cneMergedListSalSalVsDanRer <- lapply(cneSalSalVsDanRerList, cneMerge)

cneFinalListDSalSalVsDanRer <- lapply(cneMergedListSalSalVsDanRer, blatCNE)

saveRDS(cneFinalListDSalSalVsDanRer, 
        paste0("/mnt/orca/projects/AQUA-FAANG/species/CNEs/", 
               "Salmon_", querySpecies, ".RDS"))
