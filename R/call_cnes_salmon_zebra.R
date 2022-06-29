library(CNEr)
library(rtracklayer)
library(plyranges)
library(dplyr)

# this will be a test run for Salmon - Zebrafish

root_folders <- list.files(path = "/mnt/orca/projects/AQUA-FAANG/species/", 
                           recursive = T, pattern = "^GCA_\\d+\\.\\d$",
                           include.dirs = T, full.names = T)
names(root_folders) <-  c("Carp", "Zebrafish", "Seabass", "Pike", "Trout", 
                          "Salmon", "Turbot", "Seabream")

axtSalmonFile <- list.files(file.path(root_folders["Salmon"], "vsZebrafish"),
                        pattern = "net\\.axt", full.names = T)
axtZebrafishFile <- list.files(file.path(root_folders["Zebrafish"], "vsSalmon"),
                           pattern = "net\\.axt", full.names = T)
salmonFastaFile <- list.files(file.path(root_folders["Salmon"],
                                        "genome"),
                              pattern = "softmasked.fa$",
                              full.names = T)
zebrafishFastaFile <- list.files(file.path(root_folders["Zebrafish"],
                                           "genome"),
                                 pattern = "softmasked.fa$",
                                 full.names = T)

axtSalmon <- readAxt(axtSalmonFile, tAssemblyFn = salmonFastaFile,
                     qAssemblyFn = zebrafishFastaFile)
axtZebrafish <- readAxt(axtZebrafishFile, tAssemblyFn = zebrafishFastaFile,
                        qAssemblyFn = salmonFastaFile)

salmonRepeats <- read.rmskFasta(salmonFastaFile)
zebrafishRepeats <- read.rmskFasta(zebrafishFastaFile)

salmonExons <- import.gff3(list.files(file.path(root_folders["Salmon"], 
                                                "geneset"), 
                                      pattern = "gff", 
                                      recursive = T, full.names = T)) %>%
  filter(type == "exon")

zebrafishExons <- import.gff3(list.files(file.path(root_folders["Zebrafish"], 
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
zebrafish2bit <- list.files(file.path(root_folders["Zebrafish"],
                                           "genome"),
                                 pattern = "2bit",
                                 full.names = T)

cneSalSalVsDanRer <- CNE(assembly1Fn = salmon2bit,
                         assembly2Fn = zebrafish2bit, 
                         axt12Fn = axtSalmonFile,
                         axt21Fn = axtZebrafishFile,
                         cutoffs1 = 16L, cutoffs2 = 8L)

pecIdentity <- rep(c(100L, 98L, 96L, 90L, 80L, 70L), each = 2)
# identity need to be integer real values, not percentages
identity <- as.integer(window * (pecIdentity / 100))
window <- rep(c(50L, 30L), times = 6)

cneSalSalVsDanRerList <- ceScan(x=cneSalSalVsDanRer, tFilter=salmonFilter,
                                qFilter=zebrafishFilter,
                                window=window, identity=identity)

cneMergedListSalSalVsDanRer <- lapply(cneSalSalVsDanRerList, cneMerge)

cneFinalListDSalSalVsDanRer <- lapply(cneMergedListSalSalVsDanRer, blatCNE)
