library(GenomicAlignments)
library(rtracklayer)
library(genomation)
library(tidyverse)

peakSummits <- import(
  "/mnt/orca/projects/AQUA-FAANG/ATAC/Salmon/4_Mid_somitogen.mRp.clN_summits.bed")
canonical_idx <- which(seqnames(peakSummits) %in% as.character(1:30))
peakSummits <- peakSummits[canonical_idx]

narrowPeaks <- import("/mnt/orca/projects/AQUA-FAANG/ATAC/Salmon/4_Mid_somitogen.mRp.clN_peaks.narrowPeak")
suspicious_peaks <- import.bed("/mnt/orca/projects/AQUA-FAANG/ChIP/Salmon/suspicious_peaks.narrowPeak")
filter_peak <- subsetByOverlaps(peakSummits, suspicious_peaks, invert = T)
  
peakSummits750bpFlank <- filter_peak + 750

atacBam <- "/mnt/orca/projects/AQUA-FAANG/ATAC/Salmon/4_Mid_somitogen.mRp.clN.sorted.bam"

param1 <- ScanBamParam(which = peakSummits750bpFlank)
atacReads <- readGAlignments(atacBam, param = param1) 

atacReadsGR <- granges(atacReads)
atacTn5Sites <- resize(atacReadsGR, width = 1)
tn5Shift <- ifelse(strand(atacTn5Sites) == "+", 4, -5)
atacTn5Corrected <- shift(atacTn5Sites, tn5Shift)

atacBins <- ScoreMatrixBin(target = atacTn5Corrected,
                           windows = peakSummits750bpFlank,
                           bin.num = 13, bin.op = "sum")

nucleoAtacSignal <- import("/mnt/orca/projects/AQUA-FAANG/ATAC/Salmon/nucleoatac/4_Mid_somitogen.nucleoatac_signal.smooth.bw")

nucleoAtacBins <- ScoreMatrixBin(target = nucleoAtacSignal,
                           windows = peakSummits750bpFlank,
                           bin.num = 13, bin.op = "mean", weight.col = "score")

chipFiles <- list.files("/mnt/orca/projects/AQUA-FAANG/ChIP/Salmon/", 
                        pattern = "Ss4_K.+\\.bam$", full.names = T)

chipFiles <- chipFiles[-2]
chipReads <- lapply(chipFiles, function(x){
  x %>% readGAlignmentPairs() %>% granges() %>%
    IRanges::resize(width = 1, fix = "center")
})

chipBins <- lapply( chipReads, ScoreMatrixBin,
                            windows = peakSummits750bpFlank,
                            bin.num = 13, bin.op = "sum")

binList <- c(atacBins, chipBins[[1]], chipBins[[2]], 
             chipBins[[3]], nucleoAtacBins)
binList <- intersectScoreMatrixList(binList)

scaledMatrices <- lapply(binList, function(x){
    matrix(
      base::scale(
        as.vector(x@.Data),
        scale = T, center = F),
      ncol = 13)
})

umap_input <- purrr::reduce(scaledMatrices, cbind)
salmon_umap_1 <- uwot::umap(umap_input)

peaks <- filter_peak[as.numeric(row.names(binList[[1]]))]

saveRDS(list(peaks = peaks, umap = salmon_umap_1), 
        "/mnt/orca/projects/AQUA-FAANG/ChIP/Salmon/salmon_ss4_umap_v3.RDS")
