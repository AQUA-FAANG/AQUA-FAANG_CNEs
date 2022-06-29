library(Biostrings)
library(GenomicAlignments)

nanopore_reads <- readGAlignments("/mnt/orca/damir/downloads/Ssalar_Embryo_Bam_Files_OEve/All_BC_primary_alignments.sorted.bam")
salmon_genome <- readDNAStringSet("/mnt/orca/projects/AQUA-FAANG/species/Salmo_salar/GCA_905237065.2/genome/Salmo_salar-GCA_905237065.2-softmasked.fa")
names(salmon_genome) <- unlist(lapply(names(salmon_genome), function(x) str_split(x, pattern = " ")[[1]][1]))

nanopore_cage_sym <- nanopore_reads %>% granges() %>%
  resize(width = 1, fix = "start")
nanopore_initiators <- nanopore_cage_sym %>% 
  { shift(., shift = ifelse(strand(.) == "+", -1, 1)) } %>% 
  resize(width = 2) %>% trim()
nanopore_dinucl <- getSeq(salmon_genome, nanopore_initiators)
nanopore_dinucl %>% as.data.frame() %>% group_by(x) %>% 
  summarise(freq = n()) %>%  mutate(x = fct_reorder(x, freq)) %>%
  ggplot(aes(x, freq)) + geom_col() + coord_flip()

gene_model <- import("/mnt/orca/projects/AQUA-FAANG/species/Salmo_salar/GCA_905237065.2/geneset/2021_07/Salmo_salar-GCA_905237065.2-2021_07-genes.gtf.gz")
nanopore_cage_sym_prom <- subsetByOverlaps(nanopore_cage_sym, 
                                        promoters(gene_model, 
                                                  upstream = 500,
                                                  downstream = 100))
promoter_initiators <- nanopore_cage_sym_prom %>% 
  { shift(., shift = ifelse(strand(.) == "+", -1, 1)) } %>% 
  resize(width = 2) %>% trim()
nanopore_dinucl_prom <- getSeq(salmon_genome, promoter_initiators)
nanopore_dinucl_prom %>% as.data.frame() %>% group_by(x) %>% 
  summarise(freq = n()) %>%  mutate(x = fct_reorder(x, freq)) %>% 
  ggplot(aes(x, freq)) + geom_col() + coord_flip()

nanopore_ctss <- unique(nanopore_cage_sym)
nanopore_ctss$score <- countOverlaps(nanopore_ctss, nanopore_cage_sym)
nanopore_ctss_split <- split(nanopore_ctss, strand(nanopore_ctss))
score(nanopore_ctss_split$`-`) <- score(nanopore_ctss_split$`-`) * -1
export.bw(nanopore_ctss_split$`+`, "/mnt/orca/damir/downloads/Ssalar_Embryo_Bam_Files_OEve/All_cage_sym_plus.bw")
export.bw(nanopore_ctss_split$`-`, "/mnt/orca/damir/downloads/Ssalar_Embryo_Bam_Files_OEve/All_cage_sym_minus.bw")
