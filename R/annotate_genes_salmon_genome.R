library(rtracklayer)
library(AnnotationHub)
library(plyranges)

gene_model <- import("/mnt/orca/projects/AQUA-FAANG/species/Salmo_salar/GCA_905237065.2/geneset/2021_07/Salmo_salar-GCA_905237065.2-2021_07-genes.gtf.gz")

ah <- AnnotationHub()
ssal_ensembl_105 <- ah[["AH98170"]]

salmon_genes <- gene_model %>% filter(type == "gene")
ens_id <- salmon_genes$gene_id

ens_id_to_name <- AnnotationDbi::select(ssal_ensembl_105, 
                                        keytype = "GENEID",
                                        keys = ens_id, 
                                        columns = c("GENENAME", "SYMBOL"))

ens_id_to_name_vector <- ens_id_to_name$GENENAME
names(ens_id_to_name_vector) <- ens_id_to_name$GENEID

genes_with_names <- salmon_genes
mcols(genes_with_names) <- NULL
genes_with_names$name <- salmon_genes$gene_id
genes_with_names$name <- ens_id_to_name_vector[salmon_genes$gene_id]
not_in_release <- which(is.na(genes_with_names$name))
genes_with_names$name[not_in_release] <- salmon_genes$gene_id[not_in_release]
no_name <- which(genes_with_names$name == "")
genes_with_names$name[no_name] <- names(no_name)

export.gff3(genes_with_names, "/mnt/orca/projects/AQUA-FAANG/species/Salmo_salar/GCA_905237065.2/geneset/2021_07/Salmo_salar-GCA_905237065.2-2021_07-genes_common_names_DB.gff3")
