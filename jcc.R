library(jcc)
library(dplyr)
library(BSgenome)
library("BSgenome.Oindica.ASM465v1")

biasMod <- fitAlpineBiasModel(gtf = "IndicaRTD.gtf",bam = "../12_1_paired_2ndPass_Aligned.sortedByCoord.out.bam",organism = "Oryza sativa",genome = Oindica, genomeVersion = "v1", version = 1.1, minLength = 600, maxLength = 7000, minCount = 500, maxCount = 10000, subsample = TRUE, nbrSubsample = 200, seed = 1, minSize = NULL, maxSize = NULL, verbose = TRUE)

tx2gene <- readr::read_rds("tx2gene.sub.rds")

predCovProfiles <- predictTxCoverage(biasModel = biasMod$biasModel, exonsByTx = biasMod$exonsByTx, bam = "../12_1_paired_2ndPass_Aligned.sortedByCoord.out.bam",  tx2gene = tx2gene, genome = Oindica, genes = def, nCores = 30, verbose = TRUE)

txQuants <- readr::read_rds("12_1_quant.sub.rds")

txsc <- scaleTxCoverages(txCoverageProfiles = predCovProfiles, txQuants = txQuants, tx2gene = tx2gene, strandSpecific = TRUE, methodName = "Salmon", verbose = TRUE)

jcov <- read.delim("../12_1_paired_2ndPass_SJ.out.tab", header = FALSE, as.is = TRUE) %>% setNames(c("seqnames", "start", "end", "strand", "motif", "annot", "uniqreads", "mmreads", "maxoverhang")) %>% dplyr::mutate(strand = replace(strand, strand == 1, "+")) %>% dplyr::mutate(strand = replace(strand, strand == 2, "-")) %>% dplyr::select(seqnames, start, end, strand, uniqreads, mmreads) %>% dplyr::mutate(seqnames = as.character(seqnames))

combCov <- combineCoverages(junctionCounts = jcov, junctionPredCovs = txsc$junctionPredCovs, txQuants = txsc$txQuants)jcc <- calculateJCCScores(junctionCovs = combCov$junctionCovs, geneQuants = combCov$geneQuants)

write.table((jcc$geneScores),file="12_1_jccscore_test.txt", row.names = FALSE)

save.image(file='12_1_jcc_test.Rdata‘)
