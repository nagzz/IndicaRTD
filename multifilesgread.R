#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(gread))
args = commandArgs(trailingOnly=TRUE)
gtf_file <- file.path(args[1])
gtf<-read_format(gtf_file)
ans <- construct_introns(gtf, update=FALSE)[]
write.table(ans,file="",row.names = FALSE, col.names = FALSE)
