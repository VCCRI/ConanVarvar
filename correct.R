# This file is part of ConanVarvar, a tool for detection of copy number variants
#
# Copyright (C) 2020 Victor Chang Cardiac Research Institute
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is provided without warranty of any kind.
#
# <http://www.gnu.org/licenses/>


suppressMessages(library(doParallel))
suppressMessages(library(foreach))

source(paste0(working.dir, '/correctReadcount.R'))

# This function corrects raw read counts for GC content and mappability biases and convert them to scaled copy number.
correct <- function(bins, counts, somes,
                    min.mapq, rough.span, final.span) {

  num.files <- dim(counts)[2]
  samples <- basename(colnames(counts))
  copy.number <- matrix(, nrow = length(bins), ncol = 0)
  for (i in 1:num.files) {
    bins$reads <- counts[,i]
    CN.one.sample <- list()
    for (chrom in 1:length(somes)) {
      chrom.bins <- bins[bins@seqnames == somes[chrom]]
      if (all(chrom.bins$reads == 0)) {
        CN.one.sample <- append(CN.one.sample, rep(NA, length(chrom.bins)))
      } else {
        CN.chr <- correctReadcount(x = chrom.bins,
                                   mappability = min.mapq,
                                   rough.span = rough.span,
                                   final.span = final.span,
                                   verbose = FALSE)
        CN.one.sample <- append(CN.one.sample, CN.chr$copy)
      }
    }
    if (all(is.na(CN.one.sample))) {
      stop(paste0(
        "Found a sample with zero read counts in all chromosomes: ",
        colnames(counts)[i],
        ". Please make sure the reference (hg19/hg38) and format (NCBI/UCSC) are correct."
      ))
    }
    copy.number <- cbind(copy.number, unlist(CN.one.sample))
  }

  mcols(bins) <- data.frame(copy.number)
  colnames(mcols(bins)) <- samples
  return(bins)
}
