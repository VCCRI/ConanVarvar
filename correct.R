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


base::suppressMessages(base::library(doParallel))
base::suppressMessages(base::library(foreach))

base::source(base::paste0(working.dir, '/correctReadcount.R'))

# This function corrects raw read counts for GC content and mappability biases and convert them to scaled copy number.
correct <- function(bins, counts, somes,
                    min.mapq, rough.span, final.span) {

  num.files <- base::dim(counts)[2]
  samples <- base::basename(base::colnames(counts))
  copy.number <- base::matrix(, nrow = base::length(bins), ncol = 0)
  for (i in 1:num.files) {
    bins$reads <- counts[,i]
    CN.one.sample <- base::list()
    for (chrom in 1:base::length(somes)) {
      chrom.bins <- bins[bins@seqnames == somes[chrom]]
      if (base::all(chrom.bins$reads == 0)) {
        CN.one.sample <- base::append(CN.one.sample, base::rep(NA, base::length(chrom.bins)))
      } else {
        CN.chr <- correctReadcount(x = chrom.bins,
                                   mappability = min.mapq,
                                   rough.span = rough.span,
                                   final.span = final.span,
                                   verbose = FALSE)
        CN.one.sample <- base::append(CN.one.sample, CN.chr$copy)
      }
    }
    if (base::all(base::is.na(CN.one.sample))) {
      base::stop(base::paste0(
        "Found a sample with zero read counts in all chromosomes: ",
        base::colnames(counts)[i],
        ". Please make sure the reference (hg19/hg38) and format (NCBI/UCSC) are correct."
      ))
    }
    copy.number <- base::cbind(copy.number, base::unlist(CN.one.sample))
  }

  S4Vectors::mcols(bins) <- base::data.frame(copy.number)
  base::colnames(S4Vectors::mcols(bins)) <- samples
  return(bins)
}
