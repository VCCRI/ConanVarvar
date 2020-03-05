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


suppressMessages(library(exomeCopy))
suppressMessages(library(Rsamtools))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))

# This function counts how many reads fall into each of the pre-specified bin in the input BAM files.
# Each file is analysed on a separate CPU.
parCount <- function(bam.files, bins, min.mapq = 0.8,
                     skip.indexing.BAM = FALSE, ncores = 4) {

  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  counts <- foreach(i = 1:length(bam.files)) %dopar% {
    if (!skip.indexing.BAM) {Rsamtools::indexBam(bam.files[i])}
    exomeCopy::countBamInGRanges(bam.file = bam.files[i],
                                 granges = bins,
                                 min.mapq = min.mapq,
                                 get.width = TRUE)
  }

  stopCluster(cl)
  counts <- data.frame(counts)
  colnames(counts) = bam.files
  return(counts)
}
