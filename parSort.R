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


suppressMessages(library(Rsamtools))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(magrittr))

# This function sorts input BAM files using Samtools.
# Each file is analysed on a separate CPU.
parSort <- function(bam.dir, ncores = 4) {

  bam.files <- list.files(file.path(bam.dir), pattern = '.bam$', full.names = FALSE) %>%
    grep(pattern = 'SORTED', inv = TRUE, value = TRUE)

  sorted.bam.files <- list.files(file.path(bam.dir),
                                 pattern = '^SORTED.*.bam$',
                                 full.names = FALSE)

  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  list <- foreach(i = 1:length(bam.files)) %dopar% {

    bam.file <- paste0('SORTED.', bam.files[i], '.bam')

    # ignore if the file has already been sorted
    if (bam.file %in% sorted.bam.files) {
      paste(bam.dir, bam.file, sep = '/')
    }
    else {
      Rsamtools::sortBam(paste(bam.dir, bam.files[i], sep = '/'),
                         paste0(bam.dir, '/SORTED.', bam.files[i]))
    }
  }

  stopCluster(cl)
  return(as.character(list))
}
