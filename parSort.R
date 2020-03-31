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


base::suppressMessages(base::library(Rsamtools))
base::suppressMessages(base::library(foreach))
base::suppressMessages(base::library(doParallel))
base::suppressMessages(base::library(magrittr))

# This function sorts input BAM files using Samtools.
# Each file is analysed on a separate CPU.
parSort <- function(bam.dir, ncores = 4) {

  bam.files <- base::list.files(base::file.path(bam.dir), pattern = '.bam$', full.names = FALSE) %>%
    base::grep(pattern = 'SORTED', inv = TRUE, value = TRUE)

  sorted.bam.files <- base::list.files(base::file.path(bam.dir),
                                       pattern = '^SORTED.*.bam$',
                                       full.names = FALSE)

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  list <- foreach::foreach(i = 1:base::length(bam.files)) %dopar% {

    bam.file <- base::paste0('SORTED.', bam.files[i], '.bam')

    # ignore if the file has already been sorted
    if (bam.file %in% sorted.bam.files) {
      base::paste(bam.dir, bam.file, sep = '/')
    }
    else {
      Rsamtools::sortBam(base::paste(bam.dir, bam.files[i], sep = '/'),
                         base::paste0(bam.dir, '/SORTED.', bam.files[i]))
    }
  }

  parallel::stopCluster(cl)
  return(base::as.character(list))
}
