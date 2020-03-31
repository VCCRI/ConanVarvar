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


base::suppressMessages(base::library(foreach))
base::suppressMessages(base::library(doParallel))
base::suppressMessages(base::library(magrittr))
base::suppressMessages(base::library(fastseg))
base::suppressMessages(base::library(dplyr))

# This function performs segmentation on copy-number data.
parSegment <- function(copy.number, bin.size, ncores = 4) {

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  segmentation.results = foreach::foreach(i = 1:base::length(copy.number@elementMetadata), .combine = rbind, .packages = base::c('magrittr', 'fastseg', 'dplyr')) %dopar% {

    avail.chromosomes <- copy.number[,i] %>% base::as.data.frame() %>%
      dplyr::select(-start, -end, -width, -strand) %>% `colnames<-`(base::c('seqnames','CN')) %>%
      dplyr::group_by(seqnames) %>% dplyr::summarise(na = base::all(base::is.na(CN))) %>% dplyr::filter(na == FALSE) %>%
      dplyr::select(seqnames) %>% base::unlist()

    base::suppressMessages(fastseg::fastseg(
      copy.number[,i][copy.number[,i]@seqnames %in% avail.chromosomes],
      alpha = 0.01,
      cyberWeight = 1,
      minSeg = base::as.integer(1e6 / bin.size))) %>%
      base::as.data.frame()
  }

  parallel::stopCluster(cl)
  return(segmentation.results %>% GenomicRanges::GRanges() %>% base::sort() %>% base::as.data.frame())
}
