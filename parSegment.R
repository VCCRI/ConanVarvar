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


suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(magrittr))
suppressMessages(library(fastseg))
suppressMessages(library(dplyr))

# This function performs segmentation on copy-number data.
parSegment <- function(copy.number, bin.size, ncores = 4) {

  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  segmentation.results = foreach(i = 1:length(copy.number@elementMetadata), .combine = rbind, .packages = c('magrittr', 'fastseg', 'dplyr')) %dopar% {

    avail.chromosomes <- copy.number[,i] %>% as.data.frame() %>%
      select(-start, -end, -width, -strand) %>% `colnames<-`(c('seqnames','CN')) %>%
      group_by(seqnames) %>% summarise(na = all(is.na(CN))) %>% filter(na == FALSE) %>%
      select(seqnames) %>% unlist()

    suppressMessages(fastseg(
      copy.number[,i][copy.number[,i]@seqnames %in% avail.chromosomes],
      alpha = 0.01,
      cyberWeight = 1,
      minSeg = as.integer(1e6 / bin.size))) %>%
      as.data.frame()
  }

  stopCluster(cl)
  return(segmentation.results %>% GRanges() %>% sort() %>% as.data.frame())
}
