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
suppressMessages(library(GenomicRanges))
suppressMessages(library(IRanges))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))

# This function extracts copy number values for a specified sample, chromosome and positions on the chromosome.
extract.copy.number <- function(copy.number.data, ID, seq, extract.from, extract.to) {

  segment.id <- if (grepl("^[0-9]", ID)) paste0('X', ID) else ID

  segment.CN <- select(copy.number.data, all_of(segment.id), seqnames, start, end) %>%
    filter(start >= extract.from,
           end <= extract.to,
           as.character(seqnames) == as.character(seq)) %>%
    select(all_of(segment.id)) %>% unlist()
}

# This function analyses the obtained copy-number segments, calculates per-variant statistics and summarises the results.
postSegment <- function(copy.number, segmentation.results, del.threshold, dup.threshold,
                        seg.duplications, seg.duplications.threshold, syndromes, working.dir, output.dir,
                        verbose, plot.results, ncores = 4) {

  # Calculate real mean and standard deviation for each segment
  segments.mu.sigma <- sapply(seq_along(segmentation.results$seg.mean), function(i) {
    segment.CN <-extract.copy.number(copy.number.data = as.data.frame(copy.number),
                                     ID = segmentation.results$ID[[i]],
                                     seq = segmentation.results$seqnames[[i]],
                                     extract.from = segmentation.results$start[[i]],
                                     extract.to = segmentation.results$end[[i]])
    c(mean(segment.CN, na.rm = TRUE), sd(segment.CN, na.rm = TRUE))
  })
  segmentation.results$mu <-t(segments.mu.sigma)[, 1]
  segmentation.results$sigma <-t(segments.mu.sigma)[, 2]

  # Remove segments corresponding to "blank spaces"
  segmentation.results <- segmentation.results %>% dplyr::filter(!is.na(mu) & !is.na(sigma))

  number.of.segments <- dim(segmentation.results)[1]

  if (number.of.segments < 30) {

    # Use thresholds if there are too few segments
    if (verbose) message("Too few segments. Using thresholds to identify deletions and duplications.")
    segmentation.results$class <- 'NORMAL'
    segmentation.results$class[segmentation.results$seg.mean <= del.threshold] <- 'DELETION'
    segmentation.results$class[segmentation.results$seg.mean >= dup.threshold] <- 'DUPLICATION'
    segmentation.results$class <- factor(segmentation.results$class)

    potential.CNVs <- segmentation.results %>% filter(!(class %in% c('NORMAL')))
    potential.CNVs$class <- droplevels(potential.CNVs$class)

    # Retain only deletions and duplications
    potential.CNVs <- potential.CNVs %>% mutate(original.mu = mu)

  } else {

    # Remove sigma outliers
    segmentation.results <- segmentation.results[-which(segmentation.results$sigma %in% (segmentation.results$sigma %>% boxplot(plot = FALSE) %$% out)),]

    if (all(segmentation.results$mu > del.threshold & segmentation.results$mu < dup.threshold)){

      # If there are no abnormalities, there is nothing to report
      if (verbose) message("No CNVs found!")
      return(invisible(NULL))

    } else {

      # Transform per-segment mean values so that all potential CNVs are centered around zero
      max.mu <- max(abs(segmentation.results$mu), na.rm = TRUE)
      segmentation.results <- segmentation.results %>% mutate(original.mu = mu) %>% mutate(mu = mu/max.mu) %>% mutate(mu = sign(mu)*log(abs(mu)))

      # Cluster the transformed mean values
      kmeans.clusters <- segmentation.results$mu %>% kmeans(centers = c(min(segmentation.results$mu), 0, max(segmentation.results$mu)), nstart = 10)

      # Retain only deletions and duplications
      CNV.labels <- kmeans.clusters$centers %>% as.data.frame() %>% rownames_to_column() %>% setNames(c('cluster', 'centers')) %>% arrange(centers) %>% mutate(class = c('Normal+', 'CNV', 'Normal-')) %>% arrange(cluster) %>% select(class) %>% unlist()
      segmentation.results$cluster <- kmeans.clusters %$% cluster %>% factor(levels = c(1, 2, 3), labels = CNV.labels)
      potential.CNVs <- segmentation.results %>% filter(cluster %in% c('CNV'))
      potential.CNVs$cluster <- droplevels(potential.CNVs$cluster)

      # Separate deletions from duplications based on the original mean
      potential.CNVs <- potential.CNVs %>% rowwise() %>% mutate(class = if (original.mu < 0) "DELETION" else "DUPLICATION")

    }
  }

  if (nrow(potential.CNVs) == 0) {

    if (verbose) message("No CNVs found!")
    return(invisible(NULL))

  } else {

    # For each observed segment length, generate a null distribution using bootstrap
    bins <- copy.number %>% as.data.frame() %>% select(-seqnames, -start, -end, -width, -strand) %>% sapply(function(i) i) %>% as.vector()
    all.lengths <- sort(unique(potential.CNVs$num.mark))
    cl <- makeCluster(4)
    registerDoParallel(cl)
    mu.null.distribution <- foreach(segment.length = all.lengths, .packages = c('magrittr', 'dplyr')) %dopar% {
      number.of.segments <- potential.CNVs %>%
        filter(num.mark == segment.length) %>%
        nrow()
      real.segment.statistics <- potential.CNVs %>%
        filter(num.mark == segment.length) %>%
        select(original.mu) %>%
        unlist()
      fake.segment.statistics <- replicate(1000 - number.of.segments, {
        repeat {
          fake.segment <- sample(bins, size = segment.length, replace = TRUE)
          if (sum(is.na(fake.segment) == FALSE) > 1) break
        }
        mean(fake.segment, na.rm = TRUE)
      })
      c(fake.segment.statistics, real.segment.statistics) %>% sort()
    }
    stopCluster(cl)
    names(mu.null.distribution) <- paste0('X', all.lengths)

    # Calculate one-sided p-values for the identified potential CNVs
    p.value.deletions <- function(distribution, observation) {
      return(mean(round(observation, 6) >= round(distribution, 6)))
    }
    p.value.duplications <- function(distribution, observation) {
      return(mean(round(observation, 6) <= round(distribution, 6)))
    }
    potential.CNVs <- potential.CNVs %>% rowwise() %>%
      mutate(pvalue = case_when(
        class == "DELETION" ~ p.value.deletions(mu.null.distribution[[paste0('X', num.mark)]], original.mu),
        class == "DUPLICATION" ~ p.value.duplications(mu.null.distribution[[paste0('X', num.mark)]], original.mu))) %>%
      as.data.frame()

  }

  if (verbose) message('Generating the output...')

  if (plot.results) {

    # Generate per-chromosome plots
    suppressWarnings(dir.create(paste0(output.dir, '/Plots')))
    segmentation.results <- segmentation.results %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    potential.CNVs <- potential.CNVs %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    foreach(chr = unique(segmentation.results@seqnames) %>% as.vector()) %dopar% {
      source(paste0(working.dir, '/plotChromSegments.R'))
      pdf(paste0(output.dir, '/Plots/Chromosome', chr, '.pdf'), width = 12, height = 9)
      for (sample.id in unique(segmentation.results$ID)) {
        plotChromSegments(sample.id,
                          CN.values = copy.number[copy.number@seqnames == chr]@elementMetadata[sample.id] %>% unlist() %>% as.vector(),
                          CN.positions = copy.number[copy.number@seqnames == chr]@ranges@start,
                          segmentation.results = segmentation.results[segmentation.results@seqnames == chr & segmentation.results$ID == sample.id],
                          potential.CNVs = potential.CNVs[potential.CNVs@seqnames == chr & potential.CNVs$ID == sample.id], chr = chr,
                          seg.duplications = seg.duplications, seg.duplications.threshold = seg.duplications.threshold)
      }
      dev.off()
    }
    stopCluster(cl)

  }

  # Flag CNVs that overlap with large segmental duplications
  potential.CNVs <- potential.CNVs %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  overlapping <- subsetByOverlaps(potential.CNVs, seg.duplications[seg.duplications@ranges@width > seg.duplications.threshold], ignore.strand = TRUE)
  elementMetadata(potential.CNVs)[['overlaps.seg.duplication']] <- factor(potential.CNVs %in% overlapping,
                                                                          levels = c(FALSE, TRUE), labels = c('No','Yes'))

  # For each potential CNV, list any syndromes that could be caused by that CNV
  elementMetadata(potential.CNVs)[['syndrome']] <- NA
  overlapping <- suppressWarnings(findOverlaps(syndromes, potential.CNVs, ignore.strand = TRUE)) %>%
    as.data.frame() %>% group_by(subjectHits) %>% summarise(syndrome = list(queryHits))
  potential.CNVs[overlapping$subjectHits]$syndrome <- overlapping$syndrome
  potential.CNVs <- as.data.frame(potential.CNVs) %>%
    mutate(syndromic = !is.na(syndrome)) %>%
    rowwise %>%
    mutate(syndrome = ifelse(syndromic, list(as.data.frame(syndromes[syndrome]) %>% filter(as.character(type) == as.character(class)) %>% select(syndrome)), NA)) %>%
    mutate(syndrome = paste(unlist(syndrome), collapse = ', '))
  potential.CNVs$syndrome[potential.CNVs$syndrome %in% c('character(0)', 'NA')] <- ''

  # Pre-sort all potential CNVs according to their significance
  potential.CNVs <- as.data.frame(potential.CNVs) %>%
    group_by(seqnames, start, end, class) %>%
    mutate(freq = round(n()/length(unique(potential.CNVs$ID)), 2)) %>%
    ungroup() %>%
    arrange(pvalue, desc(syndrome)) %>%
    select(seqnames, start, end, width, ID, seg.mean, class, overlaps.seg.duplication, freq, syndrome, pvalue)

  # Final output
  colnames(potential.CNVs) <- c('Chromosome',
                                'Start', 'End',
                                'Width',
                                'ID',
                                'Copy Number log2 Mean Ratio',
                                'CNV type',
                                paste0('Overlaps with a segmental duplication >', seg.duplications.threshold / 1e3, 'kb'),
                                'Inter-batch frequency',
                                'Associated syndromes',
                                'P-value')
  suppressMessages(write.xlsx(potential.CNVs, paste0(output.dir, '/Variations.xlsx')))
}
