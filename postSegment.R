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
base::suppressMessages(base::library(GenomicRanges))
base::suppressMessages(base::library(IRanges))
base::suppressMessages(base::library(tibble))
base::suppressMessages(base::library(magrittr))
base::suppressMessages(base::library(ggplot2))
base::suppressMessages(base::library(openxlsx))
base::suppressMessages(base::library(dplyr))
base::suppressMessages(base::library(igraph))

# This function extracts copy number values for a specified sample, chromosome and positions on the chromosome.
extract.copy.number <- function(copy.number.data, ID, seq, extract.from, extract.to) {

  segment.id <- if (base::grepl("^[0-9]", ID)) base::paste0('X', ID) else ID

  segment.CN <- dplyr::select(copy.number.data, dplyr::one_of(segment.id), seqnames, start, end) %>%
    dplyr::filter(start >= extract.from,
                  end <= extract.to,
                  base::as.character(seqnames) == base::as.character(seq)) %>%
    dplyr::select(dplyr::one_of(segment.id)) %>% base::unlist()
}

# This function analyses the obtained copy-number segments, calculates per-variant statistics and summarises the results.
postSegment <- function(copy.number, segmentation.results, del.threshold, dup.threshold,
                        seg.duplications, seg.duplications.threshold, syndromes, working.dir, output.dir,
                        occurrence.margin, bin.size,
                        verbose, plot.results, ncores = 4) {

  # Calculate real mean and standard deviation for each segment
  segments.mu.sigma <- base::sapply(base::seq_along(segmentation.results$seg.mean), function(i) {
    segment.CN <- extract.copy.number(copy.number.data = base::as.data.frame(copy.number),
                                      ID = segmentation.results$ID[[i]],
                                      seq = segmentation.results$seqnames[[i]],
                                      extract.from = segmentation.results$start[[i]],
                                      extract.to = segmentation.results$end[[i]])
    base::c(base::mean(segment.CN, na.rm = TRUE), stats::sd(segment.CN, na.rm = TRUE))
  })
  segmentation.results$mu <- base::t(segments.mu.sigma)[, 1]
  segmentation.results$sigma <- base::t(segments.mu.sigma)[, 2]

  # Remove segments corresponding to "blank spaces"
  segmentation.results <- segmentation.results %>% dplyr::filter(!base::is.na(mu) & !base::is.na(sigma))

  number.of.segments <- base::dim(segmentation.results)[1]

  if (number.of.segments < 30) {

    # Use thresholds if there are too few segments
    if (verbose) base::message("Too few segments. Using thresholds to identify deletions and duplications.")
    segmentation.results$class <- 'NORMAL'
    segmentation.results$class[segmentation.results$seg.mean <= del.threshold] <- 'DELETION'
    segmentation.results$class[segmentation.results$seg.mean >= dup.threshold] <- 'DUPLICATION'
    segmentation.results$class <- base::factor(segmentation.results$class)

    potential.CNVs <- segmentation.results %>% dplyr::filter(!(class %in% base::c('NORMAL')))
    potential.CNVs$class <- base::droplevels(potential.CNVs$class)

    # Retain only deletions and duplications
    potential.CNVs <- potential.CNVs %>% dplyr::mutate(original.mu = mu)

  } else {

    # Remove sigma outliers
    segmentation.results <- segmentation.results[-base::which(segmentation.results$sigma %in% (segmentation.results$sigma %>% graphics::boxplot(plot = FALSE) %$% out)),]

    if (base::all(segmentation.results$mu > del.threshold & segmentation.results$mu < dup.threshold)) {

      # If there are no abnormalities, there is nothing to report
      if (verbose) base::message("No CNVs found!")
      return(base::invisible(NULL))

    } else {

      # Transform per-segment mean values so that all potential CNVs are centered around zero
      max.mu <- base::max(base::abs(segmentation.results$mu), na.rm = TRUE)
      segmentation.results <- segmentation.results %>%
        dplyr::mutate(original.mu = mu) %>%
        dplyr::mutate(mu = mu/max.mu) %>%
        dplyr::mutate(mu = base::sign(mu) * base::log(base::abs(mu)))

      # Cluster the transformed mean values
      kmeans.clusters <- segmentation.results$mu %>%
        stats::kmeans(centers = base::c(base::min(segmentation.results$mu), 0, base::max(segmentation.results$mu)), nstart = 10)

      # Retain only deletions and duplications
      CNV.labels <- kmeans.clusters$centers %>%
        base::as.data.frame() %>%
        tibble::rownames_to_column() %>%
        stats::setNames(base::c('cluster', 'centers')) %>%
        dplyr::arrange(centers) %>%
        dplyr::mutate(class = base::c('Normal+', 'CNV', 'Normal-')) %>%
        dplyr::arrange(cluster) %>%
        dplyr::select(class) %>%
        base::unlist()
      segmentation.results$cluster <- kmeans.clusters %$% cluster %>% base::factor(levels = base::c(1, 2, 3), labels = CNV.labels)
      potential.CNVs <- segmentation.results %>% dplyr::filter(cluster %in% base::c('CNV'))
      potential.CNVs$cluster <- base::droplevels(potential.CNVs$cluster)

      # Separate deletions from duplications based on the original mean
      potential.CNVs <- potential.CNVs %>% dplyr::rowwise() %>% dplyr::mutate(class = if (original.mu < 0) "DELETION" else "DUPLICATION")

    }
  }

  if (base::nrow(potential.CNVs) == 0) {

    if (verbose) base::message("No CNVs found!")
    return(base::invisible(NULL))

  } else {

    # For each observed segment length, generate a null distribution using bootstrap
    bins <- copy.number %>% base::as.data.frame() %>%
      dplyr::select(-seqnames, -start, -end, -width, -strand) %>%
      base::sapply(function(i) i) %>% base::as.vector()
    all.lengths <- base::sort(base::unique(potential.CNVs$num.mark))
    cl <- parallel::makeCluster(4)
    doParallel::registerDoParallel(cl)
    mu.null.distribution <- foreach::foreach(segment.length = all.lengths, .packages = base::c('magrittr', 'dplyr')) %dopar% {
      number.of.segments <- potential.CNVs %>%
        dplyr::filter(num.mark == segment.length) %>%
        base::nrow()
      real.segment.statistics <- potential.CNVs %>%
        dplyr::filter(num.mark == segment.length) %>%
        dplyr::select(original.mu) %>%
        base::unlist()
      fake.segment.statistics <- base::replicate(1000 - number.of.segments, {
        repeat {
          fake.segment <- base::sample(bins, size = segment.length, replace = TRUE)
          if (base::sum(base::is.na(fake.segment) == FALSE) > 1) break
        }
        base::mean(fake.segment, na.rm = TRUE)
      })
      base::c(fake.segment.statistics, real.segment.statistics) %>% base::sort()
    }
    parallel::stopCluster(cl)
    base::names(mu.null.distribution) <- base::paste0('X', all.lengths)

    # Calculate one-sided p-values for the identified potential CNVs
    p.value.deletions <- function(distribution, observation) {
      return(base::mean(base::round(observation, 6) >= base::round(distribution, 6)))
    }
    p.value.duplications <- function(distribution, observation) {
      return(base::mean(base::round(observation, 6) <= base::round(distribution, 6)))
    }
    potential.CNVs <- potential.CNVs %>% dplyr::rowwise() %>%
      dplyr::mutate(pvalue = dplyr::case_when(
        class == "DELETION" ~ p.value.deletions(mu.null.distribution[[base::paste0('X', num.mark)]], original.mu),
        class == "DUPLICATION" ~ p.value.duplications(mu.null.distribution[[base::paste0('X', num.mark)]], original.mu))) %>%
      base::as.data.frame()

  }

  if (verbose) base::message('Generating the output...')

  if (plot.results) {

    # Generate per-chromosome plots
    base::suppressWarnings(base::dir.create(base::paste0(output.dir, '/Plots')))
    segmentation.results <- segmentation.results %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    potential.CNVs <- potential.CNVs %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    foreach::foreach(chr = base::unique(segmentation.results@seqnames) %>% base::as.vector()) %dopar% {
      base::source(base::paste0(working.dir, '/plotChromSegments.R'))
      grDevices::pdf(base::paste0(output.dir, '/Plots/Chromosome', chr, '.pdf'), width = 12, height = 9)
      for (sample.id in base::unique(segmentation.results$ID)) {
        plotChromSegments(sample.id,
                          CN.values = copy.number[copy.number@seqnames == chr]@elementMetadata[sample.id] %>% base::unlist() %>% base::as.vector(),
                          CN.positions = copy.number[copy.number@seqnames == chr]@ranges@start,
                          segmentation.results = segmentation.results[segmentation.results@seqnames == chr & segmentation.results$ID == sample.id],
                          potential.CNVs = potential.CNVs[potential.CNVs@seqnames == chr & potential.CNVs$ID == sample.id], chr = chr,
                          seg.duplications = seg.duplications, seg.duplications.threshold = seg.duplications.threshold)
      }
      grDevices::dev.off()
    }
    parallel::stopCluster(cl)

  }

  # Flag CNVs that overlap with large segmental duplications
  potential.CNVs <- potential.CNVs %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  overlapping <- IRanges::subsetByOverlaps(potential.CNVs, seg.duplications[seg.duplications@ranges@width > seg.duplications.threshold], ignore.strand = TRUE)
  GenomicRanges::elementMetadata(potential.CNVs)[['overlaps.seg.duplication']] <- base::factor(potential.CNVs %in% overlapping,
                                                                                               levels = base::c(FALSE, TRUE), labels = base::c('No','Yes'))

  # For each potential CNV, list any syndromes that could be caused by that CNV
  GenomicRanges::elementMetadata(potential.CNVs)[['syndrome']] <- NA
  overlapping <- base::suppressWarnings(IRanges::findOverlaps(syndromes, potential.CNVs, ignore.strand = TRUE)) %>%
    base::as.data.frame() %>%
    dplyr::group_by(subjectHits) %>%
    dplyr::summarise(syndrome = base::list(queryHits))
  potential.CNVs[overlapping$subjectHits]$syndrome <- overlapping$syndrome
  potential.CNVs <- base::as.data.frame(potential.CNVs) %>%
    dplyr::mutate(syndromic = !base::is.na(syndrome)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(syndrome = base::ifelse(syndromic, base::list(base::as.data.frame(syndromes[syndrome]) %>% dplyr::filter(base::as.character(type) == base::as.character(class)) %>% dplyr::select(syndrome)), NA)) %>%
    dplyr::mutate(syndrome = base::paste(base::unlist(syndrome), collapse = ', '))
  potential.CNVs$syndrome[potential.CNVs$syndrome %in% base::c('character(0)', 'NA')] <- ''

  # Calculate the occurrence of each variant
  potential.CNVs <- potential.CNVs %>%
    dplyr::group_by(seqnames) %>%
    dplyr::group_modify(~ {
      CNV.dist <- .x %>% dplyr::select(start, end) %>% stats::dist(method = 'manhattan')
      CNV.dist[] <- base::c(1, 0)[(CNV.dist > (occurrence.margin * bin.size)) + 1]
      .x %>%
        dplyr::mutate(clusters = igraph::graph_from_adjacency_matrix(CNV.dist, mode = "undirected") %>% igraph::clusters() %$% membership) %>%
        dplyr::group_by(clusters) %>%
        dplyr::mutate(occurrence = base::length(base::unique(ID))) %>%
        dplyr::ungroup()
    }) %>%
    dplyr::ungroup() %>%
    dplyr::select(-clusters)

  # Pre-sort all potential CNVs according to their significance
  potential.CNVs <- base::as.data.frame(potential.CNVs) %>%
    dplyr::arrange(pvalue, dplyr::desc(syndrome)) %>%
    dplyr::select(seqnames, start, end, width, ID, seg.mean, class, overlaps.seg.duplication, occurrence, syndrome, pvalue)

  # Final output
  base::colnames(potential.CNVs) <- base::c('Chromosome',
                                            'Start', 'End',
                                            'Width',
                                            'ID',
                                            'Copy Number log2 Mean Ratio',
                                            'CNV type',
                                            base::paste0('Overlaps with a segmental duplication >', seg.duplications.threshold / 1e3, 'kb'),
                                            'Occurrence',
                                            'Associated syndromes',
                                            'P-value')
  base::suppressMessages(openxlsx::write.xlsx(potential.CNVs, base::paste0(output.dir, '/Variations.xlsx')))
}
