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


suppressMessages(library(GenomeInfoDb))
suppressMessages(library(CopywriteR))
suppressMessages(library(fastseg))
suppressMessages(library(GenomicRanges))
suppressMessages(library(IRanges))

source(paste0(working.dir, '/parSort.R'))
source(paste0(working.dir, '/parCount.R'))
source(paste0(working.dir, '/correct.R'))
source(paste0(working.dir, '/parSegment.R'))
source(paste0(working.dir, '/postSegment.R'))

# Prepare bins with AT content, GC content and mappability metadata
invisible(capture.output(preCopywriteR(output.dir,
                                       bin.size = bin.size,
                                       ref.genome = reference)))
load(paste0(output.dir, '/', reference, '_', bin.size.text, '/GC_mappability.rda'))
seqlevelsStyle(GC.mappa.grange) <- format
bins <- sort(GC.mappa.grange[GC.mappa.grange@seqnames %in% autosomes])
colnames(bins@elementMetadata) <- c('at', 'gc', 'map')

# Sort and index the input BAM files if necessary
if (skip.sorting.BAM) {
  if (verbose) message("Skipped sorting of the input files.")
  sorted.bam.files <- list.files(file.path(bam.files.dir),
                                 pattern = '.bam$',
                                 full.names = TRUE)
} else {
    if (verbose) message("Sorting the input files...")
    sorted.bam.files <- bam.files.dir %>% parSort(ncores = ncores)
}

if (verbose) {
    message("Analysing the following BAM files:")
    invisible(capture.output(sapply(sorted.bam.files, message)))
}

# Either use existing read counts or calculate new ones
if (file.exists(path.to.counts)) {
  if (verbose) message("Using existing read counts.")
  counts <- readRDS(path.to.counts)
} else {
  if (verbose) message("Counting reads in the bins created...")
  counts <- parCount(sorted.bam.files, bins,
                     min.mapq = min.mapq,
                     skip.indexing.BAM = skip.indexing.BAM,
                     ncores = ncores)
  saveRDS(counts, paste0(output.dir, '/counts.rds'))
  if (verbose) message("Finished counting the reads.")
}

# Correct raw read counts for GC content and mappability biases and convert them to scaled copy number
if (verbose) message("Correcting the read counts...")
copy.number <- suppressWarnings(correct(bins, counts, autosomes, min.mapq, rough.span, final.span))
seqlevelsStyle(copy.number) <- 'NCBI'
if (verbose) message("Finished correcting the read counts.")

# For each chromosome, exclude all bins that overlap with the centromere
centr.with.margin <- as.data.frame(centromeres)
centr.with.margin$start <- centr.with.margin$start - centromeres.margin
centr.with.margin$end <- centr.with.margin$end + centromeres.margin
centr.with.margin <- makeGRangesFromDataFrame(centr.with.margin)
overlapping <- subsetByOverlaps(copy.number, centr.with.margin, ignore.strand = TRUE)
elementMetadata(copy.number)[['centromere']] <- copy.number %in% overlapping
copy.number <- copy.number %>% as.data.frame() %>%
  dplyr::mutate_at(.vars = vars(-seqnames, -start, -end, -width, -strand, -centromere), .funs = list(~ ifelse(centromere, NA, .))) %>%
  select(-centromere) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Obtain continuous copy-number segments
if (verbose) message('Starting segmentation...')
segmentation.results <- parSegment(copy.number, bin.size, ncores)
if (verbose) message('Segmentation finished.')

# Analyse the obtained segments and summarise the results
postSegment(copy.number, segmentation.results,
            del.threshold, dup.threshold,
            seg.duplications, seg.duplications.threshold,
            syndromes,
            working.dir, output.dir,
            verbose, plot.results, ncores)

if (verbose) message('Finished!')
writeLines(capture.output(sessionInfo()), paste0(output.dir, '/sessionInfo.txt'))
