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


base::suppressMessages(base::library(GenomeInfoDb))
base::suppressMessages(base::library(CopywriteR))
base::suppressMessages(base::library(fastseg))
base::suppressMessages(base::library(GenomicRanges))
base::suppressMessages(base::library(IRanges))

base::source(base::paste0(working.dir, '/parSort.R'))
base::source(base::paste0(working.dir, '/parCount.R'))
base::source(base::paste0(working.dir, '/correct.R'))
base::source(base::paste0(working.dir, '/parSegment.R'))
base::source(base::paste0(working.dir, '/postSegment.R'))

# Prepare bins with AT content, GC content and mappability metadata
base::invisible(utils::capture.output(CopywriteR::preCopywriteR(output.dir,
                                                                bin.size = bin.size,
                                       ref.genome = reference)))
base::load(base::paste0(output.dir, '/', reference, '_', bin.size.text, '/GC_mappability.rda'))
GenomeInfoDb::seqlevelsStyle(GC.mappa.grange) <- format
bins <- base::sort(GC.mappa.grange[GC.mappa.grange@seqnames %in% autosomes])
base::colnames(bins@elementMetadata) <- base::c('at', 'gc', 'map')

# Sort and index the input BAM files if necessary
if (skip.sorting.BAM) {
  if (verbose) base::message("Skipped sorting of the input files.")
  sorted.bam.files <- base::list.files(base::file.path(bam.files.dir),
                                       pattern = '.bam$',
                                       full.names = TRUE)
} else {
  if (verbose) base::message("Sorting the input files...")
  sorted.bam.files <- bam.files.dir %>% parSort(ncores = ncores)
}

if (verbose) {
  base::message("Analysing the following BAM files:")
  base::invisible(utils::capture.output(base::sapply(sorted.bam.files, message)))
}

# Either use existing read counts or calculate new ones
if (base::file.exists(path.to.counts)) {
  if (verbose) base::message("Using existing read counts.")
  counts <- base::readRDS(path.to.counts)
} else {
  if (verbose) base::message("Counting reads in the bins created...")
  counts <- parCount(sorted.bam.files, bins,
                     min.mapq = min.mapq,
                     skip.indexing.BAM = skip.indexing.BAM,
                     ncores = ncores)
  base::saveRDS(counts, base::paste0(output.dir, '/counts.rds'))
  if (verbose) base::message("Finished counting the reads.")
}

# Correct raw read counts for GC content and mappability biases and convert them to scaled copy number
if (verbose) base::message("Correcting the read counts...")
copy.number <- base::suppressWarnings(correct(bins, counts, autosomes, min.mapq, rough.span, final.span))
GenomeInfoDb::seqlevelsStyle(copy.number) <- 'NCBI'
if (verbose) base::message("Finished correcting the read counts.")

# For each chromosome, exclude all bins that overlap with the centromere
centr.with.margin <- base::as.data.frame(centromeres)
centr.with.margin$start <- centr.with.margin$start - centromeres.margin
centr.with.margin$end <- centr.with.margin$end + centromeres.margin
centr.with.margin <- GenomicRanges::makeGRangesFromDataFrame(centr.with.margin)
overlapping <- IRanges::subsetByOverlaps(copy.number, centr.with.margin, ignore.strand = TRUE)
GenomicRanges::elementMetadata(copy.number)[['centromere']] <- copy.number %in% overlapping
copy.number <- copy.number %>% base::as.data.frame() %>%
  dplyr::mutate_at(.vars = dplyr::vars(-seqnames, -start, -end, -width, -strand, -centromere), .funs = base::list(~ base::ifelse(centromere, NA, .))) %>%
  dplyr::select(-centromere) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Obtain continuous copy-number segments
if (verbose) base::message('Starting segmentation...')
segmentation.results <- parSegment(copy.number, bin.size, ncores)
if (verbose) base::message('Segmentation finished.')

# Analyse the obtained segments and summarise the results
postSegment(copy.number, segmentation.results,
            del.threshold, dup.threshold,
            seg.duplications, seg.duplications.threshold,
            syndromes,
            working.dir, output.dir,
            occurrence.margin, bin.size,
            verbose, plot.results, ncores)

if (verbose) base::message('Finished!')
base::writeLines(utils::capture.output(utils::sessionInfo()), base::paste0(output.dir, '/sessionInfo.txt'))
