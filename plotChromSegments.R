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
# See the GNU General Public License for more details <http://www.gnu.org/licenses/>.


base::suppressMessages(base::library(ggplot2))
base::suppressMessages(base::library(IRanges))
base::suppressMessages(base::library(GenomicRanges))
base::suppressMessages(base::library(magrittr))
base::suppressMessages(base::library(dplyr))

# This function generates per-chromosome plots.
plotChromSegments <- function(sample.id, CN.values, CN.positions, segmentation.results, potential.CNVs = NULL,
                              chr, seg.duplications, seg.duplications.threshold) {

  # List of segmental duplications for chrom with the length above the user-specified threshold
  seg.dups = seg.duplications[seg.duplications@seqnames == chr & seg.duplications@ranges@width > seg.duplications.threshold]

  segments = base::data.frame(
    x = segmentation.results@ranges@start,
    xend = segmentation.results@ranges@start + segmentation.results@ranges@width - 1,
    y = segmentation.results$seg.mean,
    yend = segmentation.results$seg.mean)

  CN.data = base::data.frame(value = CN.values, position = CN.positions)

  # Indicate which copy number values are part of a segmental duplication
  CN.data$type = 'Normal'
  CN.data$type[base::sapply(CN.data$position,
                            function(pos) GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos)) %>%
                                            IRanges::subsetByOverlaps(seg.dups) %>% base::length() > 0)] = 'Part of Segmental Duplication'

  colors = stats::setNames(base::c("#9999CC", "black", "#07FDB5", "#FD074F"),
                           base::c('Part of Segmental Duplication', 'Normal', 'DUPLICATION', 'DELETION'))

  # Plot copy number values
  plot = ggplot2::ggplot(CN.data, ggplot2::aes(x = position, y = value)) +
    ggplot2::geom_point(ggplot2::aes(colour = type)) +
    ggplot2::geom_segment(data = segments,
                          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          color = '#77686a', size = 2) +
    ggplot2::ggtitle(base::paste0(sample.id, ', chromosome ', chr)) +
    ggplot2::ylab('Copy Number log2 Mean Ratio') +
    ggplot2::xlab('Chromosome Position, Mbp') +
    ggplot2::scale_x_continuous(labels = function(x) x / 1e6)

  # If there are any deletions/duplications, add them to the plot
  if (!missing(potential.CNVs)) {

    segments = base::data.frame(
      x = potential.CNVs@ranges@start,
      xend = potential.CNVs@ranges@start + potential.CNVs@ranges@width - 1,
      y = potential.CNVs$seg.mean,
      yend = potential.CNVs$seg.mean,
      call = potential.CNVs$class)

    plot = plot +
      ggplot2::geom_segment(data = segments,
                            ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = call), size = 1)
  }

  # Add extra features to the plot
  plot = plot +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 13),
                   axis.title = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(size = 12),
                   axis.text.y = ggplot2::element_text(size = 12),
                   plot.title = ggplot2::element_text(size = 14))
  print(plot)
}
