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


suppressMessages(library(ggplot2))
suppressMessages(library(latex2exp))
suppressMessages(library(IRanges))
suppressMessages(library(GenomicRanges))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))

# This function generates per-chromosome plots.
plotChromSegments <- function(sample.id, CN.values, CN.positions, segmentation.results, potential.CNVs = NULL,
                              chr, seg.duplications, seg.duplications.threshold) {

  # List of segmental duplications for chrom with the length above the user-specified threshold
  seg.dups = seg.duplications[seg.duplications@seqnames == chr & seg.duplications@ranges@width > seg.duplications.threshold]

  segments = data.frame(
    x = segmentation.results@ranges@start,
    xend = segmentation.results@ranges@start + segmentation.results@ranges@width - 1,
    y = segmentation.results$seg.mean,
    yend = segmentation.results$seg.mean)

  CN.data = data.frame(value = CN.values, position = CN.positions)

  # Indicate which copy number values are part of a segmental duplication
  CN.data$type = 'Normal'
  CN.data$type[sapply(CN.data$position,
                      function(pos) GRanges(chr, IRanges(pos, pos)) %>%
                                      subsetByOverlaps(seg.dups) %>% length() > 0)] = 'Part of Segmental Duplication'

  colors = setNames(c("#9999CC", "black", "#07FDB5", "#FD074F"),
                    c('Part of Segmental Duplication', 'Normal', 'DUPLICATION', 'DELETION'))

  # Plot copy number values
  plot = ggplot(CN.data, aes(x = position, y = value)) +
    geom_point(aes(colour = type)) +
    geom_segment(data = segments,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = '#77686a', size = 2) +
    ggtitle(paste0(sample.id, ', chromosome ', chr)) +
    ylab(TeX('Copy Number $\\log_2$ Ratio')) +
    xlab('Chromosome Position, Mbp') +
    scale_x_continuous(labels = function(x) x / 1e6)

  # If there are any deletions/duplications, add them to the plot
  if (!missing(potential.CNVs)) {

    segments = data.frame(
      x = potential.CNVs@ranges@start,
      xend = potential.CNVs@ranges@start + potential.CNVs@ranges@width - 1,
      y = potential.CNVs$seg.mean,
      yend = potential.CNVs$seg.mean,
      call = potential.CNVs$class)

    plot = plot +
      geom_segment(data = segments,
                   aes(x = x, y = y, xend = xend, yend = yend, colour = call), size = 1)
  }

  # Add extra features to the plot
  plot = plot +
    scale_color_manual(values = colors) +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 13),
          axis.title = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(size = 14))
  print(plot)
}
