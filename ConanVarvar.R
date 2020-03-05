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


suppressMessages(library(optparse))

option_list = list(
  make_option("--bamdir",
              type = "character",
              action = "store",
              default = NA,
              help = "Full path to the directory with BAM files to be processed",
              metavar = "/path/to/data"),
  make_option("--counts",
              type = "character",
              action = "store",
              default = "/path/to/counts.rds",
              help = "Full path to read counts from previous runs",
              metavar = "/path/to/counts.rds"),
  make_option("--reference",
              type = "character",
              default = "hg38",
              help = "Reference genome: hg38 or hg19 [default %default]",
              metavar = "hg38"),
  make_option("--format",
              type = "character",
              default = "UCSC",
              help = "Format of input sequences: UCSC or NCBI [default %default]",
              metavar = "UCSC"),
  make_option("--sortbam",
              action = "store_true",
              default = FALSE,
              help = "Sort the input BAM files"),
  make_option("--indexbam",
              action = "store_true",
              default = FALSE,
              help = "Index the input BAM files"),
  make_option("--outdir",
              type = "character",
              action = "store",
              default = NA,
              help = "Output directory",
              metavar = "/output"),
  make_option("--binsize",
              default = as.integer(50e3),
              help = "Bin size [default %default]",
              metavar = "50000"),
  make_option("--minmapq",
              default = as.double(0.8),
              help = "Minimum mapping quality of reads [default %default]",
              metavar = "0.8"),
  make_option("--roughspan",
              default = as.double(0.1),
              help = "First (rough) Loess span [default %default]",
              metavar = "0.1"),
  make_option("--finalspan",
              default = as.double(0.3),
              help = "Second (final) Loess span [default %default]",
              metavar = "0.3"),
  make_option("--delthresh",
              default = as.double(-0.7),
              help = "Threshold for deletions [default %default]",
              metavar = "-0.7"),
  make_option("--dupthresh",
              default = as.double(0.5),
              help = "Threshold for duplications [default %default]",
              metavar = "0.5"),
  make_option("--centrmargin",
              default = as.integer(500e3),
              help = "Centromeres margin [default %default]",
              metavar = "500000"),
  make_option("--segdupthresh",
              default = as.integer(100e3),
              help = "Lower threshold for segmental duplications [default %default]",
              metavar = "100000"),
  make_option("--ncores",
              default = as.integer(4),
              help = "Number of cores for parallelisation [default %default]",
              metavar = "4"),
  make_option("--plotresults",
              action = "store_true",
              default = FALSE,
              help = "Plot the results"),
  make_option("--verbose",
              action = "store_true",
              default = TRUE,
              help = "Print extra output [default]")
)

opts = parse_args(OptionParser(option_list = option_list))

if (is.na(opts$bamdir)) {
  stop("Please specify the directory with input BAM files.")
}
if (is.na(opts$outdir)) {
  stop("Output directory is required.")
}
if (!(tolower(opts$reference) %in% c("hg38", "hg19"))) {
  stop("The reference has to be either hg38 or hg19.")
}
if (!(opts$format %in% c("UCSC", "NCBI"))) {
  stop("The format has to be either UCSC or NCBI.")
}
if (!(opts$binsize %in% c(1e3, 10e3, 50e3, 100e3, 200e3))) {
  stop("Please set the bin size to one of the following values: 1kb, 10kb, 50kb, 100kb, 200kb.")
}
if (!(opts$minmapq %in% seq(0.1, 1.0, 0.1))) {
  stop("Minimum mapping quality has to be a value between 0.1 and 1.0 with one decimal point.")
}
if (!(all(c(opts$roughspan, opts$finalspan) %in% seq(0.05, 0.50, 0.05)))) {
  stop("At least one of the specified Loess span values is not valid.")
}
if (opts$delthresh >= 0) {
  stop("Threshold for deletions has to be a negative value.")
}
if (opts$dupthresh <= 0) {
  stop("Threshold for duplications has to be a positive value.")
}
if (!(opts$centrmargin %in% c(0, 10e3, 100e3, 500e3))) {
  stop("Please set the centromeres margin to one of the following values: 0, 10kb, 100kb, 500kb.")
}
if (!(opts$segdupthresh %in% c(0, 10e3, 50e3, 100e3))) {
  stop("Please set the lower threshold for segmental duplications to one of the following values: 0, 10kb, 50kb, 100kb.")
}
if (!(opts$ncores %in% seq(4, 28, 4))) {
  stop("Please use 4, 8, 12, 16, 20, 24 or 28 cores.")
}

suppressMessages(library(GenomeInfoDb))
suppressMessages(library(rCGH))

working.dir <<- getwd()
output.dir <<- opts$outdir
bam.files.dir <<- opts$bamdir

path.to.counts <<- opts$counts

skip.sorting.BAM <<- !opts$sortbam
skip.indexing.BAM <<- !opts$indexbam

plot.results <<- opts$plotresults
verbose <<- opts$verbose
ncores <<- opts$ncores

bin.size <<- opts$binsize
bin.size.text <<- ifelse(bin.size == 0, "0", paste0(as.character(bin.size / 1e3), "kb"))

reference <<- opts$reference
seqnames <- readRDS(paste0(reference, '.seqnames.rds'))
assign(
  x = "format",
  value = opts$format,
  envir = .GlobalEnv
)
seqlevelsStyle(seqnames) <- opts$format
assign(
  x = "seqnames",
  value = seqnames,
  envir = .GlobalEnv
)

centromeres <<- readRDS(paste0(reference, '.centromeres.rds'))
syndromes <<- readRDS(paste0(reference, '.syndromes.rds'))

if (opts$format == 'UCSC'){
  autosomes <<- c(paste0('chr', 1:22))
} else {
  autosomes <<- c(paste0('', 1:22))
}

min.mapq <<- opts$minmapq

rough.span <<- opts$roughspan
final.span <<- opts$finalspan

del.threshold <<- opts$delthresh
dup.threshold <<- opts$dupthresh

seg.duplications <<- readRDS(paste0(reference, '.segmental.duplications.rds'))

centromeres.margin <<- opts$centrmargin
seg.duplications.threshold <<- opts$segdupthresh

source("main.R")
