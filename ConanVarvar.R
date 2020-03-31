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


base::suppressMessages(base::library(optparse))

option_list = base::list(
  optparse::make_option("--bamdir",
                        type = "character",
                        action = "store",
                        default = NA,
                        help = "Full path to the directory with BAM files to be processed",
                        metavar = "/path/to/data"),
  optparse::make_option("--counts",
                        type = "character",
                        action = "store",
                        default = "/path/to/counts.rds",
                        help = "Full path to read counts from previous runs",
                        metavar = "/path/to/counts.rds"),
  optparse::make_option("--reference",
                        type = "character",
                        default = "hg38",
                        help = "Reference genome: hg38 or hg19 [default %default]",
                        metavar = "hg38"),
  optparse::make_option("--format",
                        type = "character",
                        default = "UCSC",
                        help = "Format of input sequences: UCSC or NCBI [default %default]",
                        metavar = "UCSC"),
  optparse::make_option("--sortbam",
                        action = "store_true",
                        default = FALSE,
                        help = "Sort the input BAM files"),
  optparse::make_option("--indexbam",
                        action = "store_true",
                        default = FALSE,
                        help = "Index the input BAM files"),
  optparse::make_option("--outdir",
                        type = "character",
                        action = "store",
                        default = NA,
                        help = "Output directory",
                        metavar = "/output"),
  optparse::make_option("--binsize",
                        default = base::as.integer(50e3),
                        help = "Bin size [default %default]",
                        metavar = "50000"),
  optparse::make_option("--minmapq",
                        default = base::as.double(0.8),
                        help = "Minimum mapping quality of reads [default %default]",
                        metavar = "0.8"),
  optparse::make_option("--roughspan",
                        default = base::as.double(0.1),
                        help = "First (rough) Loess span [default %default]",
                        metavar = "0.1"),
  optparse::make_option("--finalspan",
                        default = base::as.double(0.3),
                        help = "Second (final) Loess span [default %default]",
                        metavar = "0.3"),
  optparse::make_option("--delthresh",
                        default = base::as.double(-0.7),
                        help = "Threshold for deletions [default %default]",
                        metavar = "-0.7"),
  optparse::make_option("--dupthresh",
                        default = base::as.double(0.5),
                        help = "Threshold for duplications [default %default]",
                        metavar = "0.5"),
  optparse::make_option("--centrmargin",
                        default = base::as.integer(500e3),
                        help = "Centromeres margin [default %default]",
                        metavar = "500000"),
  optparse::make_option("--segdupthresh",
                        default = base::as.integer(100e3),
                        help = "Lower threshold for segmental duplications [default %default]",
                        metavar = "100000"),
  optparse::make_option("--occurrencemargin",
                        default = base::as.integer(6),
                        help = "Start/End margin of error for the Occurrence",
                        metavar = "6"),
  optparse::make_option("--ncores",
                        default = base::as.integer(4),
                        help = "Number of cores for parallelisation [default %default]",
                        metavar = "4"),
  optparse::make_option("--plotresults",
                        action = "store_true",
                        default = FALSE,
                        help = "Plot the results"),
  optparse::make_option("--verbose",
                        action = "store_true",
                        default = TRUE,
                        help = "Print extra output [default]")
)

opts = optparse::parse_args(optparse::OptionParser(option_list = option_list))

if (base::is.na(opts$bamdir)) {
  base::stop("Please specify the directory with input BAM files.")
}
if (base::is.na(opts$outdir)) {
  base::stop("Output directory is required.")
}
if (!(base::tolower(opts$reference) %in% base::c("hg38", "hg19"))) {
  base::stop("The reference has to be either hg38 or hg19.")
}
if (!(opts$format %in% base::c("UCSC", "NCBI"))) {
  base::stop("The format has to be either UCSC or NCBI.")
}
if (!(opts$binsize %in% base::c(1e3, 10e3, 50e3, 100e3, 200e3))) {
  base::stop("Please set the bin size to one of the following values: 1kb, 10kb, 50kb, 100kb, 200kb.")
}
if (!(opts$minmapq %in% base::seq(0.1, 1.0, 0.1))) {
  base::stop("Minimum mapping quality has to be a value between 0.1 and 1.0 with one decimal point.")
}
if (!(base::all(base::c(opts$roughspan, opts$finalspan) %in% base::seq(0.05, 0.50, 0.05)))) {
  base::stop("At least one of the specified Loess span values is not valid.")
}
if (opts$delthresh >= 0) {
  base::stop("Threshold for deletions has to be a negative number.")
}
if (opts$dupthresh <= 0) {
  base::stop("Threshold for duplications has to be a positive number.")
}
if (!(opts$centrmargin %in% base::c(0, 10e3, 100e3, 500e3))) {
  base::stop("Please set the centromeres margin to one of the following values: 0, 10kb, 100kb, 500kb.")
}
if (!(opts$segdupthresh %in% base::c(0, 10e3, 50e3, 100e3))) {
  base::stop("Please set the lower threshold for segmental duplications to one of the following values: 0, 10kb, 50kb, 100kb.")
}
if (!(opts$occurrencemargin %in% base::seq(0, 10, 1))) {
  base::stop("The margin has to be an integer between 1 and 10.")
}
if (!(opts$ncores %in% base::seq(4, 28, 4))) {
  base::stop("Please use 4, 8, 12, 16, 20, 24 or 28 cores.")
}

base::suppressMessages(base::library(GenomeInfoDb))
base::suppressMessages(base::library(rCGH))

working.dir <<- base::getwd()
output.dir <<- opts$outdir
bam.files.dir <<- opts$bamdir

path.to.counts <<- opts$counts

skip.sorting.BAM <<- !opts$sortbam
skip.indexing.BAM <<- !opts$indexbam

plot.results <<- opts$plotresults
verbose <<- opts$verbose
ncores <<- opts$ncores

bin.size <<- opts$binsize
bin.size.text <<- base::ifelse(bin.size == 0, "0", base::paste0(base::as.character(bin.size / 1e3), "kb"))

reference <<- opts$reference
seqnames <<- base::readRDS(base::paste0(reference, '.seqnames.rds'))
base::assign(
  x = "format",
  value = opts$format,
  envir = .GlobalEnv
)
GenomeInfoDb::seqlevelsStyle(seqnames) <- opts$format
base::assign(
  x = "seqnames",
  value = seqnames,
  envir = .GlobalEnv
)

centromeres <<- base::readRDS(base::paste0(reference, '.centromeres.rds'))
syndromes <<- base::readRDS(base::paste0(reference, '.syndromes.rds'))

if (opts$format == 'UCSC'){
  autosomes <<- base::c(base::paste0('chr', 1:22))
} else {
  autosomes <<- base::c(base::paste0('', 1:22))
}

min.mapq <<- opts$minmapq

rough.span <<- opts$roughspan
final.span <<- opts$finalspan

del.threshold <<- opts$delthresh
dup.threshold <<- opts$dupthresh

seg.duplications <<- base::readRDS(base::paste0(reference, '.segmental.duplications.rds'))

centromeres.margin <<- opts$centrmargin
seg.duplications.threshold <<- opts$segdupthresh

occurrence.margin <<- opts$occurrencemargin

source("main.R")
