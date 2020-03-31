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

# NOTE:
# This is a modified version of the analogous function from the "HMMcopy" package

correctReadcount <- function(x, mappability = 0.8,
                             rough.span = 0.1, final.span = 0.3,
                             verbose = TRUE) {

  if (base::length(x$reads) == 0 | base::length(x$gc) == 0 | base::length(x$map) == 0) {
    base::stop("Missing one of required columns: reads, gc, map")
  }

  if (verbose) base::message("Applying filter on the data...")
  x$valid <- TRUE
  x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
  x$ideal <- TRUE
  routlier <- 0.01
  range <- stats::quantile(x$reads[x$valid], prob = base::c(0, 1 - routlier), na.rm = TRUE)
  doutlier <- 0.001
  domain <- stats::quantile(x$gc[x$valid],
                            prob = base::c(doutlier, 1 - doutlier),
                            na.rm = TRUE)
  x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] |
            x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
  if (verbose) base::message("Correcting for the GC bias...")
  set <- base::which(x$ideal)
  select <- base::sample(set, base::length(set))
  rough = stats::loess(x$reads[select] ~ x$gc[select], span = rough.span)
  i <- base::seq(0, 1, by = 0.001)
  final = stats::loess(stats::predict(rough, i) ~ i, span = final.span)
  x$cor.gc <- x$reads / stats::predict(final, x$gc)
  if (verbose) base::message("Correcting for the mappability bias...")
  coutlier <- 0.01
  range <- stats::quantile(x$cor.gc[base::which(x$valid)],
                           prob = base::c(0, 1 - coutlier),
                           na.rm = TRUE)
  set <- base::which(x$cor.gc < range[2])
  select <- base::sample(set, base::length(set))
  final = stats::approxfun(stats::lowess(x$map[select], x$cor.gc[select]))
  x$cor.map <- x$cor.gc / final(x$map)
  x$copy <- x$cor.map
  x$copy[x$copy <= 0] = NA
  x$copy <- base::log(x$copy, 2)
  return(x)
}
