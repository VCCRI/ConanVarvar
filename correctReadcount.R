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

  if (length(x$reads) == 0 | length(x$gc) == 0 | length(x$map) == 0) {
    stop("Missing one of required columns: reads, gc, map")
  }

  if (verbose) message("Applying filter on the data...")
  x$valid <- TRUE
  x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
  x$ideal <- TRUE
  routlier <- 0.01
  range <- quantile(x$reads[x$valid], prob = c(0, 1 - routlier), na.rm = TRUE)
  doutlier <- 0.001
  domain <- quantile(x$gc[x$valid],
                     prob = c(doutlier, 1 - doutlier),
                     na.rm = TRUE)
  x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] |
            x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
  if (verbose) message("Correcting for the GC bias...")
  set <- which(x$ideal)
  select <- sample(set, length(set))
  rough = loess(x$reads[select] ~ x$gc[select], span = rough.span)
  i <- seq(0, 1, by = 0.001)
  final = loess(predict(rough, i) ~ i, span = final.span)
  x$cor.gc <- x$reads/predict(final, x$gc)
  if (verbose) message("Correcting for the mappability bias...")
  coutlier <- 0.01
  range <- quantile(x$cor.gc[which(x$valid)],
                    prob = c(0, 1 - coutlier),
                    na.rm = TRUE)
  set <- which(x$cor.gc < range[2])
  select <- sample(set, length(set))
  final = approxfun(lowess(x$map[select], x$cor.gc[select]))
  x$cor.map <- x$cor.gc/final(x$map)
  x$copy <- x$cor.map
  x$copy[x$copy <= 0] = NA
  x$copy <- log(x$copy, 2)
  return(x)
}
