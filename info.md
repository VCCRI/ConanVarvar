# ConanVarvar

## Overview
This program is a versatile tool for detecting copy number variations (CNVs) from BAM files. The tool was designed specifically for the detection of *large* (>1 Mb) CNVs in batches of multiple samples.

## How to use

Please, set up all the input parameters on the left-hand side and click `Run` below.
You will be notified when execution is finished in the console below.

## Warnings
* If you wish to sort and/or index input BAM files, make sure you have writing permission on that directory.
* UCSC format has a "chr" in the representation of the chromosomes, whereas the NCBI format does not.
* Mapping Quality sets a threshold for the mapping quality of reads. That is, if this parameter is set to 1 (i.e. 100%), only fully-mapped reads will be considered.
* Please note that LOESS span values are dependent on the number of data points to be processed in the correction step. For example, changing the resolution from 100kb to 50kb will increase the total number of data points, so LOESS span should be adjusted accordingly.
* Centromeres Margin is referred to a region outside the centromere's position that is believed to, nevertheless, be affected by the centromere (e.g. abnormal copy number can be observed in that region).
The region between the centromere's left and right positions is excluded from the final CNV analysis to reduce the number of possible false positives. Therefore, if Centromeres Margin is set to some non-zero value, the effect is the same as if the left and right positions were moved further by the specified value.
* High number of cores requested is only reasonable if the number of input files is high. Generally, if *n* files need to be processed, setting the number of cores to *n* will give the best performance results in terms of the execution time.
* Consider setting `Plot the results?` to `No` if you are only interested in the variants' statistics and if you do not need per-chromosome CNV plots, as this can significantly reduce the execution time.
* If you want to re-use existing read counts, please specify a path to the `counts.rds` file and use the same bin size and mapping quality threshold as what you used to generate that file.
