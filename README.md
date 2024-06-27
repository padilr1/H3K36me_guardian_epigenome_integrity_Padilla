Created by Reinnier Padilla.

Copyright (C) 2024 Reinnier Padilla

---

* The following programs are required:
  * [deepTools](https://github.com/deeptools/deepTools)
  * [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks)
  * [featureCounts](https://subread.sourceforge.net/featureCounts.html)
  * [HOMER](http://homer.ucsd.edu/homer/ngs/peakMotifs.html)
  * [R](https://cran.r-project.org/)
    * Bioconductor
      * ChIPseeker
      * DESeq2
      * DiffBind
      * GenomicFeatures
      * GenomicRanges
      * SummarizedExperiment
      * rtracklayer
      * grid
      * gridExtra
    * CRAN
      * MASS
      * cowplot
      * data.table
      * ggrepel
      * ggsignif
      * hexbin
      * isoband
      * lwgeom
      * msigdbr
      * pals
      * patchwork
      * readxl
      * reshape2
      * scales
      * sf
      * tidyverse
      * viridis
* The scripts below should be ran to pre-process a number of files 
  * `generate_norm_bw.sh` provides examples of how to generate the normalized bigWigs used to generate the figures, including depth-normalized bigWigs <code>.cpm.bw</code> and input-normalized bigWigs <code>.input_normalized.log2cpm.bw</code>. Furthermore, it provides an example of how replicates (in bigWig format) were merged.
  * `mass_spec_norm_bw.R` demonstrates how to use mass-spectrometry values to quantitatively scale the depth-normalized bigWigs. 
  * `generate_binned_BED_file.sh` demonstrates how to generate 1kb and 10kb binned <code>BED</code> files from aligned <code>BAM</code> files.
* Figures generated using the scripts above will largely resemble the figures included in the manuscript. However, there may be minor discrepancies as a result of stochasticity in specific algorithms as well as additional manual aesthetic adjustments. Nevertheless, these minor differences do not impact the conclusions presented.
