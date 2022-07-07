# TSS extractR

[![GitHub issues](https://img.shields.io/github/issues/daniilidisK/TSS_extractR?color=green)](https://github.com/daniilidisK/TSS_extractR/issues/new)
[![License:MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FdaniilidisK%2FTSSextractR&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)

<!--
[![CRAN Version](https://www.r-pkg.org/badges/version/CMplot?color=yellow)](https://CRAN.R-project.org/package=CMplot) 

[![](https://img.shields.io/badge/GitHub-4.1.0-blueviolet.svg)]()

![](http://cranlogs.r-pkg.org/badges/grand-total/CMplot?color=red) 

[![](https://cranlogs.r-pkg.org/badges/last-month/CMplot)](https://CRAN.R-project.org/package=CMplot) 
-->


## An automated method to identify Transcription Start Sites (TSS) in microbiome genomes, using Cappable Sequencing

[Cappable-seq](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2539-z) is a high throughput Next Generation Sequencing (NGS) technique, which emphasizes the enrichment of the 5' ends of the transcripts. This method is based on the biochemical identification of Transcription Start Sites (TSS) in prokaryotic organisms. These specific regions carry at their 5' end a triphosphate group, also known as a 'cap' in prokaryotic organisms, because it is the first nucleoside, during chain polymerization. Using Vaccinia Capping Enzyme (VCE), the triphosphorylated 5' end can be capped with the eukaryotic capping template, but instead of Guanosine Triphosphates (GTPs), are using destheiobiotinylated Guanosine Triphosphates, which are simple GTPs, fused with a biotin derivative, called destheiobiotin. So, this destheiobiotinylated cap can be captured with Streptavidin Magnetic Beads and in this way, the intact 5' RNA fragments can be isolated and sequenced.

TSSextractR is an R package for transcriptome analysis. This package can detect and visualize candidate TSS from Cappable-seq Sequencing data, analyze them statistically and find conserved motifs of the promoters and Shine - Dalgarno regions. Also, this package can annotate TSS into the categories: Upstream, Internal, Sense, AntiSense and Representative. Finally, it can detect and annotate leaderless *m*RNAs, *cis* and *trans* acting small RNAs and operons.


## Installation
TSSextractR can be installed with the following R code:
```r
if(!require("devtools")) {
    install.packages("devtools")
} else {
    install_github("daniilidisK/TSSextractR")
}
```

## Contact
Bug reports, questions and suggestions are welcome and appreciated.
- **Author:** Daniilidis Konstantinos
- **Contact:** daniilidis@ieee.org
- **Affiliations:** [*Computer Science and Biomedical Informatics, University of Thessaly*](http://dib.uth.gr/) and [*Diana Lab*](http://www.dianalab.gr/)
