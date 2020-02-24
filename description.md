# About this app

***Authors: Chelsea Herdman and Bradley Demarest***

The Gene Expression Explorer app allows users to visualize differential expression results
from RNA-seq datasets.  

This Shiny app is supported by the B2B Consortium Grant ([http://www.benchtobassinet.com](http://www.benchtobassinet.com))
and hosted on http://b2b.hci.utah.edu/shiny/zebrafish_heart_timecourse/
***
## Background and Methods

### Experimental Methods

Zebrafish fish hearts were mechanically separated from embyronic zebrafish at
five ages: 24, 36, 48, 60 and 72 hours post-fertilization (hpf). The procedure was 
repeated 4 times for 24, 36 and 48hpf timepoints and 5 times for 60 and 72hpf timepoints 
for a total of 22 biological samples.  
RNA was purified from approximately n = 100 hearts per replicate per timepoint and the 
sample was split equally to produce a total RNA library (Illumina RiboZero Gold Library 
Kit) and a small RNA library (QIAseq miRNA Library Kit). All libraries for both RNA classes
were multiplexed into two single RNA-Seq libraries and sequenced on an Illumina HiSeq 2500.
The Total RNA sequencing library was sequenced paired-end over 8 lanes and the Small RNA 
sequencing library was sequenced single-end over two lanes.

Data is publicly available at https://b2b.hci.utah.edu/gnomex/gnomexFlex.jsp?requestNumber=468R

### Computational Methods

Transcript abundances were quantified using kallisto (*Bray et al. 2016*) and genome build 
GRCz10 release 89 (may2017.archive.ensembl.org). Estimated counts for all transcripts 
per gene were summed to give a gene-level abundance estimation (*Soneson et al. 2015*).  
Summed estimated counts were rounded to the nearest integer in order to run DESeq2 (*Love et al. 2014*) using a negative binomial LRT model correcting for replicate 
(counts~ replicate + timepoint). This model tests in an anova-like way whether one timepoint
differs from any other timepoint across the series.

### Visualization Features

Gene panel plots can be displayed using normalized counts or rlog values from DESeq2 or 
transcripts per million from kallisto. The light grey dot represents the mean expression 
value for that gene at each timepoint and the dark grey dots represent the value for each 
replicate. The vertical line at each timepoint depicts the range of the data and a line has 
been drawn to connect the mean at each timepoint to show the expression profile.  

The heatmap displays z-scores (computed using DESeq2 normalized counts) for the selected gene list 
using a red-blue color scale.

***
#### References
* Bray N, Pimentel H, Melsted P, Pachter L (2016), _Near-optimal probabilistic RNA-seq quantification_, Nature Biotechnology, 34, 525–527. doi:10.1038/nbt.3519  
* Soneson C, Love MI, Robinson MD (2015), _Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences_, F1000Research, 4, 1521. doi:10.12688/f1000research.7563.2  
* Love MI, Huber W, Anders S (2014), _Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2_, Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8  

***
*Please email Chelsea Herdman <cherdman@genetics.utah.edu> with any questions or issues*
