# Reference-Design

The Broad Institute GATK Best Practices pipeline has helped standardize genomic analysis 
by providing step-by-step recommendations for performing pre-processing and variant 
discovery analysis. 

Intel in collaboration with the Broad Institute put together a guide for the reference 
platform for the GATK Best Practices. The benchmarking eﬀort involves automating the 
GATK Best Practices workﬂow from the Broad Institute as well as providing system-level 
profling data for the germline short variant discovery in whole genomes and exomes. 

The GATK Best Practices Workﬂow is composed of two core pipelines that are to be
performed sequentially: 
- Pre-processing which processes the raw reads to analysisready mapped reads and 
- Variant Discovery which processes the analysis-ready reads to variants. 

The per-sample data pre-processing and variant calling segment of the workﬂow, from 
BWA to GATK Haplotype Caller is implemented as the Single-Sample Calling pipeline and 
the workﬂow steps from GenotypeGVCFs to ApplyRecalibration which operates on a cohort 
of datasets is implemented as the Joint Analysis pipeline.

The ﬂow of both these pipelines are constructed as Perl scripts and is available 
[here](https://github.com/Intel-HLS/Reference-Design/tree/master/Broad_Ref_Pip/scripts). 
In addition to running these tools, system-level data for the individual tools as well as 
the overall pipeline is gathered via [Workﬂow Profiler](https://01.org/workﬂowprofler), 
an Intel Open Source Project. 

Benchmarking results for both the Whole Genome Sequences (WGS) and Exome Sequences (Ex) 
have been collected and presented in the WIKI for the end-to-end workﬂow. 
