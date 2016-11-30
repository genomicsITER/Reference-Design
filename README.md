# Reference-Design

The Broad Institute GATK Best Practices pipeline has helped standardize genomic analysis 
by providing step-by-step recommendations for performing pre-processing and variant 
discovery analysis. 
Pre-processing refers to generating analysis-ready mapped reads from raw reads using tools 
like BWA*, Picard* tools, and the Genome Analysis Tool Kit. These analysis-ready reads are 
passed through the Variant Calling step of Variant Discovery analysis to generate variants 
per-sample. 
The first part of the GATK Best Practices pipeline takes two FASTQ files, a reference genome, 
and dbSNP and 1000g_indels VCF files as input and outputs a gVCF file per-sample. These gVCF 
files are then further analyzed using Joint Genotyping and Variant Filtering steps of the 
Variant Discovery analysis.

The tools mentioned in the GATK Best Practices Pipeline require enormous computational power 
and long periods of time to complete. Benchmarking such a pipeline allows users to better 
determine the recommended hardware and optimize parameters to help reduce execution time. 
In an effort to advance the standardization and optimization of genomic pipelines, Intel has 
benchmarked the GATK Best Practices pipeline using Workflow Profiler, an open-source tool 
that provides insight into system resources (such as CPU/Disk Utilization, Committed 
Memory, etc.) and helps eliminate resource bottlenecks.
