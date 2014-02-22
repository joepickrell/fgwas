fgwas
=====

Functional genomics and genome-wide association studies

fwgas is a command line tool for integrating functional genomic information into a genome-wide association study (GWAS). The basic setup is as follows: you have performed a GWAS or a meta-analysis of many GWAS, and have identified tens of loci that influence the disease or trait (our approach works best if there are at least ~20 independent loci with p-values less than 5e-8). We set out to address the following questions:

1. Are these associations enriched in particular types of regions of the genome (coding exons, DNAse hypersensitive sites, etc.)?
2. Can we use these enrichments (if they exist) to identify novel loci influencing the trait?
The underlying model is described in: Pickrell JK (2013) Joint analysis of functional genomic data and genome-wide association studies of 18 human traits. bioRxiv 10.1101/000752

##Installation##
The most up-to-date release is: version 0.2

A user guide is available here: fgwas v0.2 User Guide


##What's new?##

###2/22/14###
We've moved to github

###1/21/14###
Release of version 0.2. Includes support for conditional analysis of annotations using the -cond flag.

###11/8/13###
First public release, version 0.1
