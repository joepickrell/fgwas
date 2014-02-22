fgwas
=====

Functional genomics and genome-wide association studies

fgwas is a command line tool for integrating functional genomic information into a genome-wide association study (GWAS). The basic setup is as follows: you have performed a GWAS or a meta-analysis of many GWAS, and have identified tens of loci that influence the disease or trait (our approach works best if there are at least ~20 independent loci with p-values less than 5e-8). We set out to address the following questions:

1. Are these associations enriched in particular types of regions of the genome (coding exons, DNAse hypersensitive sites, etc.)?
2. Can we use these enrichments (if they exist) to identify novel loci influencing the trait?

The underlying model is described in: Pickrell JK (2014) [Joint analysis of functional genomic data and genome-wide association studies of 18 human traits](http://biorxiv.org/content/early/2014/01/22/000752). bioRxiv 10.1101/000752

##Installation##

###Dependencies###
fgwas depends on:

-the [GNU Scientific Library](http://www.gnu.org/software/gsl/)

-the [Boost Libraries](http://www.boost.org)

###Quick Start###
The most up-to-date release is: version 0.3. See ["Releases"](https://github.com/joepickrell/fgwas/releases) above.
After downloading fgwas-0.3.tar.gz at the link above, run:

>tar -xvf fgwas-0.3.tar.gz

>cd fgwas-0.3

>./configure

>make

A user guide is available here: [fgwas v0.2 User Guide](https://github.com/joepickrell/fgwas/blob/master/man/fgwas_manual.pdf)

Previous versions are available from the [Google Code repository](https://code.google.com/p/gwas/).

##What's new?##

###2/22/14###
Release of version 0.3. Minor updates including sanity checks on input files.
We've moved to github!

###1/21/14###
Release of version 0.2. Includes support for conditional analysis of annotations using the -cond flag.

###11/8/13###
First public release, version 0.1
